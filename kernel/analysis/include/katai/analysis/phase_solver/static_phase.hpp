#pragma once
// Static phase strategy (Stage B9): the K0 / gravity-loading / staged-Plastic
// solve -- one Newton ramp about an internal-force baseline, then the result
// assembly (committed stress recovery, structural diagrams, chain carrier).
//
// Ramp semantics (PLAXIS parity). K0 procedure: the geostatic state (seeded
// soil stress + interface sigma_n0) IS the equilibrium; its internal force is
// held as the constant baseline B and ONLY external loads ramp, so
// residual(0) = B - f_int(0) = 0 on any mesh and self-weight is NOT
// double-counted during the ramp (essential once the soil is plastic --
// ramping gravity would yield spuriously). Gravity loading: B = 0, ramp the
// full body force + loads from a stress-free start. PLAXIS K0 caveat
// (Reference Manual): the column-overburden K0 field is in equilibrium ONLY
// when the ground surface, layer boundaries and the water table are all
// horizontal. On a slope the seeded field leaves a GENUINE out-of-balance
// force d = f_body - f_int(seed) (principal directions must rotate; shear must
// develop near the face). PLAXIS resolves it with a plastic nil-step; here
// that step is built into the K0 phase by ramping the imbalance together with
// the external loads, so target(1) = f exactly (true equilibrium, comparable
// to gravity loading). It is triggered by GEOMETRY (the caller's nil_step
// flag), not by |d|: on level geometry d is only the quadrature residual and
// the seed IS the intended answer -- keeping the undisturbed-case identity
// u = 0 exactly. Flow coupling always takes the nil-step (the hydrostatic K0
// seed cannot represent a non-hydrostatic seepage pore field), and a CHAINED
// staged phase is the same rule writ large: ramp = f(active) - B(committed)
// is exactly the PLAXIS SumMstage (excavation unloading / fill weight / new
// loads).
//
// Same contract as the other phase strategies: neutral inputs resolved at the
// caller's seam, refusal/collapse messages engine-owned and byte-identical,
// the linear solver injected by the composition root. Consolidation-class
// return discipline: on success this strategy does NOT set R.ok or R.message
// -- the phase falls through to the driver's common result tail.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>     // solve_nonlinear, Structures, StructuralInit
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/analysis/results.hpp>
#include <katai/analysis/structural_diagrams.hpp>  // DiagSpec/IfaceDiag + per-line force_diagram
#include <katai/analysis/structural_forces.hpp>    // force_envelope
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// The phase's neutral configuration.
struct StaticPhase {
    bool baseline = false;        // ramp external loads about the internal-force baseline B
    bool nil_step = false;        // add the configuration imbalance f - f_loads - B to the ramp
    bool axisymmetric = false;    // r-z kinematics
    int load_steps = 1;           // Newton load stepping (grows with nonlinearity)
    double tolerance = 1e-10;     // Newton tolerated error
    // Maximum Newton iterations per load increment; 0 = this strategy's own 80. A step that
    // needs more than this is not solved: the solver cuts the increment back and tries again,
    // so the number is a patience setting, not an accuracy one (PLAXIS "Max iterations", 60).
    int max_iterations = 0;
    double time_interval_day = 0.0;  // SSC creep time of a chained Plastic phase [days]
    std::vector<char> active;     // element activity; empty = everything active
    // Prescribed displacements (nonzero Dirichlet): full-DOF vector, ramped 0 -> value
    // together with the load (the caller fixed the corresponding dofs). Empty = none.
    Eigen::VectorXd presc;
};

// Solve the phase. `carry_full` is the parent's converged full-DOF displacement
// datum (null = no carriable parent; carry_init is then empty and the
// structures start from zero, the pre-Track-1a behaviour). On collapse
// (non-convergence) returns false with the honest limit-load message in
// R.message; on success fills R (nil_step, load factor, timings, displacement,
// recovered stress, structural diagrams, interface results, chain carrier) and
// writes the committed Gauss states to `out_states`, then returns true for the
// caller's common result tail to complete R.
inline bool solve_static_phase(
    const mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<MaterialProfile>& profiles,
    const std::vector<GaussState>& init,
    const Eigen::VectorXd& f, const Eigen::VectorXd& f_loads, const Eigen::VectorXd& B,
    const LinearSolve& solver, const Structures& structures,
    const std::vector<DiagSpec>& diag_specs, const std::vector<IfaceDiag>& iface_diags,
    const StructuralInit& carry_init, const Eigen::VectorXd* carry_full,
    const StaticPhase& in, SolveResult& R, std::vector<GaussState>* out_states) {
    const std::vector<char>& act = in.active;
    const bool static_carry = carry_full != nullptr;
    Eigen::VectorXd ramp = in.baseline ? f_loads : f;
    if (in.nil_step) ramp += f - f_loads - B;   // ramp the configuration imbalance
    NewtonOptions nopt{in.load_steps, in.max_iterations > 0 ? in.max_iterations : 80,
                       in.tolerance};
    nopt.kinematics = in.axisymmetric ? Kinematics::Axisymmetric
                                      : Kinematics::PlaneStrain;
    // SSC TIME (Stage 3): a staged Plastic phase's Time interval [days] enters the constitutive
    // path -- the solver apportions it over the increments by the delta-lambda ratio (the Stage-2
    // contract, SumMstage parity). Only SoftSoilCreep reads it; every other model ignores it, so
    // projects without SSC are bit-identical. The caller passes a nonzero value ONLY for chained
    // (staged) phases: in PLAXIS the initial phase is TIMELESS -- an unconditional duration here
    // once leaked 1 day of creep into the K0 phase and broke the geostatic identity (measured:
    // K0 max |u| = 26 mm = exactly mu* ln2 H on an SSC column).
    nopt.time_interval = in.time_interval_day;
    // The nil-step imbalance is plastic redistribution work: give the solver load steps to do it
    // gradually even when the soil model alone would not have asked for stepping.
    if (in.nil_step) nopt.load_steps = std::max(nopt.load_steps, 10);
    // Structural init only on the carry path: passing a datum whose f_s0 is NOT in B would leave
    // residual(0) != 0 and ramp a spurious imbalance (the caller's B augmentation and this
    // argument must switch together).
    const StructuralInit no_init;
    const NewtonResult nr = solve_nonlinear(
        mesh, dofs, models, ramp, solver, nopt, init, act, structures, in.presc, B, profiles,
        static_carry ? carry_init : no_init);
    R.nil_step = in.nil_step;
    R.load_factor = nr.load_factor;
    R.timings = nr.timings;
    R.iterations = nr.total_iterations;
    if (!nr.converged) {
        R.message = "Did not fully converge: equilibrated " +
                    std::to_string((int)std::lround(100.0 * nr.load_factor)) +
                    "% of the applied load before a collapse mechanism formed (the remaining load "
                    "exceeds the soil capacity). This equilibrated fraction is the incremental "
                    "limit (collapse) load.";
        return false;
    }
    if (out_states) *out_states = nr.gauss_states;   // committed -> next phase
    // TOTAL displacement for the structural elements: parent datum + this phase's increment
    // (identical to what the solver evaluated them at). Without carry it is just this phase's
    // field, as before. R.disp stays the PHASE INCREMENT (PLAXIS convention) -- only the
    // structural diagrams and the chain carrier need the total.
    Eigen::VectorXd disp_total = nr.displacement;
    if (static_carry) disp_total += *carry_full;
    // Track 1a: raw structural state for the next chained phase / nonlinear Dynamic child (see
    // StructCarryState). The full displacement keeps the extra structural DOFs that R.disp
    // truncates away -- the interface's structure side lives exactly there. With carry it is the
    // ACCUMULATED total, so grandchildren continue from the right datum; the plastic vectors are
    // the solver's committed (seeded + evolved = totals either way).
    R.struct_state.full_disp = disp_total;
    R.struct_state.anchor_plastic = nr.anchor_plastic;
    R.struct_state.geogrid_plastic = nr.geogrid_plastic;
    R.struct_state.interface_slip = nr.interface_slip;
    R.struct_state.interface5_slip = nr.interface5_slip;
    R.struct_state.embedded_skin_slip = nr.embedded_skin_slip;
    R.struct_state.embedded_foot_slip = nr.embedded_foot_slip;
    R.struct_state.plate_plastic = nr.plate_plastic;
    R.struct_state.plate5_plastic = nr.plate5_plastic;
    R.disp = nr.displacement.head(mesh.node_count * 2);
    // Recover the EFFECTIVE stress from the solver's committed Gauss states (includes the K0
    // initial field + any increments) -- not re-derived from displacement, which would miss K0.
    // Passive (excavated) elements are excluded from the nodal averaging.
    R.stress = recover_nodal_stresses_from_gauss(mesh, nr.gauss_states, act);
    // Support reactions at the fixed dofs, from the SAME committed states (discrete B^T sigma,
    // no recovery smoothing): what a support -- or an imposed displacement -- must exert. Soil
    // contribution (structural end forces on a fixed node are not included in v1; stated in
    // results.hpp). Free dofs are zeroed so the field reads as "reactions", not "f_int".
    {
        Eigen::VectorXd f_int =
            nodal_internal_force_from_gauss(mesh, nr.gauss_states, in.axisymmetric, act);
        R.reaction = Eigen::VectorXd::Zero(2 * mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n)
            for (int c = 0; c < 2; ++c)
                if (dofs.is_fixed(dofs.global_dof(n, c)))
                    R.reaction[2 * n + c] = f_int[2 * n + c];
    }

    // Structural force diagrams (PLAXIS Output -> Structures): from the converged TOTAL solution
    // (disp_total covers the rotation/independent extra DOFs AND the carried parent datum) + the
    // committed plastic state, via the validated post-processors in structural_forces.hpp --
    // exactly the state the solver applied its caps to (the D6b rule).
    for (const auto& sp : diag_specs) {
        StructForce d = force_diagram(sp, structures, mesh, dofs, disp_total,
                                      nr.anchor_plastic, nr.geogrid_plastic,
                                      /*elastic=*/false, nr.plate_plastic, nr.plate5_plastic);
        // Plate M-N hinge formed (this phase or carried in): committed plastic state nonzero
        // anywhere on this chain -- return mapping only ever writes a nonzero value when it fires.
        if (sp.kind == 0 || sp.kind == 5) {
            const auto& psv = sp.kind == 0 ? nr.plate_plastic : nr.plate5_plastic;
            const size_t S = sp.kind == 0 ? (size_t)plate::kPlasticStateSize
                                          : (size_t)plate::kPlasticStateSize5;
            for (size_t k = sp.begin * S; k < sp.end * S && k < psv.size(); ++k)
                if (psv[k] != 0.0) { d.yielded = true; break; }
        }
        if (static_carry) {
            // Station ux/uy drive the deformed overlay, and the GUI deforms the mesh by the
            // PHASE INCREMENT (R.disp) -- so the overlay must use the increment too, or the
            // structure would visually detach from its soil. Forces stay evaluated at the total.
            const StructForce dinc = force_diagram(sp, structures, mesh, dofs,
                                                   nr.displacement, {}, {});
            for (size_t k = 0; k < d.stations.size() && k < dinc.stations.size(); ++k) {
                d.stations[k].ux = dinc.stations[k].ux;
                d.stations[k].uy = dinc.stations[k].uy;
            }
        }
        const auto env = force_envelope(d.stations);
        d.max_N = env.max_abs_N; d.max_Q = env.max_abs_Q; d.max_M = env.max_abs_M;
        R.struct_forces.push_back(std::move(d));
    }
    // Interface results (PLAXIS Output -> Interfaces): tau / sigma_n / relative slip along each
    // Coulomb joint, recomputed from the converged TOTAL displacement + committed slip (mirrors
    // the solver).
    for (const auto& is : iface_diags) {
        InterfaceResult ir = force_diagram(is, structures, mesh, dofs, disp_total,
                                           nr.interface_slip, nr.interface5_slip, false);
        ir.slip_checked = true;   // the static solve applies the Coulomb return -- a real check
        for (const auto& st : ir.stations) {
            ir.max_abs_tau = std::max(ir.max_abs_tau, std::fabs(st.tau));
            ir.max_abs_sigma_n = std::max(ir.max_abs_sigma_n, std::fabs(st.sigma_n));
            ir.max_abs_slip = std::max(ir.max_abs_slip, std::fabs(st.slip));
            if (st.slipping) ir.any_slip = true;
        }
        R.interface_forces.push_back(std::move(ir));
    }
    return true;
}

} // namespace katai::core
