#pragma once
// Consolidation phase strategy (Stage B9). Time-dependent Biot consolidation,
// PLAXIS "Consolidation": the configuration's load increment dF = f - B (the
// SumMstage imbalance) is applied at t = 0+, generating an undrained excess
// pore pressure that dissipates over the phase's time interval -- the classic
// settlement-time (Terzaghi U-t) development. v1: soil-only.
//
// Same contract as the transient-flow pilot: neutral inputs resolved at the
// caller's seam, refusal messages engine-owned and byte-identical, material
// names carried as diagnostic labels. Two differences worth stating:
//   - the linear solver enters as factories built by the composition root
//     (one factorization, back-substitution every step), never named here;
//   - on success this strategy does NOT set R.ok or R.message -- unlike the
//     flow-only pilot, a consolidation phase falls through to the driver's
//     common result tail, exactly as it always has.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/constants.hpp>
#include <katai/analysis/hydraulic_boundary.hpp>
#include <katai/analysis/post/consolidation_recovery.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/analysis/results.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// One material's consolidation description, resolved at the caller's seam.
// Constitutive behaviour (incl. the Soft Soil Creep detection) comes from the
// engine's own `models` table, not from here.
struct ConsolidationPhaseMaterial {
    std::string name;              // diagnostic label for honest refusals
    double kx = 1.0, ky = 1.0;     // permeability [m/day]
    double porosity = 0.3;         // n = e / (1 + e)
    bool nonporous = false;        // refused: a non-porous region would silently behave water-filled
    // Undrained (C) is a TOTAL stress material: it has no pore pressures to consolidate.
    // PLAXIS states the same fact as "a Consolidation calculation does not affect Undrained (C)
    // materials"; this build refuses the phase rather than solving part of the mesh in total
    // stress and the rest in effective stress with one Kw/n between them.
    bool total_stress = false;
};

// The phase's neutral configuration.
struct ConsolidationPhase {
    std::vector<ConsolidationPhaseMaterial> materials;  // by material id
    std::vector<FlowEdge> flow_edges;                   // B4 vocabulary, Closed edges included
    bool have_flow_bcs = false;                         // unfiltered declaration scan (B4 rule)
    std::vector<char> active;                           // element activity; empty = everything active
    // Nodes of the DRAINS active in this phase. "In consolidation analysis, drains reduce the
    // excess pore pressure to zero and the specified head is ignored" (PLAXIS Ref sec. 5.9.2) --
    // which is exactly what this solver's drained-node mask means, so a drain enters here rather
    // than as a new kind of boundary condition.
    std::vector<int> drain_nodes;
    bool has_structural_elements = false;               // plates/anchors/geogrids/walls/embedded present
    double duration_day = 1.0;                          // <= 0 falls back to 1 day
    int time_steps = 25;                                // clamped to [1, 2000]
    double yscale = 1.0;                                // model height, for the top-drain tolerance
};

// Solve the phase. Fills the settlement-time series, the final displacement and
// the recovered nodal effective stresses in R, and returns the final committed
// Gauss state in `committed_out` (the phase chain carries it). Returns false on
// an honest refusal or non-convergence, with R.message set; on success the
// caller's common result tail completes R.
inline bool solve_consolidation_phase(
    const katai::mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<MaterialProfile>& profiles,
    const std::vector<GaussState>& init, const Eigen::VectorXd& dF, bool nonlinear_soil,
    const ConsolidationPhase& in,
    const ConsolidationSolveFactory& le_factory, const ConsolidationSolveFactory& plastic_factory,
    SolveResult& R, std::vector<GaussState>& committed_out) {
    if (in.has_structural_elements) {
        R.message = "Consolidation v1 is soil-only: structural elements (plates / anchors / "
                    "geogrids / walls / embedded beams) are not supported in a consolidation "
                    "phase yet -- use a Plastic phase for them.";
        return false;
    }
    const size_t nmat = in.materials.size();
    // Soft Soil Creep x consolidation: honest refusal. End-to-end verification WAS attempted and
    // failed (2026-07-21, test_soft_soil_gui): on a fast-draining column the settlement DEPENDS on
    // the time-step size and the creep tail came out 37% short of the exact mu* H ln(100) oracle
    // that the verified Plastic-phase path holds to +0.1%. Until the Biot x creep interaction is
    // resolved: apply loads in a timed Plastic phase, model holding as a timed no-change Plastic
    // phase (both pinned against the isotach closed forms).
    for (size_t mi = 0; mi < nmat; ++mi) {
        bool used_m = false;
        for (int e = 0; e < mesh.element_count && !used_m; ++e)
            if ((in.active.empty() || in.active[e]) && mesh.element_material[e] == (int)mi) used_m = true;
        if (used_m && models[mi].type == MaterialType::SoftSoilCreep) {
            R.message = "Material '" + in.materials[mi].name + "' is Soft Soil Creep: the "
                        "consolidation (Biot) x creep interaction is not verified yet -- "
                        "measured settlements DEPEND on the time-step size on this path, so "
                        "it is refused rather than shipped silently wrong. Apply loads in a "
                        "Plastic phase WITH its real time interval, and model holding "
                        "periods as a chained no-change Plastic phase with the holding time "
                        "(both verified against the isotach closed forms).";
            return false;
        }
    }
    // Non-porous x consolidation: honest refusal. The Biot solver gives every active element a
    // pore-pressure DOF and ONE scalar Kw/n -- a non-porous region would silently behave
    // water-filled (concrete gaining pore stiffness, drainage paths corrupted). Region-based
    // pore-DOF exclusion is a separate work item (with per-material Kw/n).
    for (size_t mi = 0; mi < nmat; ++mi) {
        bool used_np = false;
        for (int e = 0; e < mesh.element_count && !used_np; ++e)
            if ((in.active.empty() || in.active[e]) && mesh.element_material[e] == (int)mi) used_np = true;
        if (used_np && in.materials[mi].nonporous) {
            R.message = "Material '" + in.materials[mi].name + "' is Non-porous: consolidation "
                        "/ fully-coupled phases give every element a pore-pressure DOF and a "
                        "single fluid stiffness, so a non-porous region would silently behave "
                        "water-filled. Model the concrete with Drained + high stiffness in "
                        "consolidation phases, or keep Non-porous to static/dynamic phases.";
            return false;
        }
        if (used_np && in.materials[mi].total_stress) {
            R.message = "Material '" + in.materials[mi].name + "' is Undrained (C), a total "
                        "stress analysis: it carries no pore pressure, so there is nothing in it "
                        "to consolidate (PLAXIS: \"a Consolidation calculation does not affect "
                        "Undrained (C) materials\"). Solving the phase would put part of the mesh "
                        "in total stress and the rest in effective stress. Give the material "
                        "effective parameters with Drained or Undrained (A)/(B) for the "
                        "consolidating phases.";
            return false;
        }
    }
    // Permeability per material [m/day]; every used material needs k > 0.
    std::vector<Permeability> cperm(nmat, {1.0, 1.0});
    std::vector<char> used(nmat, 0);
    for (int e = 0; e < mesh.element_count; ++e)
        if (in.active.empty() || in.active[e]) used[mesh.element_material[e]] = 1;
    for (size_t mi = 0; mi < nmat; ++mi) {
        cperm[mi] = {in.materials[mi].kx, in.materials[mi].ky};
        if (used[mi] && (cperm[mi].kx <= 0.0 || cperm[mi].ky <= 0.0)) {
            R.message = "Material '" + in.materials[mi].name + "' has no permeability -- set "
                        "kx, ky (> 0) on its Groundwater tab before a consolidation phase.";
            return false;
        }
    }
    // Pore-fluid stiffness Kw/n from the real water bulk modulus (Verruijt; PLAXIS Sci.Man sec. 4):
    // near-incompressible -> cv = k Eoed / gamma_w. v1 uses one representative porosity.
    constexpr double kWaterBulk = 2.0e6;   // bulk modulus of water [kPa]
    double porosity = 0.3;
    for (size_t mi = 0; mi < nmat; ++mi)
        if (used[mi]) { porosity = in.materials[mi].porosity; break; }
    const double kw_over_n = kWaterBulk / std::max(0.05, porosity);

    // Drainage boundary (engine service, B4): prescribed-head / seepage edges drain; with no
    // declared flow BCs the model top drains; inactive-only nodes carry no pore DOF.
    // (not const: an active drain adds its own nodes to the same mask -- the drain and the
    // draining boundary say the same thing to this solver, that the excess pore pressure there
    // is zero.)
    std::vector<char> drained =
        flow_drained_nodes(in.flow_edges, in.have_flow_bcs, mesh, in.active, in.yscale);
    for (int n : in.drain_nodes)
        if (n >= 0 && n < mesh.node_count) drained[n] = 1;

    // Time stepping from the phase duration / step count.
    const double duration = in.duration_day > 0.0 ? in.duration_day : 1.0;
    const int nsteps = std::clamp(in.time_steps, 1, 2000);
    const double dt = duration / nsteps;
    const std::vector<double> p0(mesh.node_count, 0.0);   // load generates the excess pore at t=0+

    // Solve with the linear-elastic OR the elastoplastic (MC/HS) consolidation core. The
    // series + final committed EFFECTIVE Gauss state come back through a common pointer.
    const ConsolidationResult* series = nullptr;
    std::vector<GaussState> committed;
    ConsolidationResult le_res;
    ConsolidationPlasticResult pl_res;
    if (nonlinear_soil) {
        // Elastoplastic skeleton: effective stress from the constitutive return mapping; each
        // time step is a monolithic coupled Newton. K_T is nonsymmetric for non-associated
        // flow -> the caller's plastic factory selects the matching solver. The previous
        // phase's committed effective stresses (`init`) seed the state; dF drives it.
        pl_res = solve_consolidation_plastic(mesh, dofs, models, cperm, kGammaWater,
            kw_over_n, drained, init, {}, dt, nsteps, in.active, &dF, plastic_factory,
            40, 1e-6, profiles);
        if (!pl_res.converged) {
            R.message = "Elastoplastic consolidation did not converge in a time step (the load "
                        "increment may exceed the soil capacity, or the time step is too large "
                        "-- use a smaller load increment or more, smaller steps).";
            return false;
        }
        series = &pl_res.series;
        committed = std::move(pl_res.committed);   // final effective Gauss state (from return mapping)
    } else {
        le_res = solve_consolidation(mesh, dofs, models, cperm, kGammaWater,
                                     kw_over_n, drained, p0, dt, nsteps, in.active, &dF,
                                     le_factory, profiles);
        if (le_res.displacement.empty()) { R.message = "Consolidation produced no result."; return false; }
        series = &le_res;
        committed = init;   // sigma' = init + D B v_final (linear-elastic skeleton)
        recover_consolidation_stress(mesh, dofs, models, le_res.displacement.back(), in.active, committed);
    }

    // Settlement-time curve: max vertical settlement |u_y| + max |excess pore| per step.
    for (size_t s = 0; s < series->times.size(); ++s) {
        double smax = 0.0, pmax = 0.0;
        for (int n = 0; n < mesh.node_count; ++n) {
            smax = std::fmax(smax, std::fabs(series->displacement[s][dofs.global_dof(n, 1)]));
            pmax = std::fmax(pmax, std::fabs(series->pore[s][n]));
        }
        R.consol_time.push_back(series->times[s]);
        R.consol_settlement.push_back(smax);
        R.consol_excess_pore.push_back(pmax);
    }
    R.disp = series->displacement.back().head(mesh.node_count * 2);
    R.stress = recover_nodal_stresses_from_gauss(mesh, committed, in.active);
    R.load_factor = 1.0;
    R.iterations = nsteps;
    committed_out = std::move(committed);
    return true;
}

} // namespace katai::core
