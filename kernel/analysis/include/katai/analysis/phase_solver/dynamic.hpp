#pragma once
// Dynamic (seismic) phase strategy (Stage B9). PLAXIS "Dynamic":
// M u'' + C u' + K u = -M r a_g(t), Newmark-beta (gamma=1/2, beta=1/4) +
// Rayleigh damping (dynamics.hpp; docs/references/dynamic-seismic-formulation.md).
// v1: linear-elastic skeleton, HORIZONTAL base acceleration (site response /
// SSI), rigid base (the user fixes the base) or the compliant (absorbing) base,
// free or Lysmer free-field lateral sides, opt-in fully nonlinear time
// integration (per-step Newton), seismic structural envelopes, response
// spectrum + TBDY/EC8 design-spectrum overlays.
//
// Same contract as the other phase strategies: neutral inputs resolved at the
// caller's seam (per-material unit weights arrive drainage-resolved; the
// model-family stiffness selection reads the MOTOR's models table, never a
// schema enum), refusal messages engine-owned and byte-identical, and every
// factorization enters as a callback built by the composition root. Return
// discipline is TransientFlow-class: this strategy fills R.ok and R.message
// itself and never falls through to the driver's common result tail.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/constants.hpp>            // kGravity, kUtilCap
#include <katai/analysis/dynamics.hpp>             // solve_newmark, DynamicsSolveFactory, rayleigh_from_modes
#include <katai/analysis/dynamics_nonlinear.hpp>   // solve_newmark_nonlinear (opt-in per-step Newton)
#include <katai/analysis/free_field.hpp>           // 1D free-field column (lateral boundary, D3b)
#include <katai/analysis/nonlinear_solver.hpp>     // Structures, StructuralInit, LinearSolve
#include <katai/analysis/post/stress_recovery.hpp> // nodal recovery of the committed Gauss state
#include <katai/analysis/response_spectrum.hpp>    // response spectrum + TBDY/EC8 design spectra
#include <katai/analysis/results.hpp>              // SolveResult, StructForce, InterfaceResult
#include <katai/analysis/structural_diagrams.hpp>  // DiagSpec/IfaceDiag + per-line force_diagram
#include <katai/analysis/structural_dynamics.hpp>  // structural K/M, seismic_influence_x
#include <katai/analysis/structural_forces.hpp>    // force_envelope
#include <katai/fem/assembly/assembler.hpp>        // mass/stiffness/dashpot assembly, expand_to_full
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/element_traits.hpp>   // Tri6Element/Tri15Element (Gauss counts)
#include <katai/materials/linear_elastic.hpp>
#include <katai/materials/material_model.hpp>      // MaterialModel/Type, hs frozen stiffness, profiles
#include <katai/materials/soft_soil.hpp>           // softsoil::kPmin
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/boundary_extraction.hpp>      // extract_boundary_edges, BoundaryEdgeChain
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Engine-owned seismic wave kinds. Mirrors the schema's enumerators with
// file-stable values (the B1 DesignApproach precedent); the schema enum is
// mapped ONCE, at the driver seam.
enum class SeismicWave { Harmonic = 0, Ricker = 1, Record = 2 };

// The phase's neutral configuration. Field defaults are the values the driver
// used when no phase configuration was given, so a null-config caller builds
// this struct and gets the identical run.
struct DynamicPhase {
    // -- time stepping
    double duration = 1.0;                        // shaking duration [s]
    int time_steps = 200;
    // -- seismic input
    SeismicWave wave = SeismicWave::Harmonic;
    double amplitude = 1.0;                       // a_g scale [m/s^2] (Record: multiplies the samples)
    double frequency = 2.0;                       // harmonic/Ricker frequency [Hz]
    std::vector<double> accel_record;             // uniformly sampled a_g [m/s^2] (Record)
    double record_dt = 0.02;                      // record sampling step [s]
    double rayleigh_f1 = 1.0, rayleigh_f2 = 8.0;  // Rayleigh target band [Hz]
    double damping_ratio = 0.05;                  // xi at both targets
    bool free_field = false;                      // Lysmer free-field lateral boundaries (D3b)
    bool compliant_base = false;                  // absorbing base (must match the BC setup's flag)
    bool nonlinear = false;                       // per-step-Newton time integration (opt-in)
    // -- code-spectrum overlays
    double tbdy_ss = 1.0, tbdy_s1 = 0.4;          // TBDY 2018 map coefficients
    int site_class = 2;                           // ZA..ZE as 0..4
    bool ec8_enabled = false;
    double ec8_gamma = 1.0, ec8_agr = 0.0;        // gamma_I, a_gR [g]
    int ec8_ground = 0, ec8_type = 0;             // ground type A..E as 0..4; 0 = Type 1
    // -- model
    bool axisymmetric = false;                    // r-z mode; refused (not available yet)
    std::vector<char> active;                     // element activity; empty = everything active
    bool has_water = false;                       // a water table exists (saturation checks)
    bool phreatic_mass = false;                   // declared table with >= 2 points (mass assembly)
    std::function<double(double)> water_table_y;  // y_w(x); consulted only under the flags above
    std::vector<double> gamma;                    // per-material unit weight above the table [kN/m^3]
    std::vector<double> gamma_sat;                // per-material saturated unit weight, drainage-
                                                  // resolved at the seam (NonPorous holds no water:
                                                  // its saturated weight IS gamma_unsat)
};

// Solve the phase. `carry_full` is the parent's converged full-DOF displacement
// (null when no carriable parent; carry_init is then empty), `parent` the parent
// phase's result for the static superposition (may be null), `init_states` the
// parent's committed Gauss states (stress-dependent frozen stiffness + the
// nonlinear path's continuation). `spd_factory` factors an SPD system once and
// back-solves per call; `make_nonsym_solver` builds the nonsymmetric Newton
// solver and is invoked only when the phase asks for the nonlinear path.
// Returns false on an honest refusal with R.message set; on success R is
// COMPLETE (R.ok, R.message, fields) except R.mesh, which the caller moves in.
inline bool solve_dynamic_phase(
    const mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<LinearElastic>& mats,
    const std::vector<MaterialProfile>& profiles, const Structures& structures,
    const std::vector<DiagSpec>& diag_specs, const std::vector<IfaceDiag>& iface_diags,
    const StructuralInit& carry_init, const Eigen::VectorXd* carry_full,
    const SolveResult* parent, const std::vector<GaussState>* init_states,
    const DynamicPhase& in, const DynamicsSolveFactory& spd_factory,
    const std::function<LinearSolve()>& make_nonsym_solver, SolveResult& R) {
    if (in.axisymmetric) {
        R.message = "Dynamic (seismic) analysis is not available in axisymmetric mode yet.";
        return false;
    }
    const std::vector<char>& act = in.active;
    // NONLINEAR opt-in: route the time integration to the per-step-Newton solver so the soil can
    // plastify and the structures can slip/yield DURING shaking. Default off = the LINEAR path
    // below, bit-identical to before. When on, the interface post-processor must use the SAME
    // (Coulomb) constitutive the solver applied -- an elastic tau would report a stress the solver
    // never produced (the D6b lesson) -- so iface_elastic tracks it.
    const bool dyn_nl = in.nonlinear;
    const bool iface_elastic = !dyn_nl;
    // COMPLIANT BASE combinations (both former v1 refusals are now IMPLEMENTED):
    // + free-field sides: the 1D side columns are solved on a compliant base too
    //   (solve_free_field_column_compliant -- same rho Vs halfspace, same half-of-within
    //   drive), so sides and base speak the same total-motion physics.
    // + nonlinear dynamic: the solver's rest-start baseline argument holds verbatim in the
    //   total-motion frame (F(t) carries no static part; committed sigma balances it;
    //   v_up(0) = 0 -> F(0) = 0 -> exact rest start). V&V: test_compliant_base (e)-(g).
    const bool compliant_base = in.compliant_base;
    // SOIL-STRUCTURE INTERACTION: plates/walls, anchors, geogrids and interfaces take part in the
    // dynamic system -- their ELASTIC stiffness goes into K and the plate inertia into M
    // (analysis/structural_dynamics.hpp), so the wall and the soil are solved together. Every
    // structural DOF the DofMap holds must be assembled: a rot/trans DOF left with an all-zero
    // row makes K_eff singular and the linear solver corrupts the heap (an ACCESS VIOLATION, not
    // a C++ exception -- /EHsc catch(...) does NOT catch it). Embedded beams (pile rows) are in
    // too now: their skin/foot coupling is mesh-NONCONFORMING, so it needs the soil element's
    // shape functions -> assemble_structural_stiffness dispatches on mesh.nodes_per_element
    // (the assembler.cpp pattern) instead of refusing them.
    // HARDENING SOIL / HSsmall in a Dynamic phase. The dynamic system is linear-elastic and takes
    // its stiffness from `mats` (E', nu) -- but the Mechanical tab never exposes E' for an HS
    // material (it offers E50_ref/Eoed_ref/Eur_ref), so E' would keep its struct default and the
    // ENTIRE seismic run would ride a stiffness the user never entered. The right answer is not to
    // refuse: an HS material DOES carry the information, it is just not called E'.
    //   HSsmall: E_0,ref = 2(1+nu_ur) G_0,ref -- literally the small-strain stiffness dynamics needs.
    //   plain HS: E_ur,ref -- the unload/reload stiffness, the closest thing it has (an under-
    //             estimate of E_0, so Vs and f_1 come out low; named in the phase message).
    // Both are then STRESS-DEPENDENT: E_ur(sigma_3) = E_ur,ref ((c cos(phi) + sigma_3 sin(phi)) /
    // (c cos(phi) + p_ref sin(phi)))^m. Freezing that at p_ref would be wrong by ~2x over a deep
    // deposit, so evaluate it at each stress point from the PARENT phase's committed stress --
    // which this phase now has (init_states). Reuses the validated hs_small_strain_params (with
    // gamma_hist = 0 = small strain, the dynamic regime) + hs_frozen_Eur, so the stiffness is the
    // model's own, not a correlation of ours. Model-family detection reads the MOTOR's models
    // table (MaterialType), never a schema enum (the B9 rule).
    std::vector<LinearElastic> gauss_elastic;
    bool has_stress_dep = false;   // HS/HSsmall/SoftSoil(Creep): stiffness from the committed stress
    for (const auto& mm : models)
        if (mm.type == MaterialType::HardeningSoil || mm.type == MaterialType::SoftSoil ||
            mm.type == MaterialType::SoftSoilCreep)
            has_stress_dep = true;
    if (has_stress_dep) {
        const int ngp = mesh.nodes_per_element == 15 ? Tri15Element::kGaussCount
                                                     : Tri6Element::kGaussCount;
        const std::vector<GaussState>* st = init_states;
        if (!st || st->size() != (size_t)mesh.element_count * ngp) {
            R.message = "Dynamic (seismic) analysis with a Hardening Soil / HSsmall / Soft Soil "
                        "material needs the stress state of a previous phase to evaluate its "
                        "stress-dependent stiffness. Put the dynamic phase after a static (K0 or "
                        "Plastic) phase.";
            return false;
        }
        gauss_elastic.resize((size_t)mesh.element_count * ngp);
        for (int e = 0; e < mesh.element_count; ++e) {
            const MaterialType mtype = models[mesh.element_material[e]].type;
            for (int g = 0; g < ngp; ++g) {
                const size_t gi = (size_t)e * ngp + g;
                LinearElastic le = mats[mesh.element_material[e]];
                if (mtype == MaterialType::HardeningSoil) {
                    // gamma_hist = 0 -> the small-strain branch (HSsmall: E_0; plain HS: E_ur).
                    const auto pe = hs_small_strain_params(
                        models[mesh.element_material[e]].hs, 0.0);
                    le.youngs_modulus =
                        hs_frozen_Eur(pe, (*st)[gi].stress, (*st)[gi].stress_zz);
                    le.poisson_ratio = pe.nu_ur;
                } else if (mtype == MaterialType::SoftSoil ||
                           mtype == MaterialType::SoftSoilCreep) {
                    // SS/SSC unload/reload stiffness frozen at the committed mean stress:
                    // K_ur = p'/kappa* -> E_ur = 3K(1 - 2 nu_ur). The linear dynamic path sees
                    // small oscillations on the u/r branch -- the same class of decision as the
                    // HS frozen-Eur above; shaking lasts seconds, so creep is negligible here.
                    const auto& mmod = models[mesh.element_material[e]];
                    const double kap = mtype == MaterialType::SoftSoil
                                           ? mmod.ssoil.kap_star
                                           : mmod.ssc.kap_star;
                    const double nu_ur = mtype == MaterialType::SoftSoil
                                             ? mmod.ssoil.nu_ur
                                             : mmod.ssc.nu_ur;
                    const auto& gsi = (*st)[gi];
                    const double p_c = std::max(
                        -(gsi.stress(0) + gsi.stress(1) + gsi.stress_zz) / 3.0,
                        softsoil::kPmin);
                    le.youngs_modulus = 3.0 * (p_c / kap) * (1.0 - 2.0 * nu_ur);
                    le.poisson_ratio = nu_ur;
                }
                gauss_elastic[gi] = le;
            }
        }
    }
    constexpr double kG = kGravity;   // gravity [m/s^2] -> mass density rho = gamma/g [Mg/m^3]
    // Mass density. Below the water table the soil is SATURATED: its total unit weight (grains +
    // pore water) is gamma_sat, and it is the TOTAL mass that carries the seismic inertia force
    // -M r a_g. Using gamma_unsat everywhere would under-predict the mass (typ. 17 vs 20 kN/m^3
    // = 15% low) -> a non-conservative inertia force. Handled per Gauss point below, exactly like
    // assemble_gravity_phreatic does for the body force. gamma_sat arrives drainage-resolved from
    // the seam, so the seismic inertia and the boundary impedance (mat_props reads rho_sat from
    // here) stay correct for a NonPorous (concrete) block too.
    std::vector<double> density(in.gamma.size(), 0.0);
    std::vector<double> density_sat(in.gamma.size(), 0.0);
    for (size_t mi = 0; mi < in.gamma.size(); ++mi) {
        density[mi] = std::fmax(0.0, in.gamma[mi]) / kG;
        const double gs = in.gamma_sat[mi] > 0.0 ? in.gamma_sat[mi] : in.gamma[mi];
        density_sat[mi] = std::fmax(0.0, gs) / kG;
    }
    const int neq = dofs.equation_count();
    if (neq == 0) {
        R.message = "Dynamic analysis: the model is fully fixed (no free DOFs). Fix only the "
                    "base and leave the sides/top free.";
        return false;
    }
    // Assemble M (consistent) and K (linear-elastic), then Rayleigh C = alpha M + beta K.
    // Soil first, then the structures on TOP of the same builders (shared translation DOFs add
    // up -> soil-structure coupling; the interface elements carry the wall<->soil joint).
    // `act` is the staged-construction mask: an excavated (or not-yet-placed) element must NOT
    // contribute stiffness OR mass. Passing it is not optional -- without it the excavation the
    // GUI draws is still full of soil in the solver, bracing the wall and carrying inertia, and
    // the reported seismic wall forces come out too LOW (unconservative) with no warning.
    katai::math::SparseMatrixBuilder bM(neq), bK(neq);
    if (in.phreatic_mass)
        assemble_mass_phreatic(mesh, dofs, density, density_sat, in.water_table_y, bM, act);
    else
        assemble_mass(mesh, dofs, density, bM, act);
    assemble_stiffness(mesh, dofs, mats, bK, act, profiles, gauss_elastic);
    assemble_structural_mass(mesh, dofs, structures, bM);
    assemble_structural_stiffness(mesh, dofs, structures, bK);
    const katai::math::CsrMatrix Mmat = bM.build(), Kmat = bK.build();
    // Phase configuration, with the same guards the driver always applied: a nonpositive or
    // degenerate value falls back exactly as before, so a null-config caller and a defaulted
    // struct produce the identical run.
    const double f1 = in.rayleigh_f1 > 0.0 ? in.rayleigh_f1 : 1.0;
    const double f2 = in.rayleigh_f2 > f1 ? in.rayleigh_f2 : (f1 * 5.0);
    const double xi = std::clamp(in.damping_ratio, 0.0, 0.5);
    const double dur = in.duration > 0.0 ? in.duration : 1.0;
    const int nst = std::clamp(in.time_steps, 1, 20000);
    const double amp = in.amplitude;
    const double freq = in.frequency > 0.0 ? in.frequency : 2.0;
    const SeismicWave wave = in.wave;
    const auto ray = rayleigh_from_modes(f1, xi, f2, xi);
    const double dt = dur / nst;
    constexpr double kPi = 3.14159265358979323846;
    const double wf = 2.0 * kPi * freq, t0 = 1.2 / freq;   // Ricker delay
    auto ag = [&](double t) {
        if (wave == SeismicWave::Record && in.accel_record.size() >= 2) {
            // Accelerogram: linear interpolation on the uniformly sampled record; ZERO after
            // the record ends (free vibration -- honest, not a wrapped/held signal).
            const auto& rec = in.accel_record;
            const double rdt = in.record_dt > 0.0 ? in.record_dt : 0.02;
            const double x = t / rdt;
            if (x <= 0.0) return amp * rec.front();
            const size_t i = (size_t)x;
            if (i + 1 >= rec.size()) return 0.0;
            const double f2 = x - (double)i;
            return amp * (rec[i] * (1.0 - f2) + rec[i + 1] * f2);
        }
        if (wave == SeismicWave::Ricker) { const double a = kPi * freq * (t - t0), a2 = a * a; return amp * (1 - 2 * a2) * std::exp(-a2); }
        return amp * std::sin(wf * t);
    };

    // --- Free-field lateral boundaries (D3b, optional) ------------------------------------
    // The sides follow the 1D free-field site response (a Lysmer dashpot C_b + a driving force
    // C_b v_ff) instead of reflecting the interior's scattered waves back in -- which otherwise
    // turns a free/roller side into a spurious box resonator. The base stays rigid (relative
    // formulation). v_ff is the 1D shear-column response of each side's own soil profile
    // (free_field.hpp), so a laterally homogeneous domain reproduces the 1D site response
    // exactly. Only FREE lateral ux DOFs are driven (a node the user fixed keeps its support).
    // Net side force is C_b (v_ff - v_2D). Verified in test_dynamics (i) + test_dynamic_gui.
    katai::math::SparseMatrixBuilder bCb(neq);
    bool has_ff = false;
    struct FFSide { FreeFieldColumnResult col; std::vector<int> ux_eq; };
    std::vector<FFSide> ff_sides;
    if (in.free_field && !mats.empty()) {
        const auto edges = katai::mesh::extract_boundary_edges(mesh, act);
        // Adjacent-material dashpot coefficients c_n = rho Vp, c_t = rho Vs (Lysmer-Kuhlemeyer).
        // Density is WATER-AWARE: below the water table the far field is saturated too -- with
        // the interior mass matrix on rho_sat, leaving the boundary impedance at gamma_unsat put
        // Vs ~8% high and the rho-Vs damping low, and the impedance mismatch with the saturated
        // interior produced spurious reflections (audit finding). (xe, ye) = edge midpoint.
        auto mat_props = [&](int mat, double xe, double ye, double& cn, double& ct,
                             double& G, double& rho) {
            mat = std::clamp(mat, 0, (int)mats.size() - 1);
            const double E = std::fmax(0.0, mats[mat].youngs_modulus);
            const double nu = std::clamp(mats[mat].poisson_ratio, 0.0, 0.49);
            const bool sat = in.has_water && ye < in.water_table_y(xe);
            rho = std::fmax(1e-12, sat ? density_sat[mat] : density[mat]);
            G = E / (2.0 * (1.0 + nu));
            const double Mmod = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));  // constrained modulus
            cn = rho * std::sqrt(Mmod / rho); ct = rho * std::sqrt(G / rho);
        };
        // Collect each side's column: unique nodes (indexed), y, and one shear segment per edge.
        struct SideBuild { std::vector<int> nodes; std::unordered_map<int, int> idx;
                           std::vector<double> y; std::vector<ShearColumnSegment> segs; };
        SideBuild left, right;
        auto add_node = [](SideBuild& s, int n2d, double yv) {
            auto it = s.idx.find(n2d); if (it != s.idx.end()) return it->second;
            const int id = (int)s.nodes.size(); s.idx[n2d] = id; s.nodes.push_back(n2d); s.y.push_back(yv);
            return id;
        };
        // A free-field boundary is only meaningful on the TRUE lateral extremes of the model --
        // the planes beyond which the soil is assumed to continue laterally uniform to infinity.
        // extract_boundary_edges returns EVERY edge owned by exactly one active element, which
        // also includes: an embedded wall's seam (split_mesh duplicates the nodes, so both faces
        // look like a boundary), an excavation face, and any slope steeper than 30 deg (|nx| =
        // sin(theta) >= 0.5). Attaching Lysmer dashpots + a 1D free-field driver to those puts an
        // absorbing "infinity" in the MIDDLE of the model, and the left/right bucketing below
        // would merge the real side chain with the interior chain into one column with
        // overlapping y and a floating second branch -> garbage v_ff, no error. So require the
        // edge to lie ON the domain's own left/right extreme.
        double xmin_d = 1e300, xmax_d = -1e300;
        for (const auto& be : edges)
            for (int i = 0; i < be.npe; ++i) {
                xmin_d = std::fmin(xmin_d, mesh.x[be.node[i]]);
                xmax_d = std::fmax(xmax_d, mesh.x[be.node[i]]);
            }
        const double xtol = 1e-6 * std::fmax(1.0, xmax_d - xmin_d);
        auto on_extreme = [&](const katai::mesh::BoundaryEdgeChain& be, double xe) {
            for (int i = 0; i < be.npe; ++i)
                if (std::fabs(mesh.x[be.node[i]] - xe) > xtol) return false;
            return true;
        };
        for (const auto& be : edges) {
            if (std::fabs(be.nx) < 0.5) continue;         // not a lateral edge (base / top / slope)
            // Every node of the edge must sit on the same extreme plane -> a genuine far-field side.
            if (!on_extreme(be, xmin_d) && !on_extreme(be, xmax_d)) continue;
            double xm = 0.0, ym = 0.0;
            for (int i = 0; i < be.npe; ++i) { xm += mesh.x[be.node[i]]; ym += mesh.y[be.node[i]]; }
            xm /= be.npe; ym /= be.npe;
            double cn, ct, G, rho; mat_props(be.material, xm, ym, cn, ct, G, rho);
            std::vector<int> chain(be.node.begin(), be.node.begin() + be.npe);
            assemble_boundary_dashpot(mesh, dofs, chain, cn, ct, bCb);   // lateral dashpot
            SideBuild& s = be.nx < 0.0 ? left : right;
            ShearColumnSegment seg; seg.nn = be.npe; seg.G = G; seg.rho = rho;
            for (int i = 0; i < be.npe; ++i)
                seg.node[i] = add_node(s, be.node[i], mesh.y[be.node[i]]);
            s.segs.push_back(seg);
        }
        for (SideBuild* s : {&left, &right}) {
            if (s->segs.empty()) continue;
            FFSide side;
            if (compliant_base) {
                // The 2D base absorbs -- the side columns must too, or they would impose the
                // rigid-base free field on an absorbing model. Halfspace impedance = the
                // side column's own deepest segment (same convention as the 2D base).
                double ymin_c = 1e300, rhovs = 0.0;
                for (const auto& sg : s->segs)
                    for (int i = 0; i < sg.nn; ++i)
                        if (s->y[sg.node[i]] < ymin_c) {
                            ymin_c = s->y[sg.node[i]];
                            rhovs = sg.rho * std::sqrt(sg.G / sg.rho);
                        }
                side.col = solve_free_field_column_compliant(
                    s->y, s->segs, ray.alpha, ray.beta, ag, dt, nst, rhovs);
            } else {
                side.col = solve_free_field_column(s->y, s->segs, ray.alpha,
                                                   ray.beta, ag, dt, nst);
            }
            if (side.col.v.empty()) continue;
            side.ux_eq.resize(s->nodes.size());
            for (size_t k = 0; k < s->nodes.size(); ++k)
                side.ux_eq[k] = dofs.equation(dofs.global_dof(s->nodes[k], 0));
            ff_sides.push_back(std::move(side));
        }
        has_ff = !ff_sides.empty();
    }
    const katai::math::CsrMatrix Cb = bCb.build();

    // --- COMPLIANT (absorbing) BASE (Joyner & Chen 1975; Sci 6.3.2 Eqn 6-12; formulation
    // locked in dynamic-seismic-formulation.md sec 11). Lysmer dashpots (rho Vs of the
    // deepest layer -- the halfspace CONTINUES it) along the bottom extreme plane, and the
    // input applied THERE through the same dashpot matrix: F(t) = C_base (2 v_up r_base),
    // with v_up = 0.5 * integral(a_g) -- the phase's a_g is the bedrock (within) motion and
    // the upward wave is HALF of it (Tut 17.8.5). Reusing C_base for the drive makes the
    // traction integration IDENTICAL to the absorption term by construction.
    katai::math::SparseMatrixBuilder bCbase(neq);
    Eigen::VectorXd rbase;
    bool has_cbase = false;
    if (compliant_base) {
        const auto edges = katai::mesh::extract_boundary_edges(mesh, act);
        double ymin_d = 1e300, ymax_d = -1e300;
        for (const auto& be : edges)
            for (int i = 0; i < be.npe; ++i) {
                ymin_d = std::fmin(ymin_d, mesh.y[be.node[i]]);
                ymax_d = std::fmax(ymax_d, mesh.y[be.node[i]]);
            }
        const double ytol2 = 1e-6 * std::fmax(1.0, ymax_d - ymin_d);
        rbase = Eigen::VectorXd::Zero(neq);
        for (const auto& be : edges) {
            if (std::fabs(be.ny) < 0.5) continue;   // not a horizontal edge (side / slope)
            bool on_base_edge = true;               // EVERY node on the bottom extreme plane
            for (int i = 0; i < be.npe; ++i)
                if (mesh.y[be.node[i]] > ymin_d + ytol2) { on_base_edge = false; break; }
            if (!on_base_edge) continue;
            const int mi2 = std::clamp(be.material, 0, (int)mats.size() - 1);
            const double E2 = std::fmax(0.0, mats[mi2].youngs_modulus);
            const double nu2 = std::clamp(mats[mi2].poisson_ratio, 0.0, 0.49);
            // A base edge below the water table uses the saturated density (same rationale as
            // the lateral dashpots: the rho-Vs impedance must match the interior mass).
            double xb = 0.0;
            for (int i = 0; i < be.npe; ++i) xb += mesh.x[be.node[i]];
            xb /= be.npe;
            const bool sat2 = in.has_water && ymin_d < in.water_table_y(xb);
            const double rho2 = std::fmax(1e-12, sat2 ? density_sat[mi2] : density[mi2]);
            const double G2 = E2 / (2.0 * (1.0 + nu2));
            const double M2 = E2 * (1.0 - nu2) / ((1.0 + nu2) * (1.0 - 2.0 * nu2));
            const double cn2 = rho2 * std::sqrt(M2 / rho2);   // rho Vp (u_y fixed -> inert in v1)
            const double ct2 = rho2 * std::sqrt(G2 / rho2);   // rho Vs (the SH absorber)
            std::vector<int> chain(be.node.begin(), be.node.begin() + be.npe);
            assemble_boundary_dashpot(mesh, dofs, chain, cn2, ct2, bCbase);
            for (int i = 0; i < be.npe; ++i) {
                const int eq2 = dofs.equation(dofs.global_dof(be.node[i], 0));
                if (eq2 >= 0) rbase(eq2) = 1.0;
            }
            has_cbase = true;
        }
        if (!has_cbase) {
            R.message = "Compliant base: no bottom boundary edges were found on the model's "
                        "lowest plane (is the base a single horizontal line?).";
            return false;
        }
    }
    const katai::math::CsrMatrix Cbase = bCbase.build();

    katai::math::SparseMatrixBuilder bC(neq);
    auto add_scaled = [&](const katai::math::CsrMatrix& A, double s) {
        for (int rr = 0; rr < A.rows; ++rr)
            for (int p = A.row_ptr[rr]; p < A.row_ptr[rr + 1]; ++p)
                bC.add_entry(rr, A.col_indices[p], s * A.values[p]);
    };
    add_scaled(Mmat, ray.alpha); add_scaled(Kmat, ray.beta);
    if (has_ff) add_scaled(Cb, 1.0);                     // lateral dashpots into the global C
    if (has_cbase) add_scaled(Cbase, 1.0);               // base absorbers into the global C
    const katai::math::CsrMatrix Cmat = bC.build();
    // Influence vector r for horizontal base shaking: 1 on every free HORIZONTAL TRANSLATION
    // DOF -- soil nodes AND a wall's independent translation DOFs (a rigid base translation
    // moves them too; skipping them would leave the wall without its inertia force).
    const Eigen::VectorXd rvec = seismic_influence_x(mesh, dofs, structures);
    const Eigen::VectorXd Mr = Mmat * rvec;
    // A model with no mass has no seismic force at all: F = -M r a_g = 0 -> the run would
    // "succeed" and report peak |u| = 0, a perfectly plausible-looking zero. That happens as
    // soon as every active material has unit weight 0 (assemble_mass skips rho == 0 elements).
    // r^T M r is the total translational mass (verified exactly in test_ssi_dynamics).
    if (!(rvec.dot(Mr) > 0.0)) {
        R.message = "Dynamic analysis: the model has no mass, so the earthquake exerts no force "
                    "(F = -M r a_g = 0). Give every active soil material a unit weight (Material "
                    "> General > gamma), and a plate its weight w.";
        return false;
    }
    // --- Model fundamental frequency f1 (small-strain elastic K, fixed base) ---------------
    // Inverse-power iteration on K x = lambda M x: y = K^{-1} M x converges to the lowest
    // mode; the influence vector r is the ideal start (the fundamental of a deposit on a
    // rigid base IS the lateral shear mode). One extra factorization + ~tens of solves --
    // negligible next to the time integration. WHY: engineers otherwise pick the Rayleigh
    // targets from the quarter-wave travel-time hand rule, which mis-places f1 by 20%+ on a
    // strong impedance-contrast profile (measured, site-response benchmark); computing the
    // model's own f1 removes the estimate entirely. Semi-definite M (massless rotations,
    // weightless plates) is harmless here: those modes have lambda -> infinity and the
    // iteration suppresses them. Failure is graceful (f1 stays 0, the note is absent).
    // With a compliant base the horizontal rigid-body mode is restrained only by dashpots:
    // K alone is SINGULAR (the SPD factorization would be garbage), so the f1 estimate is
    // skipped there -- the informational note is simply absent (honest, not wrong).
    double model_f1 = 0.0;
    if (!compliant_base) try {
        const auto ksolve = spd_factory(Kmat);
        Eigen::VectorXd x = rvec;
        x /= std::sqrt(std::fmax(x.dot(Mmat * x), 1e-300));
        double lam = 0.0;
        for (int it = 0; it < 40; ++it) {
            const Eigen::VectorXd y = ksolve(Mmat * x);
            const double ymy = y.dot(Mmat * y);
            if (!(ymy > 0.0)) { lam = 0.0; break; }
            const Eigen::VectorXd xn = y / std::sqrt(ymy);
            const double lam_new = xn.dot(Kmat * xn);   // Rayleigh quotient (M-normalized)
            const bool done = it > 2 && std::fabs(lam_new - lam) < 1e-8 * std::fabs(lam_new);
            lam = lam_new; x = xn;
            if (done) break;
        }
        if (lam > 0.0) model_f1 = std::sqrt(lam) / (2.0 * kPi);
    } catch (...) { model_f1 = 0.0; }
    R.dyn_model_f1 = model_f1;
    // Within-motion velocity v_g(t) = integral of a_g (trapezoid at the step times) -- the
    // compliant-base drive needs the VELOCITY of the upward wave, not the acceleration.
    std::vector<double> vg(nst + 1, 0.0);
    if (has_cbase)
        for (int k = 1; k <= nst; ++k)
            vg[k] = vg[k - 1] + 0.5 * dt * (ag((k - 1) * dt) + ag(k * dt));
    // RIGID base (default): F = -M r a_g, relative frame (+ free-field C_b v_ff if enabled).
    // COMPLIANT base: TOTAL-motion frame -- no -M r a_g inertia drive; the input enters as
    // the base traction 2 rho Vs v_up through the SAME dashpot matrix. Both factors are kept
    // explicit (2 x 0.5 -- Joyner-Chen doubling x within->upward halving) so neither can be
    // silently "simplified away" by a future edit.
    auto force = [&, has_ff, has_cbase](int step) {
        Eigen::VectorXd f;
        if (has_cbase) {
            const double v_up = 0.5 * vg[std::min(step, nst)];    // Tut 17.8.5: half of within
            f = Cbase * ((2.0 * v_up) * rbase);                   // Sci Eqn 6-12: factor 2
        } else
            f = -ag(step * dt) * Mr;
        if (has_ff) {
            Eigen::VectorXd vff = Eigen::VectorXd::Zero(neq);
            const int st = std::min(step, nst);
            for (const auto& side : ff_sides) {
                const Eigen::VectorXd& vc = side.col.v[st];
                for (size_t k = 0; k < side.ux_eq.size(); ++k)
                    if (side.ux_eq[k] >= 0) vff[side.ux_eq[k]] = vc[(int)k];
            }
            f += Cb * vff;
        }
        return f;
    };
    // K_eff = K + a0 M + a1 C is symmetric positive definite -> factor once, solve each step;
    // the factory is injected by the composition root (backend never named here).
    // Surface monitor node: highest (top) node closest to the horizontal centre (used by the
    // step observer below).
    // The monitor must be a node that is (a) on an ACTIVE element -- in an excavated model the
    // original ground surface may be gone -- and (b) FREE to move horizontally. If it is fixed or
    // orphaned, a_surf below would collapse to exactly a_g(t) and we would publish the INPUT
    // motion as the "surface response spectrum", overlaid on the TBDY design spectrum. That is a
    // fully plausible-looking plot of the wrong thing, so search only among admissible nodes and
    // refuse if none exists.
    const int nc = mesh.node_count;
    std::vector<char> on_active(nc, act.empty() ? 1 : 0);
    if (!act.empty())
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!act[e]) continue;
            for (int k = 0; k < mesh.nodes_per_element; ++k) on_active[mesh.node_of(e, k)] = 1;
        }
    double ytop = -1e300, xlo = 1e300, xhi = -1e300;
    bool any_surf = false;
    for (int n = 0; n < nc; ++n) {
        if (!on_active[n] || dofs.equation(dofs.global_dof(n, 0)) < 0) continue;
        ytop = std::fmax(ytop, mesh.y[n]); xlo = std::fmin(xlo, mesh.x[n]); xhi = std::fmax(xhi, mesh.x[n]);
        any_surf = true;
    }
    if (!any_surf) {
        R.message = "Dynamic analysis: no free surface node to monitor -- every node of the active "
                    "soil has its horizontal displacement fixed. Leave the ground surface free "
                    "(fix only the base) so the site response can be measured.";
        return false;
    }
    const double xmid = 0.5 * (xlo + xhi);
    int surf = -1; double bestd = 1e300;
    for (int n = 0; n < nc; ++n) {
        if (!on_active[n] || dofs.equation(dofs.global_dof(n, 0)) < 0) continue;
        if (mesh.y[n] < ytop - 1e-6 * std::fmax(1.0, std::fabs(ytop))) continue;
        const double d = std::fabs(mesh.x[n] - xmid);
        if (d < bestd) { bestd = d; surf = n; }
    }
    const int surf_eq = surf >= 0 ? dofs.equation(dofs.global_dof(surf, 0)) : -1;
    if (surf_eq < 0) {   // unreachable given the guard above; never publish a_g as the response
        R.message = "Dynamic analysis: could not place a surface monitor node.";
        return false;
    }
    // Stream every Newmark step (CONSTANT memory -- storing the full O(steps x dofs) history
    // would exhaust memory for a large mesh x many time steps and crash the run). Accumulate the
    // peak-|u| response-vs-time curve, the surface total horizontal acceleration a_x(t) (relative
    // + a_g, for the response spectrum), and keep only the deformed snapshot at peak response.
    // (Full-field time playback is deferred, P3.)
    // SEISMIC STRUCTURAL ENVELOPE: the design forces an engineer needs from a dynamic run are
    // the EXTREMES over the shaking, not the forces at any one instant -- the peak wall moment
    // does not generally occur at the peak-displacement instant, so a single snapshot would
    // UNDER-report it. Evaluate each element's diagram every step (the same force_diagram the
    // static tail uses) and keep max|N|,|Q|,|M| per station. Streaming, so still constant memory.
    // SIGNED extremes, not max|.|: the design action is static + dynamic, and the static offset
    // breaks the symmetry -- max_t |N_s + N_d(t)| = max(|N_s + max_t N_d|, |N_s + min_t N_d|),
    // which needs BOTH signed extremes of the dynamic part. (max|N_d| alone is only the answer
    // when N_s = 0.) PLAXIS tracks the same signed historical min/max pair (Ref sec 9.4.5).
    // --- Track 1a: carry the PARENT phase's structural state into the NONLINEAR increment
    // (validated + built ONCE, shared with the chained static tail -- carry_init/carry_full).
    // Without it the Coulomb cap / anchor capacity / geogrid slack act on the dynamic
    // increment about a ZERO static preload and the slip/yield ONSET is wrong wherever a
    // static preload exists (the class D7 measured at up to 8.7x demand/capacity).
    const StructuralInit& sinit = carry_init;   // empty when no carriable parent
    const bool nl_carry = dyn_nl && carry_full != nullptr;
    // Streamed committed structural state (nonlinear path): the envelope post-processors must
    // evaluate the SAME constitutive state the solver commits (the D6b rule) -- these follow
    // the solver's on_commit stream; null on the linear path (empty state passed, as before).
    const std::vector<double>* cur_anch = nullptr;
    const std::vector<double>* cur_geo = nullptr;
    const std::vector<double>* cur_if3 = nullptr;
    const std::vector<double>* cur_if5 = nullptr;
    const std::vector<double>* cur_plate = nullptr;
    const std::vector<double>* cur_plate5 = nullptr;
    // Real plastic EVENTS (nonlinear path), detected from the committed state itself after the
    // solve: the return mapping writes a new plastic value only when it actually fires (the
    // elastic branch returns the committed value bit-for-bit), so "final state != seed" is an
    // exact, rounding-proof detector. Re-evaluating the constitutive at the committed state
    // cannot do this: a station that just slipped sits exactly ON the yield boundary there
    // (|tau| == tau_max), and its `slipping` flag becomes an equality coin-flip.
    std::vector<char> if_slipped(iface_diags.size(), 0);
    std::vector<char> anch_yielded(structures.anchors.size(), 0);
    std::vector<char> plate_diag_yielded(diag_specs.size(), 0);   // per plate/wall diagram
    std::vector<StructForce> env_forces;
    struct SignedEnv { std::vector<double> nmin, nmax, qmin, qmax, mmin, mmax; };
    std::vector<SignedEnv> sgn(diag_specs.size());
    for (size_t i = 0; i < diag_specs.size(); ++i) {
        StructForce d = force_diagram(diag_specs[i], structures, mesh, dofs,
                                      Eigen::VectorXd::Zero(dofs.total_dofs()), {}, {});
        for (auto& st : d.stations) { st.N = st.Q = st.M = 0.0; }   // geometry kept, forces zeroed
        d.envelope = true;
        const size_t ns2 = d.stations.size();
        sgn[i].nmin.assign(ns2, 1e300); sgn[i].nmax.assign(ns2, -1e300);
        sgn[i].qmin.assign(ns2, 1e300); sgn[i].qmax.assign(ns2, -1e300);
        sgn[i].mmin.assign(ns2, 1e300); sgn[i].mmax.assign(ns2, -1e300);
        env_forces.push_back(std::move(d));
    }
    // Same envelope for the wall<->soil Coulomb joints: the seismic shear / normal stress the
    // interface actually carries. ELASTIC branch (see force_diagram) -- the dynamic system
    // solves the joint with k_n,k_s, so reporting a Coulomb-capped tau would contradict the
    // displacements it produced.
    std::vector<InterfaceResult> env_ifaces;
    // util: per-instant demand/capacity running max (nonlinear-carry path only -- both tau and
    // sigma_n from the SAME instant, so the ratio is the one the solver actually enforced;
    // pairing signed envelope extremes from different instants could exceed 1 spuriously).
    struct SignedIfEnv { std::vector<double> tmin, tmax, smin, smax, dmax, gmax, util; };
    std::vector<SignedIfEnv> sgn_if(iface_diags.size());
    for (size_t i = 0; i < iface_diags.size(); ++i) {
        InterfaceResult ir = force_diagram(iface_diags[i], structures, mesh, dofs,
                                           Eigen::VectorXd::Zero(dofs.total_dofs()),
                                           {}, {}, iface_elastic);
        for (auto& st : ir.stations) { st.tau = st.sigma_n = st.slip = st.gap = 0.0; }
        ir.envelope = true;
        const size_t ns2 = ir.stations.size();
        sgn_if[i].tmin.assign(ns2, 1e300); sgn_if[i].tmax.assign(ns2, -1e300);
        sgn_if[i].smin.assign(ns2, 1e300); sgn_if[i].smax.assign(ns2, -1e300);
        sgn_if[i].dmax.assign(ns2, 0.0);   sgn_if[i].gmax.assign(ns2, 0.0);
        sgn_if[i].util.assign(ns2, 0.0);
        env_ifaces.push_back(std::move(ir));
    }
    double peak = 0.0; Eigen::VectorXd peak_uf;
    static const std::vector<double> kNoPlastic;   // linear path: empty committed state, as before
    auto observer = [&](int, double t, const Eigen::VectorXd& u,
                        const Eigen::VectorXd&, const Eigen::VectorXd& a) {
        const Eigen::VectorXd uf = expand_to_full(dofs, u);
        double mx = 0.0;
        for (int n = 0; n < nc; ++n) mx = std::fmax(mx, std::hypot(uf[2 * n], uf[2 * n + 1]));
        R.consol_time.push_back(t);
        R.consol_settlement.push_back(mx);      // reused: peak |u| over the domain at time t
        R.consol_excess_pore.push_back(0.0);
        // Track 1a: with the parent state carried, the solver's structural elements live at
        // TOTAL displacement (parent datum + dynamic) with the streamed committed plastic
        // state -- so evaluate the diagrams there too (same constitutive as the solve, D6b).
        // The stations are then the TOTAL action; no superposition happens afterwards.
        Eigen::VectorXd uf_tot;
        if (nl_carry) uf_tot = *carry_full + uf;
        const Eigen::VectorXd& uf_s = nl_carry ? uf_tot : uf;
        const std::vector<double>& anch_c = cur_anch ? *cur_anch : kNoPlastic;
        const std::vector<double>& geo_c = cur_geo ? *cur_geo : kNoPlastic;
        const std::vector<double>& if3_c = cur_if3 ? *cur_if3 : kNoPlastic;
        const std::vector<double>& if5_c = cur_if5 ? *cur_if5 : kNoPlastic;
        const std::vector<double>& pl_c = cur_plate ? *cur_plate : kNoPlastic;
        const std::vector<double>& pl5_c = cur_plate5 ? *cur_plate5 : kNoPlastic;
        for (size_t i = 0; i < diag_specs.size(); ++i) {
            // Linear path: uncapped anchors + full-EA geogrids + elastic plates, matching the
            // linear solve (iface_elastic == !dyn_nl is exactly that flag). Nonlinear path: the
            // streamed committed state + the real caps, matching the per-step Newton solve.
            const StructForce cur = force_diagram(diag_specs[i], structures, mesh, dofs, uf_s,
                                                  anch_c, geo_c, iface_elastic, pl_c, pl5_c);
            SignedEnv& s = sgn[i];
            const size_t n2 = std::min(cur.stations.size(), s.nmin.size());
            for (size_t k = 0; k < n2; ++k) {
                s.nmin[k] = std::fmin(s.nmin[k], cur.stations[k].N);
                s.nmax[k] = std::fmax(s.nmax[k], cur.stations[k].N);
                s.qmin[k] = std::fmin(s.qmin[k], cur.stations[k].Q);
                s.qmax[k] = std::fmax(s.qmax[k], cur.stations[k].Q);
                s.mmin[k] = std::fmin(s.mmin[k], cur.stations[k].M);
                s.mmax[k] = std::fmax(s.mmax[k], cur.stations[k].M);
            }
        }
        for (size_t i = 0; i < iface_diags.size(); ++i) {
            const InterfaceResult cur = force_diagram(iface_diags[i], structures, mesh, dofs,
                                                      uf_s, if3_c, if5_c, iface_elastic);
            SignedIfEnv& s = sgn_if[i];
            const auto& props = (iface_diags[i].order == 15)
                                    ? structures.interfaces5[iface_diags[i].begin].props
                                    : structures.interfaces[iface_diags[i].begin].props;
            const size_t n2 = std::min(cur.stations.size(), s.tmin.size());
            for (size_t k = 0; k < n2; ++k) {
                s.tmin[k] = std::fmin(s.tmin[k], cur.stations[k].tau);
                s.tmax[k] = std::fmax(s.tmax[k], cur.stations[k].tau);
                s.smin[k] = std::fmin(s.smin[k], cur.stations[k].sigma_n);
                s.smax[k] = std::fmax(s.smax[k], cur.stations[k].sigma_n);
                s.dmax[k] = std::fmax(s.dmax[k], std::fabs(cur.stations[k].slip));
                s.gmax[k] = std::fmax(s.gmax[k], std::fabs(cur.stations[k].gap));
                if (nl_carry) {
                    // Same-instant demand/capacity: the solver capped |tau| at tau_max(sigma_n)
                    // THIS instant, so the ratio is <= ~1 by construction wherever Coulomb
                    // governed -- it shows where the joint rode its capacity, and the margin
                    // elsewhere. (tau_max = 0 = separated: tau is 0 too, so report 0.)
                    const double tmax2 =
                        std::fmax(0.0, props.c_i - cur.stations[k].sigma_n * std::tan(props.phi_i));
                    const double ut = tmax2 > 1e-12
                                          ? std::fabs(cur.stations[k].tau) / tmax2
                                          : (std::fabs(cur.stations[k].tau) > 1e-12 ? kUtilCap : 0.0);
                    s.util[k] = std::fmax(s.util[k], std::fmin(ut, kUtilCap));
                }
            }
        }
        if (mx > peak) { peak = mx; peak_uf = uf; }
        // Total surface acceleration: the rigid-base path solves RELATIVE motion (add a_g);
        // the compliant-base path solves TOTAL motion already (adding a_g would silently
        // double-count the ground motion -- the formulation-lock's trap #4).
        const double a_surf = (surf_eq >= 0 ? a[surf_eq] : 0.0) + (has_cbase ? 0.0 : ag(t));
        R.dyn_time.push_back(t);
        R.dyn_surface_ax.push_back(a_surf);
        R.dyn_peak_surface_a = std::fmax(R.dyn_peak_surface_a, std::fabs(a_surf));
    };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
    // Initial acceleration, in closed form: starting from rest the whole system rides the base
    // rigidly, a(0) = -r a_g(0)  ->  M a(0) = -M r a_g(0) = F(0) EXACTLY, for any M. Passing it
    // means M is never factorized -- and with structures M IS generally singular (a weightless
    // plate, and every plate rotation DOF, leave zero rows). Factorizing that would corrupt the
    // solver. (Free-field: v_ff(0) = 0 from rest, so F(0) = -M r a_g(0) still holds exactly.)
    // Compliant base: TOTAL-motion frame from rest -- v_g(0) = 0 so F(0) = 0 and a(0) = 0
    // exactly (no rigid ride at t = 0). Rigid base: the closed-form relative-frame start.
    const Eigen::VectorXd a_init =
        has_cbase ? Eigen::VectorXd(Eigen::VectorXd::Zero(neq))
                  : Eigen::VectorXd(-ag(0.0) * rvec);
    // LINEAR (default) vs NONLINEAR (opt-in) time integration. Both STREAM to the same observer,
    // so every envelope / response-spectrum post-processing below is shared -- only the internal-
    // force model differs. The nonlinear path assembles K_T(u) from the constitutive each Newton
    // iteration (generally NONsymmetric -> a nonsymmetric solve), continues from the parent
    // stress (init_states) with f_int0 subtracted, and can plastify soil / slip interfaces /
    // yield structures during shaking. It is much slower (the tangent refactors every iteration).
    bool nl_ok = true; int nl_steps = nst;
    std::vector<GaussState> nl_gauss;   // nonlinear path: final committed effective stress
    if (dyn_nl) {
        const LinearSolve nl_solver = make_nonsym_solver();
        const std::vector<GaussState> nl_init =
            init_states ? *init_states : std::vector<GaussState>{};
        NewmarkNonlinearOptions nlopt;   // gamma=1/2, beta=1/4, 30 Newton, tol 1e-6
        // Stream the committed structural state to the envelope observer (fires right before
        // each on_step, incl. step 0 = the seeded parent state) -- the diagrams then use the
        // solver's REAL constitutive state, never a stale or zero one (D6b rule).
        nlopt.on_commit = [&](int, const NewmarkNonlinearResult& s) {
            cur_anch = &s.anchor_plastic; cur_geo = &s.geogrid_plastic;
            cur_if3 = &s.interface_slip; cur_if5 = &s.interface5_slip;
            cur_plate = &s.plate_plastic; cur_plate5 = &s.plate5_plastic;
        };
        auto nlres = solve_newmark_nonlinear(
            mesh, dofs, models, Mmat, Cmat, force, dt, nst, z, z, a_init,
            nl_solver, observer, nl_init, act, structures, profiles, sinit, nlopt);
        nl_ok = nlres.converged; nl_steps = nlres.steps_completed;
        // Real plastic events = the committed state moved off its seed. NOT bit-for-bit: a
        // station the PARENT left exactly ON the yield surface re-evaluates at |tau| ==
        // tau_max, and last-bit rounding can push it to the plastic branch, which rewrites
        // slip_p with ~ulp noise while tau stays at the cap -- no physical slip, but a
        // bitwise detector would brand the undriven joint [SLIPPING]. A picometer threshold
        // (1e-12 m, far below engineering meaning, far above ulp of any real displacement)
        // kills that noise and keeps every physical event.
        auto state_changed = [](const std::vector<double>& fin, const std::vector<double>& seed,
                                size_t b, size_t e2) {
            for (size_t k = b; k < e2 && k < fin.size(); ++k) {
                const double s0 = k < seed.size() ? seed[k] : 0.0;
                if (std::fabs(fin[k] - s0) > 1e-12 + 1e-9 * std::fabs(s0)) return true;
            }
            return false;
        };
        for (size_t i = 0; i < iface_diags.size(); ++i) {
            const auto& is = iface_diags[i];
            if_slipped[i] = is.order == 15
                ? state_changed(nlres.interface5_slip, sinit.interface5_slip,
                                is.begin * (size_t)iface::kPointCount5,
                                is.end * (size_t)iface::kPointCount5)
                : state_changed(nlres.interface_slip, sinit.interface_slip,
                                is.begin * (size_t)iface::kPointCount,
                                is.end * (size_t)iface::kPointCount);
        }
        for (size_t ai = 0; ai < structures.anchors.size(); ++ai)
            anch_yielded[ai] = state_changed(nlres.anchor_plastic, sinit.anchor_plastic, ai, ai + 1);
        // Plate M-N hinge really formed during shaking (same committed-state event detector;
        // the eps_p/kap_p values are strains ~1e-3, so the 1e-12 threshold is equally safe).
        for (size_t i = 0; i < diag_specs.size(); ++i) {
            const auto& sp2 = diag_specs[i];
            if (sp2.kind == 0)
                plate_diag_yielded[i] = state_changed(
                    nlres.plate_plastic, sinit.plate_plastic,
                    sp2.begin * (size_t)plate::kPlasticStateSize,
                    sp2.end * (size_t)plate::kPlasticStateSize);
            else if (sp2.kind == 5)
                plate_diag_yielded[i] = state_changed(
                    nlres.plate5_plastic, sinit.plate5_plastic,
                    sp2.begin * (size_t)plate::kPlasticStateSize5,
                    sp2.end * (size_t)plate::kPlasticStateSize5);
        }
        nl_gauss = std::move(nlres.gauss_states);
    } else {
        solve_newmark(Mmat, Cmat, Kmat, force, dt, nst, z, z, 0.5, 0.25, spd_factory,
                      observer, a_init);
    }
    R.disp = peak_uf.size() ? peak_uf : expand_to_full(dofs, z);
    R.max_disp = peak;
    // --- Superpose the parent phase's STATIC action -------------------------------------------
    // The dynamic solve is the LINEAR increment about the static state (small-strain unload /
    // reload), so the action a section is designed for is static + dynamic. Superposition is
    // exact for the increment: K u_s = f_s and M u_d'' + C u_d' + K u_d = -M r a_g sum to the
    // total equation, because u_s is time-invariant. The parent must be a real static state --
    // an envelope (a Dynamic parent) is not one, and the diagram must line up station for
    // station (same element, same count); otherwise report the dynamic action alone and say so.
    const SolveResult* base = parent;
    const bool base_ok = base && base->ok && !base->struct_forces.empty() &&
                         !base->struct_forces.front().envelope;
    auto static_of = [&](const std::string& name, size_t nst2) -> const StructForce* {
        if (!base_ok) return nullptr;
        for (const auto& sf : base->struct_forces)
            if (sf.name == name && sf.stations.size() == nst2) return &sf;
        return nullptr;
    };
    auto static_if_of = [&](const std::string& name, size_t nst2) -> const InterfaceResult* {
        if (!base || !base->ok) return nullptr;
        for (const auto& ir : base->interface_forces)
            if (ir.name == name && ir.stations.size() == nst2 && !ir.envelope) return &ir;
        return nullptr;
    };
    bool superposed_any = false, superposed_all = true;
    // Design extreme of a static offset plus a signed dynamic range.
    auto design = [](double s, double lo, double hi) {
        return std::fmax(std::fabs(s + lo), std::fabs(s + hi));
    };
    for (size_t i = 0; i < env_forces.size(); ++i) {
        StructForce& e = env_forces[i];
        const StructForce at_peak = force_diagram(diag_specs[i], structures, mesh, dofs,
                                                  R.disp, {}, {});
        // Track 1a: with the parent state CARRIED the stations already hold the total action
        // (the solver evaluated the structures at parent datum + dynamic increment) --
        // superposing the parent's statics again would double-count them. `superposed` keeps
        // its meaning ("the stations are the TOTAL design action"), just via the solver now.
        const StructForce* st0 = nl_carry ? nullptr : static_of(e.name, e.stations.size());
        if (st0 || nl_carry) superposed_any = true; else superposed_all = false;
        e.superposed = st0 != nullptr || nl_carry;
        // Anchor really yielded during shaking (nonlinear path; exact committed-state event).
        if (dyn_nl && diag_specs[i].kind == 1 && diag_specs[i].begin < anch_yielded.size() &&
            anch_yielded[diag_specs[i].begin])
            e.yielded = true;
        // Plate M-N hinge formed during shaking (nonlinear path; committed-state event).
        if (dyn_nl && (diag_specs[i].kind == 0 || diag_specs[i].kind == 5) &&
            plate_diag_yielded[i])
            e.yielded = true;
        for (size_t k = 0; k < e.stations.size(); ++k) {
            if (k < at_peak.stations.size()) {   // deformed-overlay placement only
                e.stations[k].ux = at_peak.stations[k].ux;
                e.stations[k].uy = at_peak.stations[k].uy;
            }
            const double sN = st0 ? st0->stations[k].N : 0.0;
            const double sQ = st0 ? st0->stations[k].Q : 0.0;
            const double sM = st0 ? st0->stations[k].M : 0.0;
            e.stations[k].N = design(sN, sgn[i].nmin[k], sgn[i].nmax[k]);
            e.stations[k].Q = design(sQ, sgn[i].qmin[k], sgn[i].qmax[k]);
            e.stations[k].M = design(sM, sgn[i].mmin[k], sgn[i].mmax[k]);
        }
        const auto env = force_envelope(e.stations);
        e.max_N = env.max_abs_N; e.max_Q = env.max_abs_Q; e.max_M = env.max_abs_M;
        R.struct_forces.push_back(std::move(e));
    }
    // --- Joints: total tau / sigma_n + the Coulomb DEMAND/CAPACITY the linear solve never checks.
    // The dynamic joint is elastic, so it cannot slip and never reports slip. That is exactly why
    // the utilisation must be reported: |tau_total| / tau_max > 1 means the real joint WOULD have
    // slipped and the linear result is unconservative there. With the static state superposed the
    // check is rigorous -- tau_max = max(0, c_i - sigma_n_total tan(phi_i)) needs the TOTAL normal
    // stress (static incl. sigma_n0 + dynamic), which is precisely what the parent phase supplies.
    for (size_t i = 0; i < env_ifaces.size(); ++i) {
        InterfaceResult& e = env_ifaces[i];
        const InterfaceResult* st0 = nl_carry ? nullptr : static_if_of(e.name, e.stations.size());
        if (st0 || nl_carry) superposed_any = true; else superposed_all = false;
        e.superposed = st0 != nullptr || nl_carry;
        const auto& is = iface_diags[i];
        double over = 0.0;
        for (size_t k = 0; k < e.stations.size(); ++k) {
            const double sT = st0 ? st0->stations[k].tau : 0.0;
            const double sS = st0 ? st0->stations[k].sigma_n : 0.0;
            e.stations[k].tau = design(sT, sgn_if[i].tmin[k], sgn_if[i].tmax[k]);
            e.stations[k].sigma_n = design(sS, sgn_if[i].smin[k], sgn_if[i].smax[k]);
            e.stations[k].slip = sgn_if[i].dmax[k];
            e.stations[k].gap = sgn_if[i].gmax[k];
            // This superposed demand/capacity ratio is a LINEAR-path concept: it measures WHERE
            // an elastic joint that cannot slip would be unconservative. The NONLINEAR path
            // never uses it: WITH the parent state carried (Track 1a) the per-instant ratio is
            // accumulated in the observer instead (the cap acted on the total, so it is exact
            // and <= ~1 by construction); WITHOUT a carriable parent the increment is capped
            // about a zero static-shear start, so superposing the parent's static tau onto the
            // capped increment would double-count -- no ratio there, the slip is the measure.
            if (!st0 || dyn_nl) continue;
            // Worst utilisation over the two signed extremes, each with its own normal stress
            // (tension-positive: compression raises tau_max).
            const auto& props = (is.order == 15) ? structures.interfaces5[is.begin].props
                                                 : structures.interfaces[is.begin].props;
            double util = 0.0;
            for (int c = 0; c < 2; ++c) {
                const double tt = sT + (c == 0 ? sgn_if[i].tmin[k] : sgn_if[i].tmax[k]);
                const double sn = sS + (c == 0 ? sgn_if[i].smin[k] : sgn_if[i].smax[k]);
                const double tmax = std::fmax(0.0, props.c_i - sn * std::tan(props.phi_i));
                const double u = tmax > 1e-12 ? std::fabs(tt) / tmax
                                              : (std::fabs(tt) > 1e-12 ? kUtilCap : 0.0);
                util = std::fmax(util, std::fmin(u, kUtilCap));
            }
            e.stations[k].utilisation = util;
            e.max_utilisation = std::max(e.max_utilisation, util);
            if (util > 1.0) over += 1.0;
        }
        if (nl_carry) {
            // Track 1a: per-instant demand/capacity from the observer (tau and sigma_n of the
            // SAME instant, the ratio the solver actually enforced on the TOTAL action). It is
            // <= ~1 wherever Coulomb governed -- it shows where the joint rode its capacity
            // and the margin everywhere else; over_fraction stays ~0 because the solver
            // RESOLVES the exceedance the linear path could only measure.
            double over_nl = 0.0;
            for (size_t k = 0; k < e.stations.size() && k < sgn_if[i].util.size(); ++k) {
                e.stations[k].utilisation = sgn_if[i].util[k];
                e.max_utilisation = std::max(e.max_utilisation, sgn_if[i].util[k]);
                if (sgn_if[i].util[k] > 1.0) over_nl += 1.0;
            }
            if (!e.stations.empty()) e.over_fraction = over_nl / (double)e.stations.size();
        } else if (st0 && !e.stations.empty()) {
            e.over_fraction = over / (double)e.stations.size();
        }
        for (const auto& st : e.stations) {   // stations now hold the TOTAL design values
            e.max_abs_tau = std::max(e.max_abs_tau, st.tau);
            e.max_abs_sigma_n = std::max(e.max_abs_sigma_n, st.sigma_n);
            e.max_abs_slip = std::max(e.max_abs_slip, st.slip);
        }
        // LINEAR path: the dynamic joint is elastic -> it cannot slip by construction, so
        // any_slip stays false and the UI shows "[elastic, no slip check]". NONLINEAR path:
        // any_slip = the joint's committed plastic slip really moved during the shaking (an
        // exact state event). It is NOT derived from the relative shear displacement (`slip`
        // stations): an elastic, fully bonded joint always has a nonzero du_s (= tau / k_s),
        // and calling that "slipping" would brand every bonded joint with a false [SLIPPING].
        e.any_slip = dyn_nl && i < if_slipped.size() && if_slipped[i];
        e.slip_checked = dyn_nl;   // Coulomb branch ran -> [SLIPPING]/[bonded] are real findings
        R.interface_forces.push_back(std::move(e));
    }
    // 5%-damped elastic response spectrum of the surface motion + TBDY 2018 design spectrum
    // [both m/s^2 over T = 0.05..4 s], for the code-comparison overlay (docs/references/
    // tbdy-2018-seismic.md). Design spectrum from the phase's map coefficients + site class.
    for (int i = 0; i <= 40; ++i) R.dyn_period.push_back(0.05 * std::pow(4.0 / 0.05, i / 40.0));
    R.dyn_response_sa = response_spectrum(R.dyn_surface_ax, dt, R.dyn_period, 0.05);
    {
        const double ss = in.tbdy_ss, s1 = in.tbdy_s1;
        const int scls = std::clamp(in.site_class, 0, 4);
        const auto dcoef = tbdy_design_coefficients((SiteClass)scls, ss, s1);
        // EC8 overlay (optional): Se(T) in m/s^2, a_g = gamma_I * a_gR [g] * 9.81.
        if (in.ec8_enabled) {
            const double ag = std::fmax(0.0, in.ec8_gamma) *
                              std::fmax(0.0, in.ec8_agr) * kG;
            const auto gnd = (Ec8GroundType)std::clamp(in.ec8_ground, 0, 4);
            const auto typ = in.ec8_type == 1 ? Ec8SpectrumType::Type2
                                              : Ec8SpectrumType::Type1;
            for (double T : R.dyn_period)
                R.dyn_design_sa_ec8.push_back(
                    ec8_elastic_spectrum(ag, gnd, typ, T));
        }
        for (double T : R.dyn_period)
            R.dyn_design_sa.push_back(tbdy_elastic_spectrum(dcoef.F_S, dcoef.F_1, T) * kG);
    }
    // Seismic effective-stress recovery. The LINEAR path solves K u about a zero-stress datum, so
    // it does not evolve the stress field -> deferred (zero). The NONLINEAR path carries the
    // committed effective stress through the shaking (parent stress + plastic evolution), so
    // recover the nodal field from its final committed Gauss state: the post-earthquake stress
    // state, showing where the soil plastified. (Uses the same recovery as the static branches.)
    if (dyn_nl && !nl_gauss.empty()) {
        R.stress = recover_nodal_stresses_from_gauss(mesh, nl_gauss, act);
    } else {
        R.stress.stress.assign(nc, Eigen::Vector3d::Zero());   // linear path: stress recovery deferred
        // Zeros are a legitimate datum here and an illegitimate ANSWER: a reader who plots the
        // stress field of a linear Dynamic phase sees an unstressed soil, which is not what the
        // analysis found -- it is what the analysis never computed.
        add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A001", "",
                       "Linear dynamic phase: the stress field is NOT recovered and is reported "
                       "as zero everywhere. Displacements, accelerations and structural forces "
                       "are the results of this phase; for stresses run it nonlinear.");
    }
    R.pore.assign(nc, 0.0);
    R.load_factor = 1.0; R.iterations = nst;
    R.active = act;
    R.ok = true;
    char buf[220];
    if (dyn_nl)
        std::snprintf(buf, sizeof(buf),
                      "Nonlinear dynamic analysis: %d/%d steps over %.3g s, peak |u| = %.4g m.%s",
                      nl_steps, nst, dur, peak,
                      nl_ok ? "" : " STOPPED early -- a step did not converge; the result is partial.");
    else
        std::snprintf(buf, sizeof(buf), "Dynamic analysis: %d steps over %.3g s, peak |u| = %.4g m.",
                      nst, dur, peak);
    R.message = buf;
    if (has_cbase)
        R.message += " COMPLIANT BASE: the bottom boundary absorbs outgoing waves (Lysmer "
                     "rho-Vs dashpots of the deepest layer -- the halfspace continues it) and "
                     "the motion is TOTAL (the base moves). Input convention: a_g(t) is the "
                     "bedrock (within) motion; the upward wave = HALF of it is applied "
                     "(Joyner-Chen, factor 2 on the dashpot reaction). v1 scope: horizontal "
                     "SH only (base u_y stays fixed).";
    // f1 + Rayleigh-band guardrail. The informational sentence always appears (the engineer
    // gets the model's own f1 instead of estimating it by hand); the NOTE appears only when
    // the band genuinely mis-damps the fundamental -- a properly bracketing band (f_R1 ~ f1,
    // f_R2 ~ 3 f1) keeps the in-band deviation <= ~13% and never trips the 20% threshold,
    // so this is a measure, not paranoia.
    if (model_f1 > 0.0) {
        const double xi_eff = ray.alpha / (4.0 * kPi * model_f1) + ray.beta * kPi * model_f1;
        char fb[360];
        std::snprintf(fb, sizeof(fb),
                      " Model fundamental frequency (small-strain, fixed base): f1 = %.3g Hz; "
                      "effective Rayleigh damping there %.2f%% (target %.2f%% at %.3g-%.3g Hz).",
                      model_f1, 100.0 * xi_eff, 100.0 * xi, f1, f2);
        R.message += fb;
        const bool outside = model_f1 < f1 || model_f1 > f2;
        if (xi > 0.0 && (outside || std::fabs(xi_eff - xi) > 0.20 * xi)) {
            std::snprintf(fb, sizeof(fb),
                          " NOTE: the Rayleigh band does not properly cover f1 -- the model's own "
                          "resonance is damped at %.2f%% instead of the intended %.2f%%. Set "
                          "rayleigh_f1 near %.3g Hz (and rayleigh_f2 ~ 3x that); do not rely on the "
                          "quarter-wave travel-time estimate, it can be 20%%+ off on layered profiles.",
                          100.0 * xi_eff, 100.0 * xi, model_f1);
            R.message += fb;
        }
    }
    // The dynamic system is LINEAR, so the structural elements are linearised about the taut /
    // non-slipping state (their incremental stiffness under a static preload). Name the
    // simplifications that actually apply to THIS model instead of leaving them in the manual:
    // a SLACK geogrid is the sharp case -- the static path drops it entirely while the linear
    // branch keeps EA (test_ssi_dynamics (e) measures that gap).
    if (!R.struct_forces.empty() || !R.interface_forces.empty())
        R.message += superposed_all ? " Structural forces are the TOTAL design action (static + dynamic)."
                   : superposed_any ? " Some structural forces are TOTAL (static + dynamic); the rest are"
                                      " the dynamic action alone -- superpose them yourself."
                                    : " Structural forces are the DYNAMIC action alone (no static parent"
                                      " phase to add) -- superpose your static phase yourself.";
    if (!dyn_nl) {
        std::string lin;
        auto add = [&](const char* s) { lin += lin.empty() ? " Linearised: " : ", "; lin += s; };
        if (!structures.geogrids.empty()) add("geogrids stay elastic in compression (no tension cut-off)");
        if (!structures.anchors.empty()) add("anchors do not yield");
        if (!structures.interfaces.empty() || !structures.interfaces5.empty())
            add("interfaces do not slip");
        if (!lin.empty()) {
            R.message += lin + ".";
            // The same sentence, machine-tagged: a front end or a script must be able to act on
            // "this run linearised structural behaviour" without parsing the message.
            add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A002", "",
                           "Linear dynamic phase, with structural elements present." + lin +
                               ". Their capacity is not checked during the shaking; run the phase "
                               "nonlinear where yielding or slip governs.");
        }
    } else if (!structures.geogrids.empty() || !structures.anchors.empty() ||
               !structures.interfaces.empty() || !structures.interfaces5.empty()) {
        // NONLINEAR path: the very simplifications the linear branch had to disclose are now solved
        // for -- interfaces slip, geogrids go slack, anchors yield, soil plastifies during shaking.
        // Whether the STATIC preload took part in the caps is a safety-relevant fact: say it.
        R.message += nl_carry
            ? " Fully nonlinear, continuing from the parent phase's structural state: the"
              " Coulomb / yield / slack checks act on the TOTAL (static + dynamic) action --"
              " interfaces may slip, geogrids may go slack, anchors may yield and the soil may"
              " plastify during shaking (per-step Newton)."
            : " Fully nonlinear: interfaces may slip, geogrids may go slack, anchors may"
              " yield and the soil may plastify during shaking (per-step Newton). NOTE: no"
              " parent structural state was available to this phase, so the structural checks"
              " act on the dynamic increment about a ZERO static preload -- slip/yield onset"
              " ignores any static preload. Put the dynamic phase directly after the static"
              " phase it continues from.";
    }
    // The one linearisation we can now MEASURE instead of merely disclosing: where the total
    // Coulomb demand exceeds the capacity, the real joint slips and this elastic result is
    // unconservative. Report it as a measure, not a verdict -- sigma_n -> 0 at the ground
    // surface, so a short zone at the top of any wall goes over under real shaking. Give the
    // worst ratio AND how much of the joint is over, and let the engineer judge. LINEAR path only:
    // the nonlinear solver actually slips the joint (tau is Coulomb-capped), so it is not
    // unconservative there -- the exceedance is resolved, not merely disclosed.
    for (const auto& ir : R.interface_forces)
        if (!dyn_nl && ir.superposed && ir.max_utilisation > 1.0) {
            char ub[256];
            if (ir.max_utilisation >= kUtilCap)
                std::snprintf(ub, sizeof(ub),
                              " Interface \"%s\": %.0f%% of the joint exceeds its Coulomb capacity "
                              "(it separates / has no shear capacity somewhere). The real joint would "
                              "slip there; this linear analysis cannot, so it is unconservative there.",
                              ir.name.c_str(), 100.0 * ir.over_fraction);
            else
                std::snprintf(ub, sizeof(ub),
                              " Interface \"%s\": peak Coulomb demand/capacity = %.2fx over %.0f%% of "
                              "the joint. The real joint would slip there; this linear analysis cannot, "
                              "so it is unconservative there.",
                              ir.name.c_str(), ir.max_utilisation, 100.0 * ir.over_fraction);
            R.message += ub;
        }
    // Same honesty for the plate M-N hinge: the LINEAR path solves plates elastic, so where
    // the envelope's |N|/Np + |M|/Mp exceeds 1 the real plate would form a plastic hinge and
    // redistribute -- measured and disclosed (envelope extremes pair different instants, so
    // this is an upper-bound measure, not a verdict). The nonlinear path RESOLVES it (the
    // return mapping caps M/N step by step), so no warning there.
    if (!dyn_nl)
        for (size_t i = 0; i < diag_specs.size() && i < R.struct_forces.size(); ++i) {
            const auto& sp2 = diag_specs[i];
            if (sp2.kind != 0 && sp2.kind != 5) continue;
            const auto& pp2 = sp2.kind == 0 ? structures.plates[sp2.begin].props
                                            : structures.plates5[sp2.begin].props;
            if (!pp2.plastic()) continue;
            const double iNp2 = pp2.Np > 0.0 ? 1.0 / pp2.Np : 0.0;
            const double iMp2 = pp2.Mp > 0.0 ? 1.0 / pp2.Mp : 0.0;
            double umax = 0.0;
            for (const auto& stn : R.struct_forces[i].stations)
                umax = std::max(umax, std::fabs(stn.N) * iNp2 + std::fabs(stn.M) * iMp2);
            if (umax > 1.0) {
                char ub[256];
                std::snprintf(ub, sizeof(ub),
                              " Plate \"%s\": envelope M-N demand reaches %.2fx the Mp/Np yield "
                              "diamond. The real plate would form a plastic hinge there; this "
                              "linear dynamic run keeps it elastic (enable nonlinear dynamics to "
                              "resolve the hinge).",
                              sp2.name.c_str(), umax);
                R.message += ub;
            }
        }
    return true;
}

} // namespace katai::core
