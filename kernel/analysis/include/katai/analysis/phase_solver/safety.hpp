#pragma once
// Safety phase strategy (Stage B9). Phi-c reduction / strength reduction
// method, the "Safety" phase: bisect the strength reduction factor to the slope's
// factor of safety, and show the FAILURE MECHANISM (displacement localized
// along the slip surface). Needs shear strength (Mohr-Coulomb); linear-elastic
// soil never fails.
//
// Same contract as the other phase strategies: neutral inputs resolved at the
// caller's seam (the model-family flags are registry-derived products of the
// common setup), refusal messages engine-owned and byte-identical, and the
// linear solver enters as a callback built by the composition root, never
// named here. On success this strategy does NOT set R.ok or R.message -- the
// phase falls through to the driver's common result tail, exactly as it
// always has.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/analysis/results.hpp>
#include <katai/analysis/strength_reduction.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// The phase's neutral configuration.
struct SafetyPhase {
    bool nonlinear_soil = false;    // any nonlinear constitutive model present
    bool has_hardening = false;     // hardening family present
    bool has_softsoil = false;      // soft-soil family present
    bool axisymmetric = false;      // r-z mode; refused (not available yet)
    std::vector<char> active;       // element activity; empty = everything active
    // Numerical controls for the TRIAL solves inside the strength-reduction search; 0 = this
    // strategy's own defaults below. Each trial is a means to an end (does this reduced-strength
    // slope still stand?), so its stopping rule is deliberately looser than a static phase's --
    // but it is still a stopping rule, and a study asking whether the factor of safety depends
    // on it must be able to CHANGE it. Until 2026-08-09 these were hard-coded here, and the
    // caller's override was silently ignored: a sweep over the tolerance then returned three
    // identical numbers because nothing had changed, which reads exactly like independence.
    double tolerance = 0.0;
    int load_steps = 0;
    int max_iterations = 0;
    // THE PARENT PHASE'S STATE, when the search starts from it (see safety_analysis): its
    // committed stresses, its internal force (soil + structures, equation space) and the
    // structures' carried state. All empty = the search re-solves the ground from the unstressed
    // state, as a Safety run with no phase before it must.
    std::vector<GaussState> parent_state;
    Eigen::VectorXd parent_force;
    StructuralInit parent_structures;
    bool from_parent() const { return !parent_state.empty() && parent_force.size() > 0; }
};

// THE ADVANCED MODELS ENTER A SAFETY SEARCH AS MOHR-COULOMB. Reducing the strength of a
// hardening model moves its yield surfaces and its hardening together, and the factor the search
// ends on then depends on the path it happened to take: measured on an embankment on Hardening
// Soil rebuilt from a published tutorial, the same model reported 1.864 along one sequence of
// reduction steps and 1.676 along another, depending only on whether one step converged. The
// strength reduction is defined on the Mohr-Coulomb strength these models share (c', phi', psi),
// so in a Safety phase each Hardening Soil or Soft Soil material becomes the Mohr-Coulomb material
// with that strength and nu_ur, with the stiffness it has in the state the phase starts from --
// one Young's modulus per material, from the mean stress of its stress points there (E50 at the
// mean minor principal stress for Hardening Soil, 3(1 - 2 nu_ur) p' / kappa* for Soft Soil). The
// stiffness does not decide a limit state; the strength does, and that is the soil's own.
inline std::vector<MaterialModel> safety_equivalents(const katai::mesh::Mesh& mesh,
                                                     const std::vector<MaterialModel>& models,
                                                     const std::vector<GaussState>& state,
                                                     const std::vector<char>& active,
                                                     std::vector<std::string>* converted = nullptr) {
    std::vector<MaterialModel> out = models;
    const int ne = mesh.element_count;
    const int ng = ne > 0 ? (int)(state.size() / (size_t)ne) : 0;
    for (size_t mi = 0; mi < models.size(); ++mi) {
        const MaterialModel& m = models[mi];
        const bool hs = m.type == MaterialType::HardeningSoil;
        const bool ss = m.type == MaterialType::SoftSoil || m.type == MaterialType::SoftSoilCreep;
        if (!hs && !ss) continue;
        // Mean minor principal and mean effective stress (compression positive) over the material.
        double s3 = 0.0, pm = 0.0;
        int n = 0;
        for (int e = 0; e < ne && ng > 0; ++e) {
            if (mesh.element_material[e] != (int)mi || (!active.empty() && !active[e])) continue;
            for (int g = 0; g < ng; ++g) {
                const GaussState& gs = state[(size_t)e * ng + g];
                const double sx = -gs.stress(0), sy = -gs.stress(1), txy = -gs.stress(2);
                const double sz = -gs.stress_zz;
                const double c = 0.5 * (sx + sy), r = std::hypot(0.5 * (sx - sy), txy);
                s3 += std::min({c - r, sz});
                pm += (sx + sy + sz) / 3.0;
                ++n;
            }
        }
        if (n > 0) { s3 /= n; pm /= n; }
        MaterialModel mc = m;
        mc.type = MaterialType::MohrCoulomb;
        if (hs) {
            const auto& h = m.hs;
            const double sphi = std::sin(h.friction), cphi = std::cos(h.friction);
            const double num = h.cohesion * cphi + std::max(s3, h.p_limit()) * sphi;
            const double den = h.cohesion * cphi + h.p_ref * sphi;
            mc.youngs_modulus = h.E50_ref * (den > 0.0 ? std::pow(num / den, h.m) : 1.0);
            mc.poisson_ratio = h.nu_ur;
            mc.cohesion = h.cohesion;
            mc.friction_angle = h.friction;
            mc.dilatancy_angle = h.dilatancy;
        } else {
            const bool creep = m.type == MaterialType::SoftSoilCreep;
            const double kap = creep ? m.ssc.kap_star : m.ssoil.kap_star;
            const double nu = creep ? m.ssc.nu_ur : m.ssoil.nu_ur;
            mc.youngs_modulus = 3.0 * (1.0 - 2.0 * nu) * std::max(pm, 1.0) / std::max(kap, 1e-9);
            mc.poisson_ratio = nu;
            mc.cohesion = creep ? m.ssc.c : m.ssoil.c;
            mc.friction_angle = creep ? m.ssc.phi : m.ssoil.phi;
            mc.dilatancy_angle = creep ? m.ssc.psi : m.ssoil.psi;
        }
        out[mi] = mc;
        if (converted) converted->push_back(std::to_string(mi));
    }
    return out;
}


// Solve the phase. Fills the factor of safety, the honest lower-bound flag,
// the mechanism displacement and the recovered nodal stresses in R. Returns
// false on an honest refusal or an unstable-at-minimum outcome, with
// R.message set; on success the caller's common result tail completes R.
//
// `structures` are the phase's active structural elements and `f` must already carry their
// self-weight: both enter every trial of the search (safety_analysis states what the reduction
// does and does not touch). With in.from_parent() the search starts from the parent phase's
// equilibrium (safety_analysis states how); otherwise it re-solves the ground from the unstressed
// state. No structural force diagram is produced for this phase either way -- the state the
// search stops at is the ground at the limit of its reduced strength, not a state to design for.
inline bool solve_safety_phase(
    const katai::mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<MaterialProfile>& profiles,
    const Eigen::VectorXd& f, const LinearSolve& solver, const Structures& structures,
    const SafetyPhase& in, SolveResult& R) {
    if (!in.nonlinear_soil && !in.has_hardening) {
        R.message = "Safety analysis (phi-c reduction) needs a Mohr-Coulomb or Hardening Soil "
                    "model -- linear-elastic soil has no shear strength to reduce.";
        return false;
    }
    if (in.axisymmetric) {
        R.message = "Safety analysis is not available in axisymmetric mode yet.";
        return false;
    }
    // Hardening Soil / Soft Soil Safety is gated honestly: phi-c reduction re-solves gravity per
    // trial, and a cap model from a stress-free state over-shoots the cap -> it either diverges
    // (hang) or returns a FALSE collapse (a meaningless factor of safety, even for a stable
    // confined block). factor_strength does reduce the soft-soil strength too (defensive), but
    // until a path-stable Safety -- strength reduction from the geostatic equilibrium, a tracked
    // follow-up -- the honest refusal stands: ask for Mohr-Coulomb strength for the Safety check.
    // From the parent phase's equilibrium the reduction is path-stable for these models too: the
    // cap is where the phases left it and nothing is re-loaded from zero, so the gate below is the
    // unstressed start's alone.
    if ((in.has_hardening || in.has_softsoil) && !in.from_parent()) {
        R.message = "Safety analysis (phi-c reduction) for Hardening Soil / Soft Soil cannot "
                    "start from the unstressed state (reducing strength from a stress-free state is path-unstable "
                    "with a cap model and would report a misleading factor of safety). Run the "
                    "Safety analysis as a phase after the initial phase (or after the phase that "
                    "builds the ground), and the search starts from the state it built; or use a "
                    "Mohr-Coulomb material with the same c' and phi'.";
        return false;
    }
    // HOEK-BROWN HAS NOTHING FOR phi-c REDUCTION TO HOLD, and the failure is silent in the
    // dangerous direction. factor_strength divides c and phi -- fields this model never reads --
    // so every trial would run the rock at FULL strength, no reduction would ever bring it down,
    // and the search would walk to its cap and report the cap itself -- "FoS > 3.0", for any
    // rock mass whatever, however weak. Nothing here would say so.
    //
    // Strength reduction IS defined for this model, and it is not a conversion: the Hoek-Brown
    // yield function is rewritten with the strength reduction factor inside it, evaluated at
    // each stress point's own sigma'_3 (Benz, Schwab, Kauther & Vermeer 2008, Int. J. Rock
    // Mech. Min. Sci. 45(2), 210-222). The factor that gives is in general not the factor of
    // safety of a Mohr-Coulomb material fitted to the same envelope, because the two strength
    // curves are reduced differently. This build does not implement that form, so the
    // analysis is refused. (Until 2026-09 this comment said no rule exists. It had been drawn from
    // a description that defines only the tension cut-off's reduction; the rule for the yield
    // function itself is the reformulation above.)
    //
    // What the refusal offers instead is a Mohr-Coulomb fit, and it has to say what that is. The
    // conversion (the equivalent c'/phi' of Hoek, Carranza-Torres & Corkum 2002, a linear fit that
    // balances the areas above and below the curve) is a fit over a confining range whose upper
    // limit sigma'_3max is left to the engineer, because it depends on the application -- and, as
    // stated above, its factor of safety is not the Hoek-Brown one. It is the engineer's
    // approximation, stated with the range it was made over, not a stand-in for the rule.
    for (const auto& m : models)
        if (m.type == MaterialType::HoekBrown) {
            R.message =
                "Safety analysis (phi-c reduction) is not available for the Hoek-Brown model in this "
                "build: its strength is the Hoek-Brown curve, not a c'/phi' pair, so reducing c' and "
                "phi' would leave the rock at full strength through every trial and the factor of "
                "safety reported would be too HIGH, not merely imprecise. The strength reduction "
                "defined for this model reformulates the Hoek-Brown yield function itself, and "
                "this build does not implement it yet. A Mohr-Coulomb material fitted to the "
                "envelope over the confining range your problem actually spans will run a Safety "
                "phase, but the factor it gives belongs to that fit: it does not correspond to the "
                "Hoek-Brown factor of safety, so state it as an approximation, with the range it "
                "was made over.";
            return false;
        }
    if (in.from_parent() && (in.has_hardening || in.has_softsoil))
        add_diagnostic(R, DiagnosticSeverity::Note, "K2D-A020", "Safety",
                       "The Hardening Soil / Soft Soil materials take part in this Safety phase "
                       "as Mohr-Coulomb materials with their own c', phi', psi and nu_ur, and a "
                       "stiffness taken from the stress state the phase starts from. Their "
                       "hardening is not reduced with the strength: reduced together, the factor "
                       "of safety depended on the path of the search.");
    StrengthReductionOptions sopt;
    sopt.srf_min = 0.4; sopt.srf_max = 3.0; sopt.bisection_iterations = 12;
    sopt.newton = NewtonOptions{in.load_steps > 0 ? in.load_steps : 8,
                                in.max_iterations > 0 ? in.max_iterations : 120,
                                in.tolerance > 0.0 ? in.tolerance : 1e-3};
    const auto sr =
        in.from_parent()
            ? safety_analysis(mesh, dofs, f,
                              safety_equivalents(mesh, models, in.parent_state, in.active),
                              solver, sopt, in.parent_state, in.active, profiles, structures,
                              in.parent_force, in.parent_structures)
            : safety_analysis(mesh, dofs, f, models, solver, sopt, {}, in.active, profiles,
                              structures);
    R.fos = sr.fos;
    if (sr.start_failed) {
        R.message = "Safety analysis: the ground did not reach equilibrium at full strength under "
                    "this phase's own changes to the phase before it (activated or removed "
                    "elements, loads or water), so no factor of safety is defined. Make those "
                    "changes in a Plastic phase first and run the Safety phase after it.";
        return false;
    }
    if (!sr.ok) {
        R.message = "Safety analysis: the slope did not reach equilibrium even at the lowest "
                    "strength factor (it may already be unstable, or under-restrained).";
        return false;
    }
    // No collapse anywhere up to srf_max: the FoS is only a lower bound (the cap). Common when
    // the geometry is not a slope (laterally-confined block under self-weight never fails in
    // shear), the slope face/crest are not Free, or there is no destabilizing load. Report it
    // honestly (the GUI shows "FoS > cap") instead of presenting the cap as a real FoS.
    R.fos_lower_bound = !sr.bracketed;
    R.disp = sr.mechanism.displacement.head(mesh.node_count * 2);
    R.stress = recover_nodal_stresses_from_gauss(mesh, sr.mechanism.gauss_states, in.active);
    return true;
}

} // namespace katai::core
