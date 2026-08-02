#pragma once
// Safety phase strategy (Stage B9). Phi-c reduction / strength reduction
// method, PLAXIS "Safety": bisect the strength reduction factor to the slope's
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
};

// Solve the phase. Fills the factor of safety, the honest lower-bound flag,
// the mechanism displacement and the recovered nodal stresses in R. Returns
// false on an honest refusal or an unstable-at-minimum outcome, with
// R.message set; on success the caller's common result tail completes R.
inline bool solve_safety_phase(
    const katai::mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<MaterialProfile>& profiles,
    const Eigen::VectorXd& f, const LinearSolve& solver,
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
    if (in.has_hardening || in.has_softsoil) {
        R.message = "Safety analysis (phi-c reduction) for Hardening Soil / Soft Soil is not "
                    "supported yet (reducing strength from a stress-free state is path-unstable "
                    "with a cap model and would report a misleading factor of safety). For a "
                    "slope factor of safety, use a Mohr-Coulomb material with the same c' and "
                    "phi' -- the standard strength-reduction model.";
        return false;
    }
    StrengthReductionOptions sopt;
    sopt.srf_min = 0.4; sopt.srf_max = 3.0; sopt.bisection_iterations = 12;
    sopt.newton = NewtonOptions{8, 120, 1e-3};
    const auto sr = safety_analysis(mesh, dofs, f, models, solver, sopt, {}, in.active, profiles);
    R.fos = sr.fos;
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
