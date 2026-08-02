#include <katai/analysis/strength_reduction.hpp>

#include <cmath>
#include <vector>

namespace katai::core {

// Factor a material's SHEAR STRENGTH by the strength reduction factor srf (phi-c reduction):
//   c_f = c / srf,   phi_f = atan(tan phi / srf)
// applied to BOTH the Mohr-Coulomb fields (used by the MC model) AND the Hardening Soil sub-struct
// (the HS failure surface uses hs.cohesion / hs.friction -- a SEPARATE set of fields, so reducing
// only the MC ones leaves an HS slope at full strength and it never collapses). Dilatancy is clamped
// to the reduced friction afterwards (psi <= phi is required; PLAXIS practice) so the reduced state
// stays physically admissible at high srf.
static void factor_strength(MaterialModel& m, double srf) {
    m.cohesion /= srf;
    m.friction_angle = std::atan(std::tan(m.friction_angle) / srf);
    if (m.dilatancy_angle > m.friction_angle) m.dilatancy_angle = m.friction_angle;
    // Tensile strength is a material strength too: PLAXIS documents the Safety
    // reduction of the tension cut-off value explicitly for Hoek-Brown (MMM sec
    // 4.3.7); the same rule is applied to the MC cap here (safe direction; the
    // c*cot(phi) clamp inside the return mapping shrinks with c/srf as well).
    m.tensile_strength /= srf;
    m.hs.cohesion /= srf;
    m.hs.friction = std::atan(std::tan(m.hs.friction) / srf);
    if (m.hs.dilatancy > m.hs.friction) m.hs.dilatancy = m.hs.friction;
    // Soft Soil carries its own c/phi fields (ssoil) -- without reducing them Safety would
    // SILENTLY report an inflated FoS. (SS Safety is currently gated off in build_problem;
    // these lines are defensive so it stays correct if the gate opens. The K0nc/M
    // calibration is stiffness structure, not strength: constant.)
    m.ssoil.c /= srf;
    m.ssoil.phi = std::atan(std::tan(m.ssoil.phi) / srf);
    if (m.ssoil.psi > m.ssoil.phi) m.ssoil.psi = m.ssoil.phi;
}

double factor_of_safety(const mesh::Mesh& mesh, const DofMap& dofs,
                        const Eigen::VectorXd& gravity_load,
                        const MaterialModel& base, const LinearSolve& linear_solve,
                        const StrengthReductionOptions& options) {
    // Whether the slope reaches equilibrium under gravity with strength factored
    // by srf. Re-solved from the unstressed state each trial.
    auto is_stable = [&](double srf) {
        MaterialModel m = base;
        factor_strength(m, srf);
        const std::vector<MaterialModel> materials = {m};
        try {
            const NewtonResult r = solve_nonlinear(
                mesh, dofs, materials, gravity_load, linear_solve, options.newton);
            return r.converged;
        } catch (...) {
            // A singular/!factorizable tangent at incipient collapse is itself
            // the failure signal -- treat it as an unstable (collapsed) slope.
            return false;
        }
    };

    double lo = options.srf_min, hi = options.srf_max;
    for (int i = 0; i < options.bisection_iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (is_stable(mid))
            lo = mid;  // still stable -> factor of safety is higher
        else
            hi = mid;  // collapsed -> factor of safety is lower
    }
    return 0.5 * (lo + hi);
}

SafetyResult safety_analysis(const mesh::Mesh& mesh, const DofMap& dofs,
                             const Eigen::VectorXd& gravity_load,
                             const std::vector<MaterialModel>& materials,
                             const LinearSolve& linear_solve,
                             const StrengthReductionOptions& options,
                             const std::vector<GaussState>& initial_state,
                             const std::vector<char>& active_element,
                             const std::vector<MaterialProfile>& profile) {
    // Factor every material's strength by srf and re-solve under gravity from the unstressed state;
    // failure to converge is the standard collapse signal. (This is robust for Mohr-Coulomb. Hardening
    // Soil from a stress-free state over-shoots the cap and is unreliable here -- HS Safety is gated
    // upstream until a path-stable scheme is in place; see build_problem.)
    auto try_srf = [&](double srf, NewtonResult& out) -> bool {
        std::vector<MaterialModel> m = materials;
        for (auto& mm : m) factor_strength(mm, srf);
        // The COHESION GRADIENT must be reduced by srf too. c(y) = c_ref + c_inc (y_ref - y) is one
        // strength; factoring only c_ref would leave the deep soil at its full gradient strength and
        // report an FoS that is too HIGH -- unconservative, and invisible. (E_inc is a stiffness, not a
        // strength: phi-c reduction never touches it, so it rides through unchanged.)
        std::vector<MaterialProfile> p = profile;
        for (auto& pp : p) pp.c_inc /= srf;
        try {
            out = solve_nonlinear(mesh, dofs, m, gravity_load, linear_solve, options.newton,
                                  initial_state, active_element, {}, {}, {}, p);
            return out.converged;
        } catch (...) {
            return false;  // singular tangent at incipient collapse = the failure signal
        }
    };

    double lo = options.srf_min, hi = options.srf_max;
    SafetyResult res;
    NewtonResult trial;
    for (int i = 0; i < options.bisection_iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (try_srf(mid, trial)) { lo = mid; res.mechanism = trial; res.ok = true; }
        else { hi = mid; res.bracketed = true; }  // a collapse was seen -> the FoS is finite/bracketed
    }
    res.fos = 0.5 * (lo + hi);
    return res;
}

} // namespace katai::core
