#pragma once
// SOFT SOIL CREEP (SSC) — STAGE 1: the material-point core (Vermeer & Neher 1999; locked
// formulation docs/references/soft-soil-creep-formulation.md). Time-dependent
// viscoplasticity: creep = time-dependent plastic strain, potential g = p_eq (the MCC
// ellipse, IDENTICAL to the SS f̄ measure), volumetric rate ε̇_v^c = (μ*/τ)(p_eq/p_p)^β,
// β = (λ*−κ*)/μ*, p_p ages exponentially; failure is a SEPARATE Mohr-Coulomb check AFTER
// the creep update (creep first, then the failure check, at every stress point).
//
// SINGLE-SOURCE reuse: the elastic law/M(K0NC)/q̃/kPmin come from softsoil; the MC return is
// done by giving softsoil::ss_step pp=∞ (cap off) — same elastic moduli, same edge-cascaded
// MC. The creep direction is also the sorted-role analytic of the SS cap normal (a·w + b·1)
// and is applied with the ORDERING-VALIDITY cascade (Koiter; the corner-branching lesson
// measured on SS).
//
// SCOPE (Stage 1, honest): principal-space, coaxial; the direction is FROZEN at the trial
// within a substep (explicit direction / implicit magnitude: L backward Euler, monotone
// R → bracketed secant+bisection guaranteed); on the dry side of the ellipse
// ∂p_eq/∂p′ < 0 → dilative creep (the model's own truth; the |·| floor is numerical
// protection); the deviatoric creep direction is not pinned by the Stage-1 closed forms
// (the FE stage does it, with constant-rate-of-strain and undrained creep cases). Sign
// convention compression-POSITIVE. FE/GUI = the next stage.

#include <katai/materials/soft_soil.hpp>

namespace katai::core::softsoilcreep {

struct Params {
    double lam_star = 0.10;   // modified compression index λ* [-]
    double kap_star = 0.02;   // modified swelling index κ* [-]
    double mu_star = 0.005;   // modified creep index μ* [-] (= Cα/(2.3(1+e)); λ*/μ* typically 15-25)
    double nu_ur = 0.15;
    double c = 0.0;           // [kPa]
    double phi = 0.0;         // [rad]
    double psi = 0.0;         // [rad]
    double K0nc = 0.5;        // → M (Brinkgreve 1994; same as SS)
    double tau_day = 1.0;     // reference time τ [days] — the 24-hour definition of the NC line

    softsoil::Params ss() const {   // shared-machinery view (elastic law, M, MC, q̃)
        softsoil::Params S;
        S.lam_star = lam_star; S.kap_star = kap_star; S.nu_ur = nu_ur;
        S.c = c; S.phi = phi; S.psi = psi; S.K0nc = K0nc;
        return S;
    }
};

struct StepResult {
    Eigen::Vector3d sig;      // principal stresses (compression-positive)
    double pp;                // current preconsolidation (in the p_eq measure)
    double devc = 0.0;        // volumetric creep strain accumulated this step
    bool mc_active = false;
    int nsub = 1;
};

// p_eq = p′ + q̃²/(M²(p′+c·cotφ)) — the measure identical to the SS cap function f̄
// (by the single-source reading the formula matches ss_initial_pp/f_cap).
inline double p_eq(const softsoil::Params& S, const Eigen::Vector3d& sig) {
    const double sphi = std::sin(S.phi), cphi = std::cos(S.phi);
    const double ccot = S.phi > 1e-12 ? S.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = softsoil::M_from_K0nc(S);
    const double p = sig.mean();
    const double qt = softsoil::detail::q_tilde(sig, delta);
    return qt * qt / (M * M * std::max(p + ccot, 1e-9)) + p;
}

// ONE time/strain substep (internal — call ssc_step).
StepResult ssc_substep(const Params& P, const softsoil::Params& S,
                       const Eigen::Vector3d& sig_c, double pp_c,
                       const Eigen::Vector3d& deps, double dt);

// Automatic substepping: the SS strain criterion + a creep-magnitude criterion (keep the
// explicit estimate L_exp ≤ 0.05·κ* per substep — this bounds both the direction freezing
// and the backward-Euler time error). Upper limit 500 (the implicit L solve is
// unconditionally stable; on extremely fast transitions the accuracy bands in the tests
// are honest).
inline int ssc_nsub(const Params& P, const softsoil::Params& S, const Eigen::Vector3d& sig_c,
                    double pp_c, const Eigen::Vector3d& deps, double dt) {
    const double beta = (P.lam_star - P.kap_star) / P.mu_star;
    const double pe = std::max(p_eq(S, sig_c), 1e-12);
    const double ratio = pe / std::max(pp_c, softsoil::kPmin);
    const double rate = (P.mu_star / P.tau_day) * std::pow(std::min(ratio, 10.0), beta);
    const double Lexp = dt * rate;
    const int n_creep = 1 + static_cast<int>(Lexp / (0.05 * P.kap_star));
    const int n_strain = softsoil::ss_nsub(S, deps);
    return std::min(std::max(n_creep, n_strain), 500);
}

// Apply a strain + TIME increment (principal, coaxial): deps total increment (positive =
// compression), dt in days. nsub_fixed > 0 → substeps pinned (for FD tangent columns; the
// HS/SS lesson).
StepResult ssc_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                    const Eigen::Vector3d& deps, double dt, int nsub_fixed = 0);

}  // namespace katai::core::softsoilcreep
