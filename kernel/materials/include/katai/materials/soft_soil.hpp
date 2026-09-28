#pragma once
// SOFT SOIL model — STAGE 1: the material-point core (locked formulation
// docs/references/soft-soil-formulation.md). Cam-Clay type: ln-law compression
// (λ*/κ*), ellipse cap (q̃, the same deviatoric measure as the HS cap) + exponential p_p
// hardening (associated), Mohr-Coulomb failure (M is NOT critical-state; derived from K0NC,
// Brinkgreve 1994).
//
// SCOPE (honest): principal-space, strain-driven, COAXIAL core; the MC return is the FLAT
// surface (triaxial compression path — edges/apex with the full MC machinery later); D_e is
// frozen at the trial p within the return (an integration error that vanishes as the step
// shrinks; the V&V bands account for it). The FE integration (Stage 2) lives in
// material_model.hpp: ss_return_core builds the Voigt trial with the exponential elastic
// law defined here, decomposes to principals and wraps this core with the EXACT INVERSE of
// the strain increment (Δε_v = κ*·ln(p_tr/p_n), Δe = Δs/2G) — elastic steps reproduce the
// trial to round-off. Initial p_p (OCR/POP) + GUI/project file = Stage 3. Sign convention
// IN THIS HEADER is compression-POSITIVE (like the HS core).

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

namespace katai::core::softsoil {

// p' never drops below unit stress (1 kPa) — a floor against the ln-law/K=p'/κ*
// singularity. ss_step and the FE wrapper (material_model.hpp ss_return_core) must use the
// SAME floor so the elastic-predictor invertibility (exact reconstruction of the trial)
// is not broken.
inline constexpr double kPmin = 1.0;

struct Params {
    double lam_star = 0.1;   // modified compression index λ* [-]
    double kap_star = 0.02;  // modified swelling index κ* [-]
    double nu_ur = 0.15;     // unloading/reloading Poisson ratio
    double c = 0.0;          // effective cohesion [kPa]
    double phi = 0.0;        // effective friction angle [rad] (0 FORBIDDEN: c·cotφ undefined)
    double psi = 0.0;        // dilatancy [rad] (SS default 0)
    double K0nc = 0.5;       // normal-consolidation lateral pressure coefficient (→ M)
};

// M(K0NC) — Brinkgreve (1994), the closed form written out below. Behaviour pin: oedometer
// primary loading must produce σ'_h/σ'_v → K0NC (test_soft_soil (e) measures this DIRECTLY;
// a transcription error in the formula blows up there).
inline double M_from_K0nc(const Params& P) {
    const double K = P.K0nc, nu = P.nu_ur, R = P.lam_star / P.kap_star;
    const double a = (1.0 - K) * (1.0 - K) / ((1.0 + 2.0 * K) * (1.0 + 2.0 * K));
    const double b = (1.0 - K) * (1.0 - 2.0 * nu) * (R - 1.0) /
                     ((1.0 + 2.0 * K) * (1.0 - 2.0 * nu) * R - (1.0 - K) * (1.0 + nu));
    return 3.0 * std::sqrt(a + b);
}

struct StepResult {
    Eigen::Vector3d sig;     // principal stresses (compression-positive)
    double pp;               // new preconsolidation pressure
    bool cap_active = false;
    bool mc_active = false;
    int nsub = 1;            // number of substeps used (reported for FD tangent pinning)
};

namespace detail {
// q̃ = σ1 + (δ−1)σ2 − δσ3, sorted σ1 ≥ σ2 ≥ σ3 (compression-positive), δ = (3+sinφ)/(3−sinφ)
// (the HS cap measure).
inline double q_tilde(const Eigen::Vector3d& s, double delta) {
    Eigen::Vector3d o = s;
    std::sort(o.data(), o.data() + 3, std::greater<double>());
    return o(0) + (delta - 1.0) * o(1) - delta * o(2);
}
}  // namespace detail

// ONE substep (internal — call ss_step): explicit return with the flow frozen at the trial.
// On a large increment the elastic trial overshoots the surface far and the frozen flow
// direction misses the target manifold (measured: at a Δε_v/κ* = 0.5 step the cap secant
// found no root and was silently skipped) — hence ss_step splits the increment into
// substeps via ss_nsub (HS's "unconditionally stable" substepping pattern).
StepResult ss_substep(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                      const Eigen::Vector3d& deps);

// FE initial preconsolidation (K0 seeding; the parallel of hs_initial_pp): pp = f̄(σ0)·OCR_eq
// — f̄ is the equivalent pressure at which the cap passes through σ0 (the equivalent
// isotropic pressure p_eq). NC (OCR_eq=1) ⇒ the state sits EXACTLY on the cap (f=0,
// admissible; a pp=0 seed would make the cap yield from the start). OCR_eq: OCR mode is the
// ratio directly; POP mode is converted at the caller to the equivalent ratio
// (σ'_v0+POP)/σ'_v0 (f̄ is NOT first-order homogeneous in σ (the c·cotφ shift), but the ratio
// scaling is the consistent counterpart of OCR and POP being defined on the vertical stress).
// Floor: pp ≥ max(c·cotφ, unit stress) — the threshold ellipse.
inline double ss_initial_pp(const Params& P, const Eigen::Vector3d& sig_comp_pos,
                            double ocr_eq = 1.0) {
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi);
    const double ccot = P.phi > 1e-12 ? P.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = M_from_K0nc(P);
    const double p = sig_comp_pos.mean();
    const double qt = detail::q_tilde(sig_comp_pos, delta);
    const double fbar = qt * qt / (M * M * std::max(p + ccot, 1e-9)) + p;
    return std::max(fbar * std::max(ocr_eq, 1.0), std::max(ccot, kPmin));
}

// Automatic substep count: the per-substep stress-ratio measure (|Δε_v| + 2·max|Δε_dev|)/κ*
// is kept ≤ 0.05 — since both the exponential mean (exp(Δε_v/κ*)) and the deviatoric trial
// (2G ∝ p/κ*) grow on the same 1/κ* scale, one criterion bounds both. Measured: at a 0.5
// measure the return misses, stable at ≤0.1; 0.05 with a safety margin.
inline int ss_nsub(const Params& P, const Eigen::Vector3d& deps) {
    const double em = deps.sum() / 3.0;
    const double dev = (deps - Eigen::Vector3d::Constant(em)).cwiseAbs().maxCoeff();
    return 1 + static_cast<int>((std::fabs(deps.sum()) + 2.0 * dev) / (0.05 * P.kap_star));
}

// Apply a strain increment (principal, coaxial): sig_c committed stress
// (compression-positive), pp_c committed preconsolidation, deps principal strain increment
// (positive = compression). nsub_fixed > 0 pins the substep count: the FD tangent's
// perturbed runs must follow the base run's substep sequence, otherwise the ceil jump
// pollutes the FD column by the integration error (a lesson measured on HS;
// material_model.hpp integrate_point SS branch).
StepResult ss_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                   const Eigen::Vector3d& deps, int nsub_fixed = 0);

}  // namespace katai::core::softsoil
