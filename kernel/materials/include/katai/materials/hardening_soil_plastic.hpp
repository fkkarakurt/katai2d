#pragma once
// Hardening Soil — P2.3d (part 1): the multiaxial SHEAR hardening return mapping.
// Strain-driven predictor-corrector (FE constitutive form) in principal stress space.
// Elastic: stress-dependent Eur (frozen at σ3 at the start of the step). Shear yield
// surface (Schanz 1999):
//   f = f̄(q) − γ^p ,   f̄(q) = (2/Ei) q/(1−q/qa) − 2q/Eur ,   q = σ1 − σ3
// Non-associated flow: mobilized dilatancy ψ_m, m_g = (1, R, R), R = ε3^p/ε1^p. ψ_m is Rowe's
// where Rowe is non-negative; below the phase-transformation line HS takes zero and HSsmall
// takes Li & Dafalias instead — one rule, `detail::hs_dilatancy` (see its comment).
//   Rowe flow rule: ε_v^p/ε_q^p = −sinψ_m (DILATION, comp-pos ⇒ ε_v^p<0) ⇒
//   R = −(1+sinψ_m)/(2−sinψ_m) ; hardening γ^p=−(2ε1^p−ε_v^p) ⇒ dγ^p=h·dλ,
//   h = 1−2R = (4+sinψ_m)/(2−sinψ_m). At ψ_m=0, R=−½, h=2 (volumetrically neutral, the
//   hyperbola is preserved).
// Consistent (asymmetric) tangent D^ep = D_e − (D_e m_g)(n_f^T D_e)/(n_f^T D_e m_g + h).
// In the ψ=0 check the elastoplastic machinery PRODUCES the drained triaxial hyperbola
// (not closed-form) → the machinery's correctness (test_hs_shear). Math:
// hardening-soil-formulation.md.
//
// CONVENTION: compression-POSITIVE principal stresses σ1 ≥ σ2 ≥ σ3 (sigma[0..2]). The sign
// flips in FE (σ_HS = −σ_solver) — that integration is wired in a later part.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <Eigen/Dense>

#include <katai/materials/hardening_soil.hpp>

namespace katai::core {

struct HsShearStep {
    Eigen::Vector3d stress;   // updated principal stresses
    double gamma_p;           // updated hardening parameter
    Eigen::Matrix3d tangent;  // consistent D^ep (principal space 3×3)
    bool plastic;
};

namespace detail {

// Isotropic elastic principal-space stiffness (E, ν).
inline Eigen::Matrix3d hs_elastic(double E, double nu) {
    const double f = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    Eigen::Matrix3d d;
    d << f * (1 - nu), f * nu, f * nu,
         f * nu, f * (1 - nu), f * nu,
         f * nu, f * nu, f * (1 - nu);
    return d;
}

// Mobilized dilatancy ψ_m as a function of the mobilized friction angle φ_m. This rule has
// FOUR callers (the two principal returns and their two cap-coupled siblings) and used to be
// written out four times; a rule that is right in three copies and stale in the fourth is a
// silently wrong answer on whichever stress path reaches the fourth. One definition, four
// bindings.
//
// HS (Rowe 1962; Schanz & Vermeer 1996):  sinψ_m = (sinφ_m − sinφ_cv)/(1 − sinφ_m·sinφ_cv),
// cut off to [0, sinψ]. Where Rowe returns a NEGATIVE value the Hardening Soil model takes ψ_m = 0.
//
// HSsmall (after Li & Dafalias (2000)): that zero cut-off can give too little plastic
// volumetric strain, so wherever Rowe is negative the small-strain model puts a small
// CONTRACTION there instead of nothing:
//   sinψ_m = (1/10)·(−M_c·exp[(1/15)·ln((M_c/M_d)·(q/q_a))] + M_d)          (i)
//   M_c    = 6·sinφ_cv/(3 − sinφ_cv)                                        (ii)
//   M_d    = 6·sinφ_m /(3 − sinφ_m )                                        (iii)
//   q/q_a  = max([(1−sinφ_cv)/sinφ_cv]·[sinφ_m/(1−sinφ_m)], 1e−4)           (iv)
//   sinφ_m ≥ sinφ/(2 − sinφ)                                                (v)
// The q/q_a of (iv) is its own definition, NOT the model's deviatoric hyperbola ratio q/q_a
// (q_a = q_f/R_f); the two are unrelated quantities that share a name.
//
// The two branches MEET EXACTLY, and that is the property worth knowing: at φ_m = φ_cv we get
// M_d = M_c and q/q_a = 1, so the logarithm vanishes and (i) returns exactly 0 — which is
// also precisely where Rowe changes sign. The composite rule is continuous by CONSTRUCTION,
// not to a tolerance. The floor (v) then bounds how contractant it can get, since sinψ_m falls
// monotonically as φ_m drops: the floor is the most negative ψ_m the model can produce.
// (test_hssmall pins both.)
//
// WORKED VALUES (φ=35°, ψ=5°): the upper branch is pure Rowe, the curve vanishes at
// φ_cv = 30.80°, and the plateau starts at the floor (v), φ_m = 23.71°, where (i) gives
// ψ_m = −1.679°. We implement the rule exactly as written above, with the leading 1/10. A
// plotted version of this rule whose amplitude is 1.29× larger (plateau −2.167°, as if the
// leading 1/10 were 1/7.74) is declared as a gap in hssmall-formulation.md; both readings are
// small contractions and the difference between them is far smaller than the difference
// between either and the ψ_m = 0 this replaces.
struct HsDilatancy {
    double sin_cs = 0.0;        // sinφ_cv (Rowe, from the input φ and ψ)
    double sin_psi = 0.0;       // sinψ — the upper cut-off
    double c_cot = 0.0;         // c·cotφ — the cohesion shift of the mobilized-φ definition
    double Mc = 0.0;            // M_c, (ii) (HSsmall branch only)
    double sphi_m_floor = 0.0;  // floor (v) (HSsmall branch only)
    bool li_dafalias = false;   // HSsmall (G0_ref>0) with a usable φ_cv

    double operator()(double sphi_m) const {
        const double rowe = (sphi_m - sin_cs) / (1.0 - sphi_m * sin_cs);
        if (rowe >= 0.0) return std::min(rowe, sin_psi);
        if (!li_dafalias) return 0.0;                                   // HS: the zero cut-off
        const double s = std::max(sphi_m, sphi_m_floor);                // floor (v)
        const double Md = 6.0 * s / (3.0 - s);                          // M_d, (iii)
        const double qqa =
            std::max((1.0 - sin_cs) / sin_cs * (s / (1.0 - s)), 1e-4);  // q/q_a, (iv)
        return 0.1 * (-Mc * std::exp(std::log((Mc / Md) * qqa) / 15.0) + Md);  // sinψ_m, (i)
    }

    // sinφ_m = q/(σ1 + σ3 + 2c·cotφ) from the mobilized stress state (σ1 = q + σ3, comp-positive).
    double from_q(double q, double s3v) const {
        const double s1 = q + s3v;
        const double denom = s1 + s3v + 2.0 * c_cot;
        return (*this)(denom > 1e-12 ? q / denom : 0.0);
    }
};

inline HsDilatancy hs_dilatancy(const HardeningSoilParams& p) {
    HsDilatancy d;
    const double sphi = std::sin(p.friction), cphi = std::cos(p.friction);
    const double sps = std::sin(p.dilatancy);
    d.sin_cs = (sphi - sps) / (1.0 - sphi * sps);  // critical state
    d.sin_psi = sps > 0.0 ? sps : 0.0;
    d.c_cot = (sphi > 1e-12) ? p.cohesion * cphi / sphi : 0.0;
    // φ_cv = 0 (a φ = 0 Tresca soil) would divide by zero in (iv) — but it also makes Rowe
    // non-negative everywhere, so the branch is unreachable there; the guard says so rather
    // than relying on it.
    d.li_dafalias = p.G0_ref > 0.0 && !p.dilatancy_cut && d.sin_cs > 1e-12 && sphi < 1.0;
    if (d.li_dafalias) {
        d.Mc = 6.0 * d.sin_cs / (3.0 - d.sin_cs);  // M_c, (ii)
        d.sphi_m_floor = sphi / (2.0 - sphi);      // floor (v)
    }
    return d;
}

}  // namespace detail

// One shear-hardening step: committed (σ_n, γ_n) + strain increment dε → updated state +
// consistent tangent. Eur/Ei/qa are frozen at the start-of-step σ3 = σ_n[2] (explicit
// stress-dependent stiffness). q is capped at qf at failure (perfectly
// plastic MC).
inline HsShearStep hs_shear_step(const HardeningSoilParams& p,
                                 const Eigen::Vector3d& sigma_n, double gamma_p_n,
                                 const Eigen::Vector3d& deps) {
    const double s3 = sigma_n(2);
    const double Eur = p.Eur(s3), nu = p.nu_ur;
    const double Ei = p.Ei(s3), qa = p.q_asymptote(s3), qf = p.q_failure(s3);
    const Eigen::Matrix3d De = detail::hs_elastic(Eur, nu);

    const Eigen::Vector3d sig_tr = sigma_n + De * deps;
    const double q_tr = sig_tr(0) - sig_tr(2);

    auto fbar = [&](double q) {
        return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur;
    };
    auto fbar_prime = [&](double q) {
        const double r = 1.0 - q / qa;
        return (2.0 / Ei) / (r * r) - 2.0 / Eur;
    };

    HsShearStep out;
    const double f_tr = fbar(q_tr) - gamma_p_n;
    if (f_tr <= 1e-12 * (1.0 + std::fabs(gamma_p_n))) {  // elastic
        out.stress = sig_tr;
        out.gamma_p = gamma_p_n;
        out.tangent = De;
        out.plastic = false;
        return out;
    }

    // Mobilized dilatancy. φ_cs constant (from input φ, ψ); ψ_m varies with q.
    const detail::HsDilatancy dil = detail::hs_dilatancy(p);
    auto sin_psi_m = [&](double q, double s3v) { return dil.from_q(q, s3v); };

    // Local Newton: dλ such that f(σ(dλ), γ+h·dλ) = 0. σ(dλ) = σ_tr − dλ·D_e·m_g.
    double dlam = 0.0;
    Eigen::Vector3d mg(1.0, -0.5, -0.5);
    double h = 2.0;
    Eigen::Vector3d sig = sig_tr;
    double q = q_tr;
    for (int it = 0; it < 50; ++it) {
        const double spm = sin_psi_m(q, s3);
        const double R = -(1.0 + spm) / (2.0 - spm);
        mg = Eigen::Vector3d(1.0, R, R);
        h = (4.0 + spm) / (2.0 - spm);  // 1−2R (γ^p=−(2ε1^p−ε_v^p)); 2 at ψ_m=0
        const Eigen::Vector3d Demg = De * mg;
        sig = sig_tr - dlam * Demg;
        q = sig(0) - sig(2);
        const double resid = fbar(q) - (gamma_p_n + h * dlam);
        const double dq_dl = -(Demg(0) - Demg(2));
        const double dr_dl = fbar_prime(q) * dq_dl - h;
        const double step = resid / dr_dl;
        dlam -= step;
        if (std::fabs(step) <= 1e-14 * (1.0 + std::fabs(dlam))) break;
    }

    // MC failure bound: cap q at qf (perfectly plastic plateau).
    if (q > qf) {
        const double scale = (q - qf);
        // q = σ1 − σ3; drop the excess from σ1 to bring q to qf (constant-σ3 triaxial approximation).
        sig(0) -= scale;
        q = sig(0) - sig(2);
    }

    const double spm = sin_psi_m(q, s3);
    const double R = -(1.0 + spm) / (2.0 - spm);
    mg = Eigen::Vector3d(1.0, R, R);
    h = (4.0 + spm) / (2.0 - spm);  // 1−2R (γ^p hardening); 2 at ψ_m=0
    Eigen::Vector3d nf(fbar_prime(q), 0.0, -fbar_prime(q));  // ∂f/∂σ = f̄'(q)(1,0,−1)

    const Eigen::Vector3d Demg = De * mg;
    const Eigen::RowVector3d nfDe = nf.transpose() * De;
    const double denom = nfDe * mg + h;
    out.tangent = De - (Demg * nfDe) / denom;
    out.stress = sig;
    out.gamma_p = gamma_p_n + h * dlam;
    out.plastic = true;
    return out;
}

// Stress-driven SHEAR return (for the FE wrapping): given the elastic trial principal
// stress (compression-positive, σ1≥σ2≥σ3) + committed γ^p, returns the corrected principal
// stress + γ^p (NO tangent — numerical at the Voigt level). Stiffness frozen at σ3
// (sigma3_stiff = committed σ3). The trial-based sibling of hs_shear_step (the predictor
// lives outside).
struct HsPrincipalReturn {
    Eigen::Vector3d stress;  // corrected principal stress (compression-positive)
    double gamma_p;
};

inline HsPrincipalReturn hs_shear_correct(const HardeningSoilParams& p,
                                          const Eigen::Vector3d& sig_tr,
                                          double gamma_p_n, double sigma3_stiff) {
    const double Eur = p.Eur(sigma3_stiff), nu = p.nu_ur;
    const double Ei = p.Ei(sigma3_stiff), qa = p.q_asymptote(sigma3_stiff);
    const double qf = p.q_failure(sigma3_stiff);
    const Eigen::Matrix3d De = detail::hs_elastic(Eur, nu);
    const double q_tr = sig_tr(0) - sig_tr(2);

    auto fbar = [&](double q) { return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur; };
    auto fbar_prime = [&](double q) {
        const double r = 1.0 - q / qa;
        return (2.0 / Ei) / (r * r) - 2.0 / Eur;
    };

    HsPrincipalReturn out{sig_tr, gamma_p_n};
    if (fbar(q_tr) - gamma_p_n <= 1e-12 * (1.0 + std::fabs(gamma_p_n))) return out;

    const detail::HsDilatancy dil = detail::hs_dilatancy(p);
    auto sin_psi_m = [&](double q, double s3v) { return dil.from_q(q, s3v); };

    double dlam = 0.0, q = q_tr;
    Eigen::Vector3d sig = sig_tr;
    for (int it = 0; it < 50; ++it) {
        const double spm = sin_psi_m(q, sigma3_stiff);
        const double R = -(1.0 + spm) / (2.0 - spm);
        const double h = (4.0 + spm) / (2.0 - spm);  // 1−2R (γ^p hardening); 2 at ψ_m=0
        const Eigen::Vector3d mg(1.0, R, R);
        const Eigen::Vector3d Demg = De * mg;
        sig = sig_tr - dlam * Demg;
        q = sig(0) - sig(2);
        const double resid = fbar(q) - (gamma_p_n + h * dlam);
        const double dq_dl = -(Demg(0) - Demg(2));
        const double step = resid / (fbar_prime(q) * dq_dl - h);
        dlam -= step;
        if (std::fabs(step) <= 1e-14 * (1.0 + std::fabs(dlam))) break;
    }
    if (q > qf) { sig(0) -= (q - qf); q = qf; }  // MC failure plateau
    const double spm_f = sin_psi_m(q, sigma3_stiff);
    out.stress = sig;
    out.gamma_p = gamma_p_n + (4.0 + spm_f) / (2.0 - spm_f) * dlam;  // γ^p += h_s·dλ
    return out;
}

// --- Cap (volumetric) yield surface -----------------------------------------------
// f_c = q̃²/α² + p² − p_p²,  p = (σ1+σ2+σ3)/3,  q̃ = σ1+(δ−1)σ2−δσ3, δ=(3+sinφ)/(3−sinφ)
// (q̃=q at triaxial σ2=σ3). ASSOCIATED flow g_c=f_c ⇒ dε^p=λ ∂f_c/∂σ, dε_v^pc=λ·2p.
// POWER-LAW hardening ε_v^pc = (β/(1−m))(p_c/p_ref)^(1−m), modulus
// H_cap=(p_ref/β)(p_c/p_ref)^m. In isotropic compression (q̃=0) the tangent is
// K_iso = K_e·H_cap/(K_e+H_cap) (springs in series). The consistent tangent (associated) is
// ~symmetric.
// Schanz, Vermeer & Bonnier (1999); hardening-soil-formulation.md §4.
struct HsCapStep {
    Eigen::Vector3d stress;
    double pp;                // updated preconsolidation pressure
    Eigen::Matrix3d tangent;  // consistent (symmetric) D^ep
    bool plastic;
};

inline HsCapStep hs_cap_step(const HardeningSoilParams& p,
                             const Eigen::Vector3d& sigma_n, double pp_n,
                             const Eigen::Vector3d& deps) {
    const double s3 = sigma_n(2);
    const double Eur = p.Eur(s3), nu = p.nu_ur, alpha = p.cap_alpha;
    const Eigen::Matrix3d De = detail::hs_elastic(Eur, nu);
    const double ev_n = p.cap_ev_from_pc(pp_n);  // committed cap volumetric plastic strain

    auto mean = [](const Eigen::Vector3d& s) { return (s(0) + s(1) + s(2)) / 3.0; };
    // von Mises cap (symmetric): f_c = 3J2/α² + p² − pc². (δ·q̃ was asymmetric and broke the oedometer K0.)
    auto fcap = [&](const Eigen::Vector3d& s, double pp_var) {
        const double pm = mean(s);
        const double j3 = 0.5 * ((s(0) - s(1)) * (s(0) - s(1)) +
                                 (s(1) - s(2)) * (s(1) - s(2)) +
                                 (s(2) - s(0)) * (s(2) - s(0)));  // 3·J2
        return j3 / (alpha * alpha) + pm * pm - pp_var * pp_var;
    };
    // The cap gradient is LINEAR in σ: n_c = H_c·σ. H_c = (3/α²)(I−⅓·11ᵀ) + (2/9)·11ᵀ (von Mises).
    const Eigen::Matrix3d Hc = (3.0 / (alpha * alpha)) *
                               (Eigen::Matrix3d::Identity() - (1.0 / 3.0) * Eigen::Matrix3d::Ones()) +
                               (2.0 / 9.0) * Eigen::Matrix3d::Ones();
    const Eigen::Matrix3d DeHc = De * Hc;

    const Eigen::Vector3d sig_tr = sigma_n + De * deps;
    HsCapStep out;
    if (fcap(sig_tr, pp_n) <= 1e-10 * (1.0 + pp_n * pp_n)) {  // elastic
        out.stress = sig_tr; out.pp = pp_n; out.tangent = De; out.plastic = false;
        return out;
    }

    // Return mapping (backward Euler): σ(λ) = (I + λ·D_e·H_c)⁻¹·σ_tr (n_c linear ⇒ EXACT),
    // pp(λ) = pp_n + Kc·2·p(σ(λ))·λ. 1D Newton for f_c(σ(λ),pp(λ)) = 0 (an analytic
    // derivative is unnecessary — σ(λ) is exact, a numerical 1D derivative suffices and is
    // correct).
    double lam = 0.0;
    Eigen::Vector3d sig = sig_tr;
    double pp = pp_n;
    auto solve_sig = [&](double l) {
        return (Eigen::Matrix3d::Identity() + l * DeHc).inverse() * sig_tr;
    };
    auto resid = [&](double l, Eigen::Vector3d& s, double& ppv) {
        s = solve_sig(l);
        ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * l);  // power-law cap hardening
        return fcap(s, ppv);
    };
    for (int it = 0; it < 60; ++it) {
        const double r = resid(lam, sig, pp);
        const double dl = 1e-7 * (1.0 + std::fabs(lam));
        Eigen::Vector3d s2; double pp2;
        const double drdl = (resid(lam + dl, s2, pp2) - r) / dl;
        const double step = r / drdl;
        lam -= step;
        if (std::fabs(step) <= 1e-14 * (1.0 + std::fabs(lam))) { resid(lam, sig, pp); break; }
    }

    // Consistent ALGORITHMIC tangent. The cap surface is curved (the n_c direction varies
    // with σ) ⇒ the continuum tangent is not enough; the linearization of the return
    // mapping carries the curvature correction Ξ=(Ce+λ·H_c)⁻¹ (H_c = ∂²f_c/∂σ², a constant
    // Hessian). The hardening pp(σ,λ) (ε_v^pc=λ·2p) is coupled in:
    //   dλ = (w·Ξ d(dε))/(w·Ξ n_c + 4·pp·p·Kc),  w = n_c − (4·pp·Kc·λ/3)·1
    //   D_alg = Ξ − (Ξ n_c)(Ξ w)ᵀ/denom.
    const Eigen::Vector3d nc = Hc * sig;
    const double pmean = mean(sig);
    const double Hcap = p.cap_hardening_modulus(pp);  // stress-dependent (power law)
    Eigen::Matrix3d Ce;  // elastic compliance (Eur, ν)
    Ce << 1.0, -nu, -nu, -nu, 1.0, -nu, -nu, -nu, 1.0;
    Ce /= Eur;
    const Eigen::Matrix3d Xi = (Ce + lam * Hc).inverse();
    const Eigen::Vector3d w = nc - (4.0 * pp * Hcap * lam / 3.0) * Eigen::Vector3d::Ones();
    const Eigen::Vector3d Xinc = Xi * nc;
    const Eigen::Vector3d Xiw = Xi * w;
    const double denom = w.dot(Xinc) + 4.0 * pp * Hcap * pmean;
    out.tangent = Xi - (Xinc * Xiw.transpose()) / denom;
    out.stress = sig;
    out.pp = pp;
    out.plastic = true;
    return out;
}

// --- Two-surface (shear + cap) return mapping -------------------------------------
// Strain-driven, principal-space (fixed principal directions — for the single-element
// triaxial/oedometer drivers and calibration). Active set: shear-only / cap-only / both
// (Koiter). In the both case (oedometer) σ, λ_s, λ_c are solved coupled (inner fixed
// point + outer 2×2 Newton). The σ3 stiffness is frozen at the committed minor. pp
// hardens with the cap volumetric λ_c·2p, γ^p with the shear λ_s.
struct HsState {
    Eigen::Vector3d stress;
    double gamma_p;
    double pp;
};

// sigma3_stiff: the minor principal that freezes the elastic Eur (committed σ3). The FE
// wrapping passes the committed σ3; negative sentinel ⇒ sigma_n(2) is used
// (old call).
inline HsState hs_return_principal(const HardeningSoilParams& p,
                                   const Eigen::Vector3d& sigma_n, double gamma_p_n,
                                   double pp_n, const Eigen::Vector3d& dstrain,
                                   double sigma3_stiff = -1e300) {
    const double s3 = (sigma3_stiff > -1e299) ? sigma3_stiff : sigma_n(2);
    const double Eur = p.Eur(s3), nu = p.nu_ur, alpha = p.cap_alpha;
    const double ev_n = p.cap_beta > 0.0 ? p.cap_ev_from_pc(pp_n) : 0.0;
    const double Ei = p.Ei(s3), qa = p.q_asymptote(s3), qf = p.q_failure(s3);
    const Eigen::Matrix3d De = detail::hs_elastic(Eur, nu);
    const detail::HsDilatancy dil = detail::hs_dilatancy(p);
    // von Mises cap (symmetric) — δ·q̃ was asymmetric and broke the oedometer K0.
    const Eigen::Matrix3d Hc = (3.0 / (alpha * alpha)) *
                               (Eigen::Matrix3d::Identity() - (1.0 / 3.0) * Eigen::Matrix3d::Ones()) +
                               (2.0 / 9.0) * Eigen::Matrix3d::Ones();

    auto mean = [](const Eigen::Vector3d& s) { return (s(0) + s(1) + s(2)) / 3.0; };
    auto fbar = [&](double q) { return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur; };
    auto fbar_p = [&](double q) {
        const double r = 1.0 - q / qa; return (2.0 / Ei) / (r * r) - 2.0 / Eur;
    };
    auto fcap = [&](const Eigen::Vector3d& s, double pp) {
        const double pm = mean(s);
        const double j3 = 0.5 * ((s(0) - s(1)) * (s(0) - s(1)) +
                                 (s(1) - s(2)) * (s(1) - s(2)) +
                                 (s(2) - s(0)) * (s(2) - s(0)));  // 3·J2
        return j3 / (alpha * alpha) + pm * pm - pp * pp;
    };
    auto spm_of = [&](double q) { return dil.from_q(q, s3); };
    auto shear_dir = [&](double q) {
        const double spm = spm_of(q), R = -(1.0 + spm) / (2.0 - spm);
        return Eigen::Vector3d(1.0, R, R);
    };
    const bool cap_on = p.cap_beta > 0.0;

    const Eigen::Vector3d sig_tr = sigma_n + De * dstrain;
    HsState out{sig_tr, gamma_p_n, pp_n};
    const double fs_tr = fbar(sig_tr(0) - sig_tr(2)) - gamma_p_n;
    const double fc_tr = cap_on ? fcap(sig_tr, pp_n) : -1.0;
    if (fs_tr <= 1e-12 * (1.0 + std::fabs(gamma_p_n)) &&
        fc_tr <= 1e-10 * (1.0 + pp_n * pp_n))
        return out;  // elastic

    // Shear-only correction (assumes the cap inactive).
    auto shear_only = [&]() {
        double dl = 0.0, q = sig_tr(0) - sig_tr(2);
        Eigen::Vector3d sig = sig_tr;
        for (int it = 0; it < 50; ++it) {
            const double spm = spm_of(q);
            const double h = (4.0 + spm) / (2.0 - spm);  // 1−2R (γ^p hardening); 2 at ψ_m=0
            const Eigen::Vector3d Demg = De * shear_dir(q);
            sig = sig_tr - dl * Demg;
            q = sig(0) - sig(2);
            const double r = fbar(q) - (gamma_p_n + h * dl);
            const double step = r / (fbar_p(q) * (-(Demg(0) - Demg(2))) - h);
            dl -= step;
            if (std::fabs(step) <= 1e-14 * (1.0 + std::fabs(dl))) break;
        }
        if (q > qf) {
            // Failure plateau: perfectly plastic MC (q=qf), the flow (1,R,R) stays dilatant
            // (ψ_m=ψ@failure). σ=σ_tr−λ De m_g, q=qf ⇒ λ in closed form (R constant).
            // Dilation continues along the plateau; a raw σ1 clamp would kill it.
            const double spmf = spm_of(qf);
            const double Rf_ = -(1.0 + spmf) / (2.0 - spmf);
            const Eigen::Vector3d Demg = De * Eigen::Vector3d(1.0, Rf_, Rf_);
            const double q_tr = sig_tr(0) - sig_tr(2);
            const double lam = (q_tr - qf) / (Demg(0) - Demg(2));
            HsState s{sig_tr - lam * Demg, std::max(gamma_p_n, fbar(qf)), pp_n};
            return s;
        }
        const double spm_f = spm_of(q);
        HsState s{sig, gamma_p_n + (4.0 + spm_f) / (2.0 - spm_f) * dl, pp_n};
        return s;
    };
    // Cap-only correction (closed-form: n_c=H_c·σ is linear).
    auto cap_only = [&]() {
        const Eigen::Matrix3d DeHc = De * Hc;
        double lam = 0.0; Eigen::Vector3d sig = sig_tr; double pp = pp_n;
        auto solve = [&](double l, Eigen::Vector3d& s, double& ppv) {
            s = (Eigen::Matrix3d::Identity() + l * DeHc).inverse() * sig_tr;
            ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * l);  // power-law cap hardening
            return fcap(s, ppv);
        };
        for (int it = 0; it < 60; ++it) {
            const double r = solve(lam, sig, pp);
            const double dl = 1e-7 * (1.0 + std::fabs(lam));
            Eigen::Vector3d s2; double pp2;
            const double drdl = (solve(lam + dl, s2, pp2) - r) / dl;
            lam -= r / drdl;
            if (std::fabs(r) <= 1e-12 * (1.0 + pp_n * pp_n)) { solve(lam, sig, pp); break; }
        }
        HsState s{sig, gamma_p_n, pp}; return s;
    };
    // Both-active (Koiter): σ = σ_tr − λs De ns − λc De Hc σ; fs=0, fc=0. Inner fixed
    // point (σ) + outer 2×2 Newton (λs,λc), numerical Jacobian.
    auto both = [&]() {
        Eigen::Vector2d lam(0.0, 0.0);
        Eigen::Vector3d sig = sig_tr; double pp = pp_n, gp = gamma_p_n;
        auto residual = [&](const Eigen::Vector2d& l, Eigen::Vector3d& s,
                            double& ppv, double& gpv) {
            const double ls = std::max(l(0), 0.0), lc = std::max(l(1), 0.0);
            const Eigen::Matrix3d A = Eigen::Matrix3d::Identity() + lc * (De * Hc);
            const Eigen::Matrix3d Ainv = A.inverse();
            s = sig_tr;
            for (int k = 0; k < 30; ++k) {  // inner fixed point: ns(σ)
                const double q = s(0) - s(2);
                const Eigen::Vector3d sn = sig_tr - ls * (De * shear_dir(q));
                const Eigen::Vector3d s_new = Ainv * sn;
                if ((s_new - s).cwiseAbs().maxCoeff() <
                    1e-13 * (1.0 + s.cwiseAbs().maxCoeff())) { s = s_new; break; }
                s = s_new;
            }
            const double q = s(0) - s(2);
            const double spm = spm_of(q);
            gpv = gamma_p_n + (4.0 + spm) / (2.0 - spm) * ls;  // 1−2R (γ^p hardening); 2 at ψ_m=0
            ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * lc);  // power-law cap hardening
            Eigen::Vector2d r;
            r(0) = fbar(q) - gpv;
            r(1) = fcap(s, ppv);
            return r;
        };
        for (int it = 0; it < 50; ++it) {
            const Eigen::Vector2d r = residual(lam, sig, pp, gp);
            if (r.cwiseAbs().maxCoeff() <
                1e-11 * (1.0 + std::fabs(gamma_p_n) + pp_n * pp_n))
                break;
            Eigen::Matrix2d Jc;  // numerical 2×2 Jacobian
            for (int j = 0; j < 2; ++j) {
                Eigen::Vector2d lp = lam;
                const double dl = 1e-8 * (1.0 + std::fabs(lam(j)));
                lp(j) += dl;
                Eigen::Vector3d s2; double pp2, gp2;
                Jc.col(j) = (residual(lp, s2, pp2, gp2) - r) / dl;
            }
            lam -= Jc.inverse() * r;
        }
        residual(lam, sig, pp, gp);
        double q = sig(0) - sig(2);
        if (q > qf) { sig(0) -= (q - qf); }
        HsState s{sig, gp, pp}; return s;
    };

    // Active-set selection: try a single surface; if the other is violated, both.
    if (fs_tr > 0.0 && fc_tr <= 0.0) {
        HsState s = shear_only();
        if (!cap_on || fcap(s.stress, s.pp) <= 1e-8 * (1.0 + pp_n * pp_n)) return s;
        return both();
    }
    if (fc_tr > 0.0 && fs_tr <= 0.0) {
        HsState s = cap_only();
        if (fbar(s.stress(0) - s.stress(2)) - s.gamma_p <= 1e-10 * (1.0 + std::fabs(gamma_p_n)))
            return s;
        return both();
    }
    return both();  // both exceeded
}

// --- ROBUST multi-surface integrator: explicit substepping + Koiter (Sloan/Potts&Gens) ---
// The path the literature takes for complex soil models: the outer strain increment is split
// into small substeps; at each
// substep the active surfaces (shear/cap) are determined, the plastic multipliers are
// solved from the 2×2 system with Koiter multi-surface flow (negative multiplier → that
// surface is dropped), stress+hardening are updated, and a drift correction pulls back to
// the surface. Unlike an implicit nested Newton it does NOT DIVERGE (explicit, small step)
// — robust including high stiffness (E=100+MPa). Strain-driven (the shared path of
// integrate_point + single-element triaxial/oedometer + calibration).
// Source: Sloan, Abbo & Sheng (2001); Potts & Zdravković; Schanz, Vermeer & Bonnier (1999).
struct HsIntegrated {
    Eigen::Vector3d stress;
    double gamma_p;
    double pp;
    Eigen::Matrix3d tangent;  // continuum elastoplastic (last active set)
    bool plastic;
    int nsub;                 // substeps the error control asked for
    int saturated = 0;        // 1 when the guard ceiling clipped the subdivision, i.e. the
                              // integration did NOT meet its tolerance on this call
};

// The subdivision one integration used. A perturbed run (the numerical consistent tangent's
// three forward differences) is REPLAYED along exactly this subdivision: with automatic
// substepping the perturbation would otherwise change the subdivision itself, and that
// difference is larger than the ~E*h signal the difference quotient is trying to read. This is
// the automatic-substepping replacement for the fixed substep count the previous rule pinned.
// The subdivision is UNIFORM (one size, chosen once, with a partial last step), so recording it
// takes two numbers rather than a list -- and the replay can never overflow a buffer.
struct HsSubstepPlan {
    double dT = 0.0;   // pseudo-time size of every substep but the remainder; 0 = nothing recorded
    int n = 0;         // substeps taken
    void clear() { dT = 0.0; n = 0; }
};

// Integrate the Hardening Soil model over `dstrain` from the committed state, with the local
// integration error held under `stol`.
//
// Scheme: explicit substepping with AUTOMATIC ERROR CONTROL (Sloan 1987; Sloan, Abbo & Sheng
// 2001, Engineering Computations 18(1):121-194). Each substep is evaluated twice -- forward
// Euler from the start of the substep, and again from its Euler end point -- and the two are
// averaged (modified Euler, second order). Their difference IS the local truncation error, and
// it is what sizes the subdivision: the size is measured before the walk and then held (see
// "HOW MANY SUBSTEPS" below for why this scheme does not accept/reject as the paper does). A
// drift correction pulls each state back onto the surfaces it was yielding on.
//
// What this replaces, and why. Until 2026-08-20 the substep count was
// `ceil(||De.deps|| / (0.01 (max|sigma| + p_ref)))` -- a fixed fraction of a reference stress,
// with no error estimate, no tolerance and no rejection. The code cited the paper above
// without implementing its contribution. Measured on the corpus oedometer (KV-CST-002),
// refining that fraction from 1% to 0.03% moved the answer by TWO PERCENTAGE POINTS and did
// not converge -- an error axis larger than the model deviation the numerics record attributes
// to the cap calibration, and one no sweep in this project had ever touched. It also made the
// answer depend on the iteration path: a different line search, or a different linear-solver
// backend, lands on a different subdivision, which is how two backends came to disagree on a
// published number (docs/validation/numerical-uncertainty.md, and the 0.9.0 audit).
//
// `stol` is a NUMERICAL CONTROL, not a material property: it belongs to the phase, like the
// tolerated error and the load-step count. 0 means "the default below".
// TEMPORARY (0.9.0 N-1, second half): the tolerance's permanent home is the phase's numerical
// controls, next to the tolerated error and the load-step count, so that a published run carries
// the integration accuracy it was computed with. Until that plumbing lands it is a build-time
// default with an environment override, so the tolerance can be SWEPT -- which is the whole
// point of having one.
// The guard ceiling on one increment's subdivision, and the measurement that set it. On the
// corpus oedometer (KV-CST-002, range 100->200 kPa, STOL 1e-5) the ceiling costs and buys:
//     100 -> -0.9539%, 24 s | 200 -> -0.9861%, 32 s | 500 -> -1.0689%, 67 s | 2000 -> -1.0554%, 67 s
// Nearly all of the extra work goes into TRIAL iterates far from equilibrium, whose stresses are
// discarded; the answer moves by about a tenth of a percentage point across a twentyfold ceiling,
// which is this fixture's path indeterminacy, not integration error. 200 is the cheapest setting
// on that plateau. Saturation is reported (`HsIntegrated::saturated`) so a run that hit it
// is never mistaken for one that met its tolerance.
// TEMPORARY measurement seam (scaffold): the environment override, so the ceiling stays swept.
inline int hs_max_substeps() {
    static const int v = [] {
        const char* e = std::getenv("KATAI_HS_MAXSUB");
        const int d = e ? std::atoi(e) : 200;
        return d > 0 ? d : 200;
    }();
    return v;
}

// The seam that used to stand here (KATAI_HS_MCPROJ, returning the failure bound along the
// consistent direction De.n_s instead of clamping sigma1) is GONE, swept and removed on
// 2026-08-24. It was kept because a measurement said it fixed a 0.24% shortfall on the failure
// plateau at a cost of one boundary-value case. Re-measured, the fix is not there: both returns
// leave the deviator at the SAME 99.686% of q_f, and the cost still is. The shortfall was never
// the return direction -- see the bound below.

// The default integration tolerance, and the measurement that set it. KV-CST-002 (100->200 kPa,
// 40+160 increments, equilibrium tolerance 1e-6); the same material walked at the STRESS POINT
// says -1.003%, so that column is the target:
//     STOL    deviation   Newton iterations (seating/staged)   wall clock
//     1e-3    -1.536%     1887 / 4883                          91 s
//     1e-4    -1.310%     1907 / 1427                          76 s
//     1e-5    -0.986%      432 / 1186                          31 s
//     1e-6    -0.982%      387 / 3379                         117 s
// A LOOSER integration is not cheaper. Below about 1e-5 the integration noise sits above the
// residual the phase is asking for, and the equilibrium iteration grinds against it -- four times
// the seating iterations at 1e-3, for an answer half a percentage point further from the model's
// own. 1e-5 is the fastest setting AND the first one that reproduces the constitutive routine to
// the digits this case is published in.
inline double hs_default_substep_tol() {
    static const double v = [] {
        const char* e = std::getenv("KATAI_HS_STOL");
        const double d = e ? std::atof(e) : 1.0e-5;
        return d > 0.0 ? d : 1.0e-5;
    }();
    return v;
}

// The same value read as a whole-run OVERRIDE rather than as a default: 0 when the variable is
// not set. Since .k2d v15 a phase can name its own integration tolerance, and a study sweeping
// this axis has to be able to re-run a file that does -- otherwise the one thing the seam exists
// for, asking what a published number owes to its integration, is exactly what it cannot do on
// the files that state it. So the environment wins over the file here, which is the same
// precedence KATAI_CONV_NOLOCAL has over the phase's convergence setting.
inline double hs_env_substep_tol() {
    static const double v = [] {
        const char* e = std::getenv("KATAI_HS_STOL");
        const double d = e ? std::atof(e) : 0.0;
        return d > 0.0 ? d : 0.0;
    }();
    return v;
}

inline HsIntegrated hs_integrate(const HardeningSoilParams& p,
                                 const Eigen::Vector3d& sigma_n, double gamma_p_n,
                                 double pp_n, const Eigen::Vector3d& dstrain,
                                 double stol = 0.0,
                                 HsSubstepPlan* plan_out = nullptr,
                                 const HsSubstepPlan* plan_in = nullptr) {
    const double pr = p.p_ref, plim = 0.1 * pr;  // p_limit: σ3 floor in the stiffness laws
    const double nu = p.nu_ur;
    const detail::HsDilatancy dil = detail::hs_dilatancy(p);
    const double alpha = p.cap_alpha;

    // --- WHAT THE CURRENT STRESS DECIDES -----------------------------------------------------
    // E_ur, E_i, q_a and q_f are all functions of the minor principal stress, and in THIS model
    // that is not a detail -- stress-dependent stiffness is what makes it the Hardening Soil
    // model rather than Mohr-Coulomb with a cap. So they are read at the state each substep
    // starts from, which is what makes the substepping worth its cost.
    //
    // Until 2026-08-20 they were evaluated ONCE, at the state the whole increment started from,
    // and held fixed across every substep: the plastic flow was refined while the stiffness it
    // flowed against stayed at the beginning of the step. That error is FIRST ORDER in the OUTER
    // increment and no substep tolerance can see it, so an error-controlled integrator would have
    // reported "converged" while the answer still moved with the load-step count. Measured on the
    // oedometer walk (2% vertical strain, restrained lateral, tolerance 1e-6) the final sigma1 was
    //     20 steps 773.27 | 40 791.47 | 80 800.46 | 160 805.02 | 320 807.27 | 640 808.48 kPa
    // -- halving the outer step halved what was left, exactly first order, extrapolating to
    // ~809.7: the "converged" answer was 4.5% low at 20 steps and still 0.6% low at 160. Read per
    // substep the same sweep sits inside 6e-5 relative, on the number that sequence extrapolates
    // to (tests/study_hs_integration.cpp `point`, which is how both columns were measured).
    struct Stiff {
        double s3, Eur, Ei, qa, qf;
        Eigen::Matrix3d De;
    };
    auto stiff_at = [&](const Eigen::Vector3d& s) {
        Stiff k;
        k.s3 = std::max(s(2), plim);
        k.Eur = p.Eur(k.s3);
        k.Ei = p.Ei(k.s3);
        k.qa = p.q_asymptote(k.s3);
        k.qf = p.q_failure(k.s3);
        k.De = detail::hs_elastic(k.Eur, nu);
        return k;
    };
    const Stiff k_n = stiff_at(sigma_n);   // the committed state, for the elastic default
    // The cap deviatoric measure: symmetric von Mises q (q^2=3J2). This is EQUIVALENT to
    // the reduced form of the asymmetric q-tilde measure (s1+(delta-1)s2-delta s3) on the OEDOMETER
    // (axisymmetric, s2=s3) path: the s2,s3 components of the q-tilde flow average out under
    // axisymmetry -> exactly von Mises flow. (Applying q-tilde raw to a de2=de3=0 probe
    // produces the s2!=s3 absurdity; see hardening-soil-formulation.md section 4f.)
    // f_c=3J2/alpha^2+p^2-pc^2, with 3J2 and p^2 quadratic => n_c=H_c.sigma LINEAR
    // (closed-form return). H_c=(3/alpha^2)(I-(1/3)11^T)+(2/9)11^T.
    const Eigen::Matrix3d dev = Eigen::Matrix3d::Identity() -
                                (1.0 / 3.0) * Eigen::Matrix3d::Ones();
    const Eigen::Matrix3d Hc = (3.0 / (alpha * alpha)) * dev +
                               (2.0 / 9.0) * Eigen::Matrix3d::Ones();
    const bool cap_on = p.cap_beta > 0.0;
    const double ev_n = cap_on ? p.cap_ev_from_pc(pp_n) : 0.0;

    auto mean = [](const Eigen::Vector3d& s) { return (s(0) + s(1) + s(2)) / 3.0; };
    auto fbar = [&](double q, const Stiff& k) {
        return (2.0 / k.Ei) * q / (1.0 - q / k.qa) - 2.0 * q / k.Eur;
    };
    auto fbar_p = [&](double q, const Stiff& k) {
        const double r = 1.0 - q / k.qa; return (2.0 / k.Ei) / (r * r) - 2.0 / k.Eur;
    };
    auto fcap = [&](const Eigen::Vector3d& s, double ppv) {
        const double pm = mean(s);
        const double j3 = 0.5 * ((s(0) - s(1)) * (s(0) - s(1)) +
                                 (s(1) - s(2)) * (s(1) - s(2)) +
                                 (s(2) - s(0)) * (s(2) - s(0)));  // 3*J2 (von Mises q^2)
        return j3 / (alpha * alpha) + pm * pm - ppv * ppv;
    };
    auto spm_of = [&](double q, const Stiff& k) { return dil.from_q(q, k.s3); };
    auto pc_of = [&](double evv) { return cap_on ? p.cap_pc_from_ev(evv) : pp_n; };

    // --- ONE EXPLICIT INCREMENT from an arbitrary state: the modified-Euler building block ---
    // Evaluated twice per substep (at its start and at its Euler end point). It reads a state
    // and returns increments, so the pair differ only by where they were evaluated -- which is
    // exactly what makes their difference an error estimate.
    struct Inc {
        Eigen::Vector3d dsig = Eigen::Vector3d::Zero();
        double dgp = 0.0;
        double dev = 0.0;
        Eigen::Matrix3d tangent = Eigen::Matrix3d::Zero();
        bool plastic = false;
        bool as = false;       // shear surface active (after Koiter drops)
        bool ac = false;       // cap active
        bool at_fail = false;  // on the perfectly-plastic MC plateau
    };
    auto increment = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                         const Eigen::Vector3d& de_in) -> Inc {
        Inc r;
        const Stiff k = stiff_at(s_in);   // the stiffness and strength THIS substep starts from
        const Eigen::Matrix3d& De = k.De;
        const double qf = k.qf;
        r.dsig = De * de_in;
        r.tangent = De;
        const double pp = pc_of(ev_in);
        const Eigen::Vector3d sig_tr = s_in + De * de_in;
        const double q_tr = sig_tr(0) - sig_tr(2);
        bool as = fbar(q_tr, k) - gp_in > 1e-12 * (1.0 + std::fabs(gp_in));
        bool ac = cap_on && fcap(sig_tr, pp) > 1e-10 * (1.0 + pp * pp);
        if (!as && !ac) return r;
        r.plastic = true;

        // Gradients/flow/hardening of the active surfaces, at the state passed in.
        const double q = s_in(0) - s_in(2);
        // Flow direction m_g=(1,R,R), R=e3^p/e1^p. Rowe dilatancy: ev^p/eq^p=-sin(psi_m)
        // (DILATION, comp-pos => ev^p<0). R=-(1+sin psi_m)/(2-sin psi_m) gives it (-1/2 at
        // psi_m=0 = volumetrically neutral).
        const double spm = spm_of(q, k), R = -(1.0 + spm) / (2.0 - spm);
        // Hardening modulus h_s = 1-2R = (4+sin psi_m)/(2-sin psi_m); 2 at psi_m=0.
        const double h_s = (4.0 + spm) / (2.0 - spm);
        const Eigen::Vector3d n_s(1.0, R, R);  // flow direction (dilatant)
        // Failure plateau: at q>=qf the shear becomes PERFECTLY-PLASTIC MC (yield f=q-qf,
        // grad (1,0,-1), hardening 0); the flow (1,R,R) stays dilatant -> dilation continues
        // along the plateau.
        const bool at_fail = q >= qf - 1e-9 * (1.0 + qf);
        const Eigen::Vector3d m_s = at_fail ? Eigen::Vector3d(1.0, 0.0, -1.0)
                                            : (fbar_p(q, k) * Eigen::Vector3d(1.0, 0.0, -1.0));
        const double Hh_s = at_fail ? 0.0 : h_s;
        const Eigen::Vector3d n_c = Hc * s_in, m_c = n_c;  // cap associated
        const double pmean = mean(s_in);
        const double Hcap = cap_on ? p.cap_hardening_modulus(pp) : 0.0;
        r.at_fail = at_fail;

        // Active-set Koiter solve (negative multiplier -> drop the surface, solve again).
        for (int pass = 0; pass < 3; ++pass) {
            const int na = (as ? 1 : 0) + (ac ? 1 : 0);
            if (na == 0) { r.dsig = De * de_in; r.tangent = De; break; }
            Eigen::MatrixXd A(na, na); Eigen::VectorXd b(na);
            std::vector<Eigen::Vector3d> ns, ms;
            std::vector<double> Hh;
            if (as) { ns.push_back(n_s); ms.push_back(m_s); Hh.push_back(Hh_s); }
            if (ac) { ns.push_back(n_c); ms.push_back(m_c);
                      Hh.push_back(4.0 * pp * pmean * Hcap); }
            for (int i = 0; i < na; ++i) {
                b(i) = ms[i].dot(De * de_in);
                for (int j = 0; j < na; ++j)
                    A(i, j) = ms[i].dot(De * ns[j]) + (i == j ? Hh[i] : 0.0);
            }
            const Eigen::VectorXd dl = A.fullPivLu().solve(b);
            bool dropped = false;
            int idx = 0;
            double dl_s = 0.0, dl_c = 0.0;
            if (as) { dl_s = dl(idx++); if (dl_s < 0.0) { as = false; dropped = true; } }
            if (ac) { dl_c = dl(idx++); if (dl_c < 0.0) { ac = false; dropped = true; } }
            if (dropped) continue;
            Eigen::Vector3d plastic_strain = Eigen::Vector3d::Zero();
            if (as) plastic_strain += dl_s * n_s;
            if (ac) plastic_strain += dl_c * n_c;
            r.dsig = De * (de_in - plastic_strain);
            if (as && !at_fail) r.dgp = h_s * dl_s;  // plateau: hardening freezes
            if (ac) r.dev = dl_c * 2.0 * pmean;
            {   // continuum tangent for this active set
                const int m = na;
                Eigen::MatrixXd N(3, m), M(3, m);
                int c = 0;
                if (as) { N.col(c) = n_s; M.col(c) = m_s; ++c; }
                if (ac) { N.col(c) = n_c; M.col(c) = m_c; ++c; }
                r.tangent = De - De * N * A.fullPivLu().solve(M.transpose() * De);
            }
            break;
        }
        r.as = as;
        r.ac = ac;
        return r;
    };

    // --- DRIFT CORRECTION, applied to an ACCEPTED substep ------------------------------------
    // SHEAR: a stress violating the surface is pulled along the CONSISTENT (Potts & Gens 1985)
    // elastoplastic direction De.n_s -- dsig=-(f/(a^T De.n_s + h_s)) De.n_s, a=df/dsigma.
    // De.n_s is SYMMETRIC in s2,s3 (n_s=(1,R,R)) => s2=s3 IS PRESERVED (the axisymmetric /
    // oedometer edge). The bare yield gradient (1,0,-1) pushes only s1,s3 and leaves s2 -> a
    // s_r != s_theta drift in triaxial. The CAP is associated (n_c=Hc.sigma symmetric) => the
    // gradient projection never disturbs s2,s3 anyway. The surfaces corrected are the ones the
    // substep was yielding on (the union of its two evaluations).
    auto correct_drift = [&](Eigen::Vector3d& s, double& gp_io, double ev_io,
                             bool as_act, bool ac_act, bool at_fail) {
        for (int it = 0; it < 5; ++it) {
            const Stiff k = stiff_at(s);   // the correction moves sigma3, so it moves these too
            const double qf = k.qf;
            const double ppc = pc_of(ev_io);
            const double qd = s(0) - s(2);
            const double fs = at_fail ? (qd - qf) : (fbar(qd, k) - gp_io);
            const double fcp = cap_on ? fcap(s, ppc) : -1.0;
            bool corr = false;
            if (as_act && fs > 1e-9 * (1.0 + std::fabs(gp_io) + qf)) {
                const double spm = spm_of(qd, k), R = -(1.0 + spm) / (2.0 - spm);
                const double h_s = (4.0 + spm) / (2.0 - spm);
                const Eigen::Vector3d n_s(1.0, R, R);
                const Eigen::Vector3d a = (at_fail ? 1.0 : fbar_p(qd, k)) *
                                          Eigen::Vector3d(1.0, 0.0, -1.0);  // df/dsigma
                const Eigen::Vector3d Den = k.De * n_s;                     // flow direction
                const double denom = a.dot(Den) + (at_fail ? 0.0 : h_s);
                const double dlam_d = fs / denom;
                s -= dlam_d * Den;
                if (!at_fail) gp_io += h_s * dlam_d;  // hardening tracks the drift step
                corr = true;
            }
            if (ac_act && std::fabs(fcp) > 1e-8 * (1.0 + ppc * ppc)) {
                const Eigen::Vector3d g = Hc * s;
                s -= (fcp / g.squaredNorm()) * g; corr = true;
            }
            if (!corr) break;
        }
        // MC failure bound (shear): q <= q_f. This is the last line -- the drift loop above only
        // pulls back onto the surface the substep was YIELDING on, and a substep that crossed into
        // failure can leave the state outside the failure surface with the hardening surface
        // satisfied.
        //
        // It is a bare clamp on sigma1, which leaves sigma3 alone, and that is ONE-SIDED: it only
        // ever lowers sigma1, and nothing gives back what it took when q_f rises again. So it
        // ratchets. WHAT it ratchets on is the point, and the record had that wrong until
        // 2026-08-24: it ratchets on the substep INTEGRATION ERROR in sigma3, not on the shape of
        // the return, so the size of the shortfall is set by STOL. Berlin Sand III drained
        // triaxial (study_hs_integration triax, sigma3 = 200, q_f = 644.85), deviator as a
        // fraction of q_f at 5% and at 20% axial strain:
        //
        //     STOL    at 5%       at 20%
        //     1e-3    98.671%      96.177%   <- FALLS with strain: the artefact imitates softening
        //     1e-5    99.686%      99.686%   (the shipped default)
        //     1e-7    99.949%      99.984%
        //     1e-9    99.993%     100.011%
        //
        // The plateau IS reached. What stands between the shipped run and it is integration error,
        // a declared and controllable quantity (KATAI_HS_STOL, whose default was set by its own
        // measurement above), not a defect in the bound. The consistent De.n_s projection that
        // used to sit behind KATAI_HS_MCPROJ was recorded as the cure and is not one: at the
        // shipped STOL it leaves the deviator at the SAME 99.686%, differing only in the fifth
        // decimal of volumetric strain, while still costing the shear-dominated low-confinement
        // strip footing (test_hs_footing stops converging at the full service load, 1.000 ->
        // 0.815, 1005 iterations). A branch with no measured benefit and a measured cost is not a
        // seam, it is a second program, so it is gone. The ratchet is also CONSERVATIVE -- it
        // under-predicts strength -- which is what makes the residual affordable at the shipped
        // tolerance rather than merely tolerable.
        {
            const Stiff k = stiff_at(s);
            const double over = (s(0) - s(2)) - k.qf;
            if (over > 1e-12 * (1.0 + k.qf)) s(0) -= over;
        }
    };

    Eigen::Vector3d sig = sigma_n;
    double gp = gamma_p_n, ev = ev_n;
    bool any_plastic = false;
    // The committed state's elastic matrix, until a plastic substep replaces it with its own
    // continuum tangent (an elastic increment returns this one).
    Eigen::Matrix3d tangent = k_n.De;

    const double env_tol = hs_env_substep_tol();
    const double tol = env_tol > 0.0 ? env_tol
                                     : (stol > 0.0 ? stol : hs_default_substep_tol());
    // Ceiling on the measured subdivision. It is a GUARD, not a policy: an increment that asks
    // for more than this is an increment the load stepping should have cut, and the counter below
    // records every time it fires so a saturated integration is never silently reported as one
    // that met its tolerance.
    const int kMaxSubsteps = hs_max_substeps();
    if (plan_out) plan_out->clear();

    const bool replaying = plan_in && plan_in->n > 0 && plan_in->dT > 0.0;

    // --- HOW MANY SUBSTEPS: measured from the error, and CONTINUOUS in the strain increment ---
    //
    // Two decisions, both forced by measurement rather than taste.
    //
    // (1) The subdivision is chosen ONCE for the increment, from a real error estimate, and then
    //     used as-is. The textbook accept/reject form of Sloan, Abbo & Sheng (2001) re-subdivides
    //     as it goes, which makes the subdivision a discontinuous function of the strain
    //     increment: two neighbouring Newton iterates integrate along different subdivisions, the
    //     internal force stops being a smooth function of the displacement, and the global
    //     iteration cannot converge below the integration noise. MEASURED on the corpus
    //     oedometer, 2026-08-20: with accept/reject the residual fell to a relative 1.2e-3 and
    //     then wandered, the line search halving to 2.4e-4 without finding descent.
    //
    // (2) The count is a REAL number, not an integer. With ceil() the map still jumps -- the same
    //     measurement, one step further: the residual reached 5.7e-6 and stalled there, because
    //     an iterate that changes n from 4 to 5 changes the answer by the difference between two
    //     subdivisions. Taking n substeps of 1/n_real with a partial last one makes the answer a
    //     CONTINUOUS function of the strain increment, which is what a Newton iteration needs to
    //     converge to the tolerance a phase asks for.
    //
    // The estimate itself: a modified-Euler pair over a trial substep measures the local error of
    // that substep. The scheme is second order, so the local error of a substep of pseudo-time dT
    // is ~C.dT^3 and the error accumulated over the 1/dT of them is ~C.dT^2; the largest dT
    // meeting the tolerance is sqrt(STOL.dT^3/e). Read at dT=1 that is the familiar
    // sqrt(err_full/STOL) -- but a whole increment is often far outside the asymptotic regime
    // that law assumes, so the estimate is taken TWICE: once over the full increment, then again
    // over the substep the first pass proposed, where the law does hold. Both passes are
    // continuous functions of the strain increment; no branch chooses between them.
    double nreal = 1.0;
    if (replaying) {
        // Replay: the subdivision is uniform, so its one step size reproduces it exactly.
        nreal = 1.0 / plan_in->dT;
    } else if (dstrain.squaredNorm() > 0.0) {
        for (int pass = 0; pass < 2; ++pass) {
            const double dT = 1.0 / nreal;
            const Eigen::Vector3d de = dT * dstrain;
            const Inc f1 = increment(sigma_n, gamma_p_n, ev_n, de);
            const Inc f2 = increment(sigma_n + f1.dsig, gamma_p_n + f1.dgp, ev_n + f1.dev, de);
            // No plasticity gate on the estimate. An elastic step is not error-free here: E_ur
            // moves with sigma3, so the elastic response is nonlinear and its own error is what
            // the pair measures. Where the stiffness really is constant the two evaluations
            // agree, the estimate is zero, and the step is taken whole -- a gate would only hide
            // the case where it is not.
            const Eigen::Vector3d s_end = sigma_n + 0.5 * (f1.dsig + f2.dsig);
            // The denominator needs a FLOOR. A purely relative measure against ||sigma|| is
            // degenerate where a run starts from (near) zero stress -- a weightless column, a
            // surface layer, the first increment of any seating phase. The floor is p_ref/200:
            // tied to the reference pressure of the stress-dependent stiffness law, so it scales
            // with the material's own stress level rather than with the units.
            const double dn = std::max(s_end.norm(), pr / 200.0);
            const double e_step = 0.5 * (f2.dsig - f1.dsig).norm() / dn;
            if (!(e_step > 0.0)) break;             // exact on this step: nothing to subdivide
            const double dT_ok = std::sqrt(tol * dT * dT * dT / e_step);
            nreal = std::min(std::max(1.0, 1.0 / dT_ok), (double)kMaxSubsteps);
            if (nreal <= 1.0) break;
        }
    }

    const double dT_full = 1.0 / nreal;
    int taken = 0;
    if (plan_out) { plan_out->clear(); plan_out->dT = dT_full; }
    for (double T = 0.0; T < 1.0 - 1e-12 && taken < kMaxSubsteps + 2; ++taken) {
        const double dT = std::min(dT_full, 1.0 - T);
        const Eigen::Vector3d de = dT * dstrain;
        const Inc k1 = increment(sig, gp, ev, de);
        const Inc k2 = increment(sig + k1.dsig, gp + k1.dgp, ev + k1.dev, de);
        sig += 0.5 * (k1.dsig + k2.dsig);
        gp += 0.5 * (k1.dgp + k2.dgp);
        ev += 0.5 * (k1.dev + k2.dev);
        if (k1.plastic || k2.plastic) {
            any_plastic = true;
            tangent = k2.plastic ? k2.tangent : k1.tangent;
            correct_drift(sig, gp, ev, k1.as || k2.as, k1.ac || k2.ac,
                          k1.at_fail || k2.at_fail);
        }
        if (plan_out) ++plan_out->n;
        T += dT;
    }
    // Saturation is REPORTED, never absorbed: if the ceiling clipped the subdivision the
    // integration did not meet its tolerance, and the caller is told so rather than handed a
    // number that looks like every other one.
    const int accepted = taken, clipped = (nreal >= (double)kMaxSubsteps) ? 1 : 0;

    HsIntegrated out;
    out.stress = sig;
    out.gamma_p = gp;
    out.pp = cap_on ? p.cap_pc_from_ev(ev) : pp_n;
    out.tangent = tangent;
    out.plastic = any_plastic;
    out.nsub = accepted;
    out.saturated = clipped;
    return out;
}

// Cap preconsolidation initialization (FE initial state): pp = p_eq · OCR,
// p_eq = √(3J2/α² + p²) (the isotropic-equivalent pressure at which the cap passes through
// σ0, the state parameter p_eq). NC (OCR=1) ⇒ the initial state sits on the
// cap (f_c=0); this keeps the cap well-defined in FE (pp_n=0 would make the cap yield
// always). sig_comp_pos = compression-positive principal stress (σ_HS = −σ_solver).
inline double hs_initial_pp(const HardeningSoilParams& p,
                            const Eigen::Vector3d& sig_comp_pos, double OCR = 1.0) {
    const double pm = (sig_comp_pos(0) + sig_comp_pos(1) + sig_comp_pos(2)) / 3.0;
    const double j3 = 0.5 * ((sig_comp_pos(0) - sig_comp_pos(1)) * (sig_comp_pos(0) - sig_comp_pos(1)) +
                             (sig_comp_pos(1) - sig_comp_pos(2)) * (sig_comp_pos(1) - sig_comp_pos(2)) +
                             (sig_comp_pos(2) - sig_comp_pos(0)) * (sig_comp_pos(2) - sig_comp_pos(0)));
    const double a = p.cap_alpha;
    return std::sqrt(j3 / (a * a) + pm * pm) * OCR;
}

// Shear hardening initialization (FE initial state): γ^p = f̄(q0). The geostatic K0 state
// carries q0>0; starting with γ^p=0 the shear surface (f_s=f̄(q0)−γ^p) is violated FROM THE
// START → the first return makes a large inconsistent correction (the solver diverges).
// γ^p=f̄(q0) seats the state on the surface (admissible pre-stress), the shear counterpart
// of the pp init (cap). sig_comp_pos = compression-positive principals.
inline double hs_initial_gamma_p(const HardeningSoilParams& p,
                                 const Eigen::Vector3d& sig_comp_pos) {
    const double s1 = sig_comp_pos.maxCoeff(), s3 = sig_comp_pos.minCoeff();
    const double Ei = p.Ei(s3), qa = p.q_asymptote(s3), Eur = p.Eur(s3), qf = p.q_failure(s3);
    const double q = std::min(std::max(s1 - s3, 0.0), qf);
    const double fbar = (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur;
    return std::max(fbar, 0.0);
}

// --- Oedometer probe + cap calibration (on the robust hs_integrate) ------------------
// Oedometer (ε_h=0) NC primary loading (with hs_integrate, robust); returns the tangent
// Eoed at σ1=p_ref and the lateral ratio K0. (α,β) → (K0,Eoed).
// The tolerance the CALIBRATION integrates at. Deliberately independent of, and tighter than,
// whatever a run uses: alpha and beta are properties of the MATERIAL (they are what makes the
// model reproduce Eoed_ref and K0_NC), so they must not move when a phase asks for a looser or
// tighter integration. Before 2026-08-20 the probe inherited the integrator's fixed 1% substep
// rule, which tied every calibrated cap -- and therefore every Hardening Soil answer this
// program has ever produced -- to a numerical constant, in a second and completely invisible
// way. This is the fix for that half of it.
inline constexpr double kHsCalibrationTol = 1.0e-6;

inline void hs_oedometer_probe(const HardeningSoilParams& p, double& Eoed_pref,
                               double& K0) {
    const double pr = p.p_ref, p0 = 0.02 * pr;
    Eigen::Vector3d sig(p0, p0, p0);
    double gp = 0.0, pp = p0, s1b = p0;
    // Two step sizes, and the reason is the whole point of an error-controlled integrator: the
    // OUTER step no longer decides accuracy, the tolerance does, so the walk up to the reference
    // pressure can be taken in few, large steps (each subdivided internally as much as it needs)
    // and only the reading itself needs a small step -- Eoed is read as a SECANT over one step,
    // so that one must stay short to approximate the tangent. Measured: the coarse approach
    // costs about a tenth of the uniform fine walk and moves neither Eoed nor K0 in the figures
    // the calibration bisects on.
    const double de_coarse = 2.0e-3, de_fine = 2.0e-4;
    Eoed_pref = 0.0; K0 = 0.0;
    for (int i = 0; i < 3000; ++i) {
        const double de1 = (sig(0) < 0.9 * pr) ? de_coarse : de_fine;
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, 0, 0),
                                            kHsCalibrationTol);
        if (s1b < pr && r.stress(0) >= pr) {
            Eoed_pref = (r.stress(0) - s1b) / de1;
            K0 = r.stress(2) / r.stress(0);
        }
        s1b = r.stress(0); sig = r.stress; gp = r.gamma_p; pp = r.pp;
        if (sig(0) > 1.25 * pr) break;
    }
}

// Cap hardening modulus K_p = K1·K2/(K1−K2) — CLOSED FORM from Eoed_ref (elastic and
// plastic bulk springs in series, 1/K2 = 1/K1 + 1/K_p). K1=Eur_ref/(3(1−2ν)) unloading bulk;
// K2=Eoed_ref(1+2K0nc)/3. This ties the cap hardening directly to Eoed (no numerical
// calibration of β needed) ⇒ the K0-Eoed coupling is resolved: β=p_ref/(k·K_p), only α
// (and a small k correction) is calibrated.
inline double hs_cap_Kp(const HardeningSoilParams& p, double K0_NC) {
    const double K1 = p.Eur_ref / (3.0 * (1.0 - 2.0 * p.nu_ur));
    const double K2 = p.Eoed_ref * (1.0 + 2.0 * K0_NC) / 3.0;
    return (K1 > K2) ? K1 * K2 / (K1 - K2) : K1;  // K1>K2 (Eur≫Eoed) typical
}

// Cap parameters (α, β) from the standard HS inputs (Eoed_ref, K0_NC). β = p_ref/(k·K_p) (K_p
// closed form, sets Eoed in CLOSED FORM → the coupling is resolved); α by outer bisection for
// K0_NC; k by inner bisection as an Eoed_ref fine correction (k≈1). Both are fitted by
// simulated oedometer tests. If K0_NC is unreachable, the nearest α.
inline void hs_calibrate_cap(HardeningSoilParams& p, double K0_NC) {
    const double Kp = hs_cap_Kp(p, K0_NC);
    // Inner: at a given α find k (β=p_ref/(k·K_p)) for Eoed_ref; return K0.
    auto solve_k_return_K0 = [&](double alpha) {
        p.cap_alpha = alpha;
        double klo = 0.2, khi = 5.0;  // k correction factor (around ≈1)
        for (int it = 0; it < 40; ++it) {
            const double km = 0.5 * (klo + khi);
            p.cap_beta = p.p_ref / (km * Kp);
            double e, k0; hs_oedometer_probe(p, e, k0);
            if (e > p.Eoed_ref) khi = km; else klo = km;  // Eoed grows with k (stiff cap)
            if (khi - klo < 1e-4) break;
        }
        p.cap_beta = p.p_ref / (0.5 * (klo + khi) * Kp);
        double e, k0; hs_oedometer_probe(p, e, k0);
        return k0;
    };
    // Outer: bisect α for K0=K0_NC (K0 falls with α).
    double alo = 0.1, ahi = 80.0;
    const double klo = solve_k_return_K0(alo), khi = solve_k_return_K0(ahi);
    if ((K0_NC - klo) * (K0_NC - khi) <= 0.0) {  // bracketed → reachable
        for (int it = 0; it < 30; ++it) {
            const double am = 0.5 * (alo + ahi), km = solve_k_return_K0(am);
            if (km > K0_NC) alo = am; else ahi = am;
            if (ahi - alo < 2e-3) break;
        }
        solve_k_return_K0(0.5 * (alo + ahi));
    } else {
        solve_k_return_K0(std::fabs(klo - K0_NC) < std::fabs(khi - K0_NC) ? alo : ahi);
    }
}

}  // namespace katai::core
