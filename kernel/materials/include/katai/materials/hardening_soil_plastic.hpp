#pragma once
// Hardening Soil — P2.3d (part 1): the multiaxial SHEAR hardening return mapping.
// Strain-driven predictor-corrector (FE constitutive form) in principal stress space.
// Elastic: stress-dependent Eur (frozen at σ3 at the start of the step). Shear yield
// surface (Schanz 1999):
//   f = f̄(q) − γ^p ,   f̄(q) = (2/Ei) q/(1−q/qa) − 2q/Eur ,   q = σ1 − σ3
// Non-associated flow: mobilized dilatancy ψ_m, m_g = (1, R, R), R = ε3^p/ε1^p. ψ_m is Rowe's
// where Rowe is non-negative; below the phase-transformation line HS takes zero and HSsmall
// takes Li & Dafalias instead — one rule, `detail::hs_dilatancy` (see its comment).
//   Rowe (PLAXIS Eq 6-14): ε_v^p/ε_q^p = −sinψ_m (DILATION, comp-pos ⇒ ε_v^p<0) ⇒
//   R = −(1+sinψ_m)/(2−sinψ_m) ; hardening (Eq 6-10) γ^p=−(2ε1^p−ε_v^p) ⇒ dγ^p=h·dλ,
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
// HS (Ref. Man. Eq 6-14, Rowe):  sinψ_m = (sinφ_m − sinφ_cv)/(1 − sinφ_m·sinφ_cv), cut off to
// [0, sinψ]. Where Rowe returns a NEGATIVE value the Hardening Soil model takes ψ_m = 0.
//
// HSsmall (Mat. Models Man. sec 7.9.1, after Li & Dafalias (2000)): that zero cut-off "may
// sometimes yield too little plastic volumetric strains", so wherever Rowe is negative the
// small-strain model puts a small CONTRACTION there instead of nothing:
//   sinψ_m = (1/10)·(−M_c·exp[(1/15)·ln((M_c/M_d)·(q/q_a))] + M_d)          (Eq 7-19)
//   M_c    = 6·sinφ_cv/(3 − sinφ_cv)                                        (Eq 7-20)
//   M_d    = 6·sinφ_m /(3 − sinφ_m )                                        (Eq 7-21)
//   q/q_a  = max([(1−sinφ_cv)/sinφ_cv]·[sinφ_m/(1−sinφ_m)], 1e−4)           (Eq 7-22)
//   sinφ_m ≥ sinφ/(2 − sinφ)                                                (Eq 7-23)
// Eq 7-22's q/q_a is its own definition, NOT the model's deviatoric hyperbola ratio q/q_a of
// Eq 6-3; the two are unrelated quantities that share a name.
//
// The two branches MEET EXACTLY, and that is the property worth knowing: at φ_m = φ_cv we get
// M_d = M_c and q/q_a = 1, so the logarithm vanishes and Eq 7-19 returns exactly 0 — which is
// also precisely where Rowe changes sign. The composite rule is continuous by CONSTRUCTION,
// not to a tolerance. Eq 7-23 then bounds how contractant it can get, since sinψ_m falls
// monotonically as φ_m drops: the floor is the most negative ψ_m the model can produce.
// (test_hssmall pins both.)
//
// MEASURED AGAINST THE MANUAL'S OWN FIGURE 7-10 — and they disagree. Digitising that plot
// (φ=35°, ψ=5°) confirms every structural feature of Eq 7-19: the upper branch is pure Rowe to
// ±0.01°, the curve vanishes at φ_cv = 30.80°, the plateau starts at the Eq 7-23 floor of
// 23.71°, and the 1/15 exponent reproduces the plot's curvature to better than 1%. Only the
// AMPLITUDE differs, by a constant factor: the figure is 1.29× the printed equation
// everywhere (plateau −2.167° drawn against −1.679° computed), i.e. the plot behaves as if the
// leading 1/10 were 1/7.74. We implement the EQUATION, because an equation printed in the
// specification outranks a constant reverse-engineered from a raster plot, and because 1.29 is
// not a number the manual states anywhere. The gap is declared in hssmall-formulation.md; both
// readings are small contractions and the difference between them is far smaller than the
// difference between either and the ψ_m = 0 this replaces.
struct HsDilatancy {
    double sin_cs = 0.0;        // sinφ_cv (Rowe, from the input φ and ψ)
    double sin_psi = 0.0;       // sinψ — the upper cut-off
    double c_cot = 0.0;         // c·cotφ — the cohesion shift of the mobilized-φ definition
    double Mc = 0.0;            // Eq 7-20 (HSsmall branch only)
    double sphi_m_floor = 0.0;  // Eq 7-23 (HSsmall branch only)
    bool li_dafalias = false;   // HSsmall (G0_ref>0) with a usable φ_cv

    double operator()(double sphi_m) const {
        const double rowe = (sphi_m - sin_cs) / (1.0 - sphi_m * sin_cs);
        if (rowe >= 0.0) return std::min(rowe, sin_psi);
        if (!li_dafalias) return 0.0;                                   // HS: the zero cut-off
        const double s = std::max(sphi_m, sphi_m_floor);                // Eq 7-23
        const double Md = 6.0 * s / (3.0 - s);                          // Eq 7-21
        const double qqa =
            std::max((1.0 - sin_cs) / sin_cs * (s / (1.0 - s)), 1e-4);  // Eq 7-22
        return 0.1 * (-Mc * std::exp(std::log((Mc / Md) * qqa) / 15.0) + Md);  // Eq 7-19
    }

    // sinφ_m from the mobilized stress state (Eq 7-18 / Eq 6-13, σ1 = q + σ3, comp-positive).
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
    // φ_cv = 0 (a φ = 0 Tresca soil) would divide by zero in Eq 7-22 — but it also makes Rowe
    // non-negative everywhere, so the branch is unreachable there; the guard says so rather
    // than relying on it.
    d.li_dafalias = p.G0_ref > 0.0 && !p.dilatancy_cut && d.sin_cs > 1e-12 && sphi < 1.0;
    if (d.li_dafalias) {
        d.Mc = 6.0 * d.sin_cs / (3.0 - d.sin_cs);  // Eq 7-20
        d.sphi_m_floor = sphi / (2.0 - sphi);      // Eq 7-23
    }
    return d;
}

}  // namespace detail

// One shear-hardening step: committed (σ_n, γ_n) + strain increment dε → updated state +
// consistent tangent. Eur/Ei/qa are frozen at the start-of-step σ3 = σ_n[2] (explicit
// stress-dependent stiffness, like PLAXIS). q is capped at qf at failure (perfectly
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
        h = (4.0 + spm) / (2.0 - spm);  // 1−2R (Eq 6-10: γ^p=−(2ε1^p−ε_v^p)); 2 at ψ_m=0
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
    h = (4.0 + spm) / (2.0 - spm);  // 1−2R (Eq 6-10); 2 at ψ_m=0
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
        const double h = (4.0 + spm) / (2.0 - spm);  // 1−2R (Eq 6-10); 2 at ψ_m=0
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
    out.gamma_p = gamma_p_n + (4.0 + spm_f) / (2.0 - spm_f) * dlam;  // γ^p += h_s·dλ (Eq 6-10)
    return out;
}

// --- Cap (volumetric) yield surface -----------------------------------------------
// f_c = q̃²/α² + p² − p_p²,  p = (σ1+σ2+σ3)/3,  q̃ = σ1+(δ−1)σ2−δσ3, δ=(3+sinφ)/(3−sinφ)
// (q̃=q at triaxial σ2=σ3). ASSOCIATED flow g_c=f_c ⇒ dε^p=λ ∂f_c/∂σ, dε_v^pc=λ·2p.
// POWER-LAW hardening (Eq 15.13): p_c↔ε_v^pc, modulus H_cap=(p_ref/β)(p_c/p_ref)^m. In
// isotropic compression (q̃=0) the tangent is K_iso = K_e·H_cap/(K_e+H_cap) (springs in
// series). The consistent tangent (associated) is ~symmetric.
// Schanz (1999) / Rocscience-PLAXIS HS; hardening-soil-formulation.md §4.
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
        ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * l);  // power-law hardening (Eq 15.13)
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
    const double Hcap = p.cap_hardening_modulus(pp);  // stress-dependent (Eq 15.13)
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
// wrapping passes the committed σ3 (like PLAXIS); negative sentinel ⇒ sigma_n(2) is used
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
            const double h = (4.0 + spm) / (2.0 - spm);  // 1−2R (Eq 6-10); 2 at ψ_m=0
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
            // Dilation continues along the plateau (Fig 15.4); a raw σ1 clamp would kill it.
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
            ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * l);  // power law (Eq 15.13)
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
            gpv = gamma_p_n + (4.0 + spm) / (2.0 - spm) * ls;  // 1−2R (Eq 6-10); 2 at ψ_m=0
            ppv = p.cap_pc_from_ev(ev_n + 2.0 * mean(s) * lc);  // power law (Eq 15.13)
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
// The path leading geotechnical codes (the PLAXIS/GEO5/Midas class) and the literature take
// for complex soil models: the outer strain increment is split into small substeps; at each
// substep the active surfaces (shear/cap) are determined, the plastic multipliers are
// solved from the 2×2 system with Koiter multi-surface flow (negative multiplier → that
// surface is dropped), stress+hardening are updated, and a drift correction pulls back to
// the surface. Unlike an implicit nested Newton it does NOT DIVERGE (explicit, small step)
// — robust including high stiffness (E=100+MPa). Strain-driven (the shared path of
// integrate_point + single-element triaxial/oedometer + calibration).
// Source: Sloan, Abbo & Sheng (2001); Potts & Zdravković; Rocscience/PLAXIS HS (Eq 15.1-15.13).
struct HsIntegrated {
    Eigen::Vector3d stress;
    double gamma_p;
    double pp;
    Eigen::Matrix3d tangent;  // continuum elastoplastic (last active set)
    bool plastic;
    int nsub;                 // substep count used (FD tangent perturbations run with the
                              // same nsub — see nsub_fixed)
};

// nsub_fixed > 0: do NOT pick the substep count automatically, use the given value. The
// numerical consistent tangent (Pérez-Foguet & Rodríguez-Ferran & Huerta) pins the
// perturbed runs to the base run's nsub: otherwise the ceil(...) jumps n→n±1 under
// perturbation and the integration-error difference breaks the FD column (larger than the
// ~E·h signal).
inline HsIntegrated hs_integrate(const HardeningSoilParams& p,
                                 const Eigen::Vector3d& sigma_n, double gamma_p_n,
                                 double pp_n, const Eigen::Vector3d& dstrain,
                                 int nsub_fixed = 0) {
    const double pr = p.p_ref, plim = 0.1 * pr;  // p_limit (the Eq 15.3 safeguard)
    const double s3_stiff = std::max(sigma_n(2), plim);
    const double Eur = p.Eur(s3_stiff), nu = p.nu_ur;
    const double Ei = p.Ei(s3_stiff), qa = p.q_asymptote(s3_stiff), qf = p.q_failure(s3_stiff);
    const Eigen::Matrix3d De = detail::hs_elastic(Eur, nu);
    const detail::HsDilatancy dil = detail::hs_dilatancy(p);
    const double alpha = p.cap_alpha;
    // The cap deviatoric measure: symmetric von Mises q (q²=3J2). This is EQUIVALENT to
    // the reduced form of PLAXIS's asymmetric q̃ measure (MMM Eq 6-26) on the OEDOMETER
    // (axisymmetric, σ2=σ3) path: the σ2,σ3 components of the q̃ flow average out under
    // axisymmetry → exactly von Mises flow. (Applying q̃ raw to a de2=de3=0 probe produces
    // the σ2≠σ3 absurdity; see hardening-soil-formulation.md §4f, the study_hs_calibration
    // experiment.) f_c=3J2/α²+p²−pc², with 3J2 and p² quadratic ⇒ n_c=H_c·σ LINEAR
    // (closed-form return). H_c=(3/α²)(I−⅓11ᵀ)+(2/9)11ᵀ.
    const Eigen::Matrix3d dev = Eigen::Matrix3d::Identity() -
                                (1.0 / 3.0) * Eigen::Matrix3d::Ones();
    const Eigen::Matrix3d Hc = (3.0 / (alpha * alpha)) * dev +
                               (2.0 / 9.0) * Eigen::Matrix3d::Ones();
    const bool cap_on = p.cap_beta > 0.0;
    const double ev_n = cap_on ? p.cap_ev_from_pc(pp_n) : 0.0;

    auto mean = [](const Eigen::Vector3d& s) { return (s(0) + s(1) + s(2)) / 3.0; };
    auto fbar = [&](double q) { return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur; };
    auto fbar_p = [&](double q) {
        const double r = 1.0 - q / qa; return (2.0 / Ei) / (r * r) - 2.0 / Eur;
    };
    auto fcap = [&](const Eigen::Vector3d& s, double ppv) {
        const double pm = mean(s);
        const double j3 = 0.5 * ((s(0) - s(1)) * (s(0) - s(1)) +
                                 (s(1) - s(2)) * (s(1) - s(2)) +
                                 (s(2) - s(0)) * (s(2) - s(0)));  // 3·J2 (von Mises q²)
        return j3 / (alpha * alpha) + pm * pm - ppv * ppv;
    };
    auto spm_of = [&](double q) { return dil.from_q(q, s3_stiff); };

    // Substep count: split the elastic stress change into ~1%·(|σ|+pref) pieces (robust, cheap).
    const Eigen::Vector3d dsig_e = De * dstrain;
    const double ref_stress = sigma_n.cwiseAbs().maxCoeff() + pr;
    int nsub = static_cast<int>(std::ceil(dsig_e.norm() / (0.01 * ref_stress)));
    nsub = std::min(std::max(nsub, 1), 2000);
    if (nsub_fixed > 0) nsub = std::min(nsub_fixed, 2000);  // FD perturbation: same as the base run
    const Eigen::Vector3d de = dstrain / nsub;

    Eigen::Vector3d sig = sigma_n;
    double gp = gamma_p_n, ev = ev_n;
    auto pp_of = [&]() { return cap_on ? p.cap_pc_from_ev(ev) : pp_n; };
    bool any_plastic = false;
    Eigen::Matrix3d tangent = De;  // the last substep's continuum tangent

    for (int sub = 0; sub < nsub; ++sub) {
        const Eigen::Vector3d sig_tr = sig + De * de;
        const double pp = pp_of();
        const double q_tr = sig_tr(0) - sig_tr(2);
        bool as = fbar(q_tr) - gp > 1e-12 * (1.0 + std::fabs(gp));
        bool ac = cap_on && fcap(sig_tr, pp) > 1e-10 * (1.0 + pp * pp);
        if (!as && !ac) { sig = sig_tr; tangent = De; continue; }
        any_plastic = true;

        // Gradients/flow/hardening of the active surfaces (at the current σ).
        const double q = sig(0) - sig(2);
        // Flow direction m_g=(1,R,R), R=ε3^p/ε1^p. Rowe dilatancy: ε_v^p/ε_q^p=−sinψ_m
        // (DILATION, comp-pos ⇒ ε_v^p<0). R=−(1+sinψ_m)/(2−sinψ_m) gives it (−½ at ψ_m=0 =
        // volumetrically neutral). [The earlier −(1−sinψ_m)/(2+sinψ_m) gave
        // ε_v^p/ε_q^p=+sinψ_m = CONTRACTION — a sign error; never caught because the
        // volumetric response was never verified, see test_hs_berlin Fig 15.4.]
        const double spm = spm_of(q), R = -(1.0 + spm) / (2.0 - spm);
        // Hardening modulus h_s = 1−2R = (4+sinψ_m)/(2−sinψ_m) (Eq 6-10: γ^p=−(2ε1^p−ε_v^p)); 2 at ψ_m=0.
        const double h_s = (4.0 + spm) / (2.0 - spm);
        const Eigen::Vector3d n_s(1.0, R, R);  // flow direction (dilatant)
        // Failure plateau: at q≥qf the shear becomes PERFECTLY-PLASTIC MC (yield f=q−qf,
        // grad (1,0,−1), hardening 0); the flow (1,R,R) stays dilatant → dilation continues
        // along the plateau (Fig 15.4). Otherwise the hardening surface cannot cross qf,
        // stays dormant, and the cap contraction dominates.
        // (PLAXIS: hardening shear → perfectly-plastic MC at failure.)
        const bool at_fail = q >= qf - 1e-9 * (1.0 + qf);
        const Eigen::Vector3d m_s = at_fail ? Eigen::Vector3d(1.0, 0.0, -1.0)
                                            : (fbar_p(q) * Eigen::Vector3d(1.0, 0.0, -1.0));
        const double Hh_s = at_fail ? 0.0 : h_s;
        const Eigen::Vector3d n_c = Hc * sig, m_c = n_c;  // cap associated
        const double pmean = mean(sig);
        const double Hcap = cap_on ? p.cap_hardening_modulus(pp) : 0.0;

        // Active-set Koiter solve (negative multiplier → drop the surface, solve again).
        Eigen::Vector3d dsig;
        for (int pass = 0; pass < 3; ++pass) {
            const int na = (as ? 1 : 0) + (ac ? 1 : 0);
            if (na == 0) { dsig = De * de; break; }
            Eigen::MatrixXd A(na, na); Eigen::VectorXd b(na);
            std::vector<Eigen::Vector3d> ns, ms;
            std::vector<double> Hh;
            if (as) { ns.push_back(n_s); ms.push_back(m_s); Hh.push_back(Hh_s); }
            if (ac) { ns.push_back(n_c); ms.push_back(m_c);
                      Hh.push_back(4.0 * pp * pmean * Hcap); }
            for (int i = 0; i < na; ++i) {
                b(i) = ms[i].dot(De * de);
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
            dsig = De * (de - plastic_strain);
            if (as && !at_fail) gp += h_s * dl_s;  // plateau: hardening freezes (perfectly plastic)
            if (ac) ev += dl_c * 2.0 * pmean;
            // continuum tangent (last substep, active set)
            Eigen::Matrix3d Dep = De;
            {
                const int m = na;
                Eigen::MatrixXd N(3, m), M(3, m);
                int c = 0;
                if (as) { N.col(c) = n_s; M.col(c) = m_s; ++c; }
                if (ac) { N.col(c) = n_c; M.col(c) = m_c; ++c; }
                Dep = De - De * N * A.fullPivLu().solve(M.transpose() * De);
            }
            tangent = Dep;
            break;
        }
        sig += dsig;

        // Drift correction: pull the stress back onto the active surfaces. SHEAR: a stress
        // violating the surface is pulled along the CONSISTENT (Potts & Gens 1985)
        // elastoplastic direction De·n_s — δσ=−(f/(aᵀDe·n_s+h_s))·De·n_s, a=∂f/∂σ. De·n_s is
        // SYMMETRIC in σ2,σ3 (n_s=(1,R,R)) ⇒ **σ2=σ3 IS PRESERVED** (the axisymmetric/
        // oedometer edge). The bare yield gradient (1,0,−1) pushes only σ1,σ3 and leaves σ2
        // → a σ_r≠σ_θ drift in triaxial (the axisym K0 error). The CAP is associated
        // (n_c=Hc·σ symmetric) → the gradient projection never disturbs σ2,σ3 anyway.
        for (int it = 0; it < 5; ++it) {
            const double ppc = pp_of();
            const double qd = sig(0) - sig(2);
            double fs = at_fail ? (qd - qf) : (fbar(qd) - gp);
            double fcp = cap_on ? fcap(sig, ppc) : -1.0;
            bool corr = false;
            if (as && fs > 1e-9 * (1.0 + std::fabs(gp) + qf)) {
                const Eigen::Vector3d a = (at_fail ? 1.0 : fbar_p(qd)) *
                                          Eigen::Vector3d(1.0, 0.0, -1.0);  // ∂f/∂σ
                const Eigen::Vector3d Den = De * n_s;                       // flow direction
                const double denom = a.dot(Den) + (at_fail ? 0.0 : h_s);
                const double dlam_d = fs / denom;
                sig -= dlam_d * Den;
                if (!at_fail) gp += h_s * dlam_d;  // hardening tracks the drift step
                corr = true;
            }
            if (ac && std::fabs(fcp) > 1e-8 * (1.0 + ppc * ppc)) {
                const Eigen::Vector3d g = Hc * sig;
                sig -= (fcp / g.squaredNorm()) * g; corr = true;
            }
            if (!corr) break;
        }
        // MC failure bound (shear): q ≤ qf.
        const double qd = sig(0) - sig(2);
        if (qd > qf) sig(0) -= (qd - qf);
    }

    HsIntegrated out;
    out.stress = sig;
    out.gamma_p = gp;
    out.pp = pp_of();
    out.tangent = tangent;
    out.plastic = any_plastic;
    out.nsub = nsub;
    return out;
}

// Cap preconsolidation initialization (FE initial state): pp = p_eq · OCR,
// p_eq = √(3J2/α² + p²) (the isotropic-equivalent pressure at which the cap passes through
// σ0, the PLAXIS MMM state param p_eq, §6.6). NC (OCR=1) ⇒ the initial state sits on the
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
inline void hs_oedometer_probe(const HardeningSoilParams& p, double& Eoed_pref,
                               double& K0) {
    const double pr = p.p_ref, p0 = 0.02 * pr;
    Eigen::Vector3d sig(p0, p0, p0);
    double gp = 0.0, pp = p0, s1b = p0;
    const double de1 = 2.0e-4;
    Eoed_pref = 0.0; K0 = 0.0;
    for (int i = 0; i < 3000; ++i) {
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, 0, 0));
        if (s1b < pr && r.stress(0) >= pr) {
            Eoed_pref = (r.stress(0) - s1b) / de1;
            K0 = r.stress(2) / r.stress(0);
        }
        s1b = r.stress(0); sig = r.stress; gp = r.gamma_p; pp = r.pp;
        if (sig(0) > 1.25 * pr) break;
    }
}

// Cap hardening modulus K_p = K1·K2/(K1−K2) — CLOSED FORM from Eoed_ref (Itasca
// Plastic-Hardening = the PLAXIS-HS equivalent). K1=Eur_ref/(3(1−2ν)) unloading bulk;
// K2=Eoed_ref(1+2K0nc)/3. This ties the cap hardening directly to Eoed (no numerical
// calibration of β needed) ⇒ the K0-Eoed coupling is resolved: β=p_ref/(k·K_p), only α
// (and a small k correction) is calibrated.
inline double hs_cap_Kp(const HardeningSoilParams& p, double K0_NC) {
    const double K1 = p.Eur_ref / (3.0 * (1.0 - 2.0 * p.nu_ur));
    const double K2 = p.Eoed_ref * (1.0 + 2.0 * K0_NC) / 3.0;
    return (K1 > K2) ? K1 * K2 / (K1 - K2) : K1;  // K1>K2 (Eur≫Eoed) typical
}

// Cap parameters (α, β) from the PLAXIS-standard inputs. β = p_ref/(k·K_p) (K_p closed
// form, sets Eoed in CLOSED FORM → the coupling is resolved); α by outer bisection for
// K0_NC; k by inner bisection as an Eoed_ref fine correction (k≈1). The Itasca PH / PLAXIS
// HS method. If K0_NC is unreachable, the nearest α.
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
