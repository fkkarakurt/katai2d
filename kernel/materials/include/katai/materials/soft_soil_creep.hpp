#pragma once
// SOFT SOIL CREEP (SSC) — STAGE 1: the material-point core (PLAXIS MMM §11 verbatim; locked
// formulation docs/references/soft-soil-creep-formulation.md). Time-dependent
// viscoplasticity: creep = time-dependent plastic strain, potential g = p_eq (the MCC
// ellipse, IDENTICAL to the SS f̄ measure), volumetric rate ε̇_v^c = (μ*/τ)(p_eq/p_p)^β,
// β = (λ*−κ*)/μ*, p_p ages exponentially; failure is a SEPARATE Mohr-Coulomb check AFTER
// the creep update (the manual §11.7 order).
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
// (the FE stage does it with MMM §17.4). Sign convention compression-POSITIVE. FE/GUI =
// the next stage.

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
    double tau_day = 1.0;     // reference time τ [days] — the 24-hour definition of the NC line (Eq 11-14)

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
// (Eq 11-20; by the single-source reading the formula matches ss_initial_pp/f_cap).
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
inline StepResult ssc_substep(const Params& P, const softsoil::Params& S,
                              const Eigen::Vector3d& sig_c, double pp_c,
                              const Eigen::Vector3d& deps, double dt) {
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi);
    const double ccot = P.phi > 1e-12 ? P.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = softsoil::M_from_K0nc(S);
    const double beta = (P.lam_star - P.kap_star) / P.mu_star;
    const double ppfloor = std::max(ccot, softsoil::kPmin);

    // --- The same exponentially exact elastic predictor as SS ---------------------------------
    const double p_c = std::max(sig_c.mean(), softsoil::kPmin);
    const double dv = deps.sum();
    const double p_tr = std::max(p_c * std::exp(dv / P.kap_star), softsoil::kPmin);
    const double K_tr = p_tr / P.kap_star;
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const Eigen::Vector3d de = deps - Eigen::Vector3d::Constant(dv / 3.0);
    const Eigen::Vector3d s_dev_c = sig_c - Eigen::Vector3d::Constant(sig_c.mean());
    const Eigen::Vector3d sig_tr =
        Eigen::Vector3d::Constant(p_tr) + s_dev_c + 2.0 * G * de;

    StepResult out{sig_tr, std::max(pp_c, ppfloor), 0.0, false, 1};

    // --- Creep correction: L = Δε_v^c backward Euler, direction frozen at the trial -----------
    // Sorted-role analytic direction n = a·w + b·1 (the SS cap_normal structure) +
    // VALIDITY cascade: if the creep step breaks the trial ordering, the tied pair's roles
    // are averaged and the pair equalized (with w2=w3, p_eq depends only on σ2+σ3 → the
    // equalization changes neither p_eq nor the hardening; the creep counterpart of the
    // corner-branching lesson measured on SS — tolerance-free, scale-independent).
    const double pe_tr = p_eq(S, sig_tr);
    if (pe_tr > 1e-12 && dt > 0.0) {
        int idx[3] = {0, 1, 2};
        std::sort(idx, idx + 3, [&](int a, int b) { return sig_tr(a) > sig_tr(b); });
        bool tied01 = false, tied12 = false;
        for (int attempt = 0; attempt < 3; ++attempt) {
            Eigen::Vector3d w(1.0, delta - 1.0, -delta);
            if (tied01 && tied12) {
                w.setConstant(0.0);   // full tie: no q̃ contribution (purely volumetric)
            } else if (tied01) {
                const double mw = 0.5 * (w(0) + w(1)); w(0) = mw; w(1) = mw;
            } else if (tied12) {
                const double mw = 0.5 * (w(1) + w(2)); w(1) = mw; w(2) = mw;
            }
            // n = a·w + b·1, tr(n) = 3b (w sums to 0). m = n/tr(n): tr(m) = 1 → the
            // volumetric component is EXACTLY the 1D rate (the Eq 11-30 normalization).
            // Around the nose, |3b| is protected by the floor.
            const double pb = std::max(sig_tr.mean() + ccot, 1e-9);
            const double qt = sig_tr(idx[0]) + (delta - 1.0) * sig_tr(idx[1]) - delta * sig_tr(idx[2]);
            const double a = 2.0 * qt / (M * M * pb);
            double b3 = 1.0 - qt * qt / (M * M * pb * pb);
            if (std::fabs(b3) < 1e-6) b3 = b3 < 0.0 ? -1e-6 : 1e-6;
            Eigen::Vector3d n;
            for (int k = 0; k < 3; ++k) n(idx[k]) = a * w(k) + b3 / 3.0;
            const Eigen::Vector3d m = n / b3;                         // tr(m) = 1
            const Eigen::Vector3d mdev = m - Eigen::Vector3d::Constant(1.0 / 3.0);

            auto state_at = [&](double L, Eigen::Vector3d& sL, double& ppL) {
                const double pL = std::max(p_tr * std::exp(-L / P.kap_star), softsoil::kPmin);
                sL = Eigen::Vector3d::Constant(pL) +
                     (sig_tr - Eigen::Vector3d::Constant(sig_tr.mean())) - L * 2.0 * G * mdev;
                ppL = std::max(out.pp * std::exp(L / (P.lam_star - P.kap_star)), ppfloor);
            };
            auto resid = [&](double L, Eigen::Vector3d& sL, double& ppL) {
                state_at(L, sL, ppL);
                const double pe = std::max(p_eq(S, sL), 1e-12);
                return L - dt * (P.mu_star / P.tau_day) * std::pow(pe / ppL, beta);
            };
            // Bracketed root: R(0) ≤ 0, R monotone increasing (L↑ ⇒ p_eq↓, p_p↑) → guaranteed convergence.
            Eigen::Vector3d sL; double ppL;
            double lo = 0.0, flo = resid(0.0, sL, ppL);
            double L = 0.0;
            if (flo < 0.0) {
                double hi = std::max(1e-6 * P.kap_star, -flo);
                Eigen::Vector3d sh; double ph;
                double fhi = resid(hi, sh, ph);
                for (int it = 0; it < 200 && fhi < 0.0; ++it) { hi *= 2.0; fhi = resid(hi, sh, ph); }
                for (int it = 0; it < 100 && (hi - lo) > 1e-16 + 1e-12 * hi; ++it) {
                    // secant candidate + bisection safeguard
                    double mid = 0.5 * (lo + hi);
                    if (fhi > flo) {
                        const double sec = lo - flo * (hi - lo) / (fhi - flo);
                        if (sec > lo && sec < hi) mid = sec;
                    }
                    Eigen::Vector3d sm; double pm;
                    const double fm = resid(mid, sm, pm);
                    if (fm < 0.0) { lo = mid; flo = fm; }
                    else { hi = mid; fhi = fm; }
                }
                L = hi;
                state_at(L, sL, ppL);
            } else {
                state_at(0.0, sL, ppL);
            }
            // equalization (tied pair → common value; p_eq/hardening invariant)
            if (tied01 && tied12) {
                sL.setConstant(sL.mean());
            } else if (tied01) {
                const double ms = 0.5 * (sL(idx[0]) + sL(idx[1]));
                sL(idx[0]) = ms; sL(idx[1]) = ms;
            } else if (tied12) {
                const double ms = 0.5 * (sL(idx[1]) + sL(idx[2]));
                sL(idx[1]) = ms; sL(idx[2]) = ms;
            }
            const double slack = 1e-9 * (1.0 + std::fabs(sL(idx[0])));
            const bool bad01 = !tied01 && sL(idx[0]) - sL(idx[1]) < -slack;
            const bool bad12 = !tied12 && sL(idx[1]) - sL(idx[2]) < -slack;
            if (!bad01 && !bad12) {
                out.sig = sL;
                out.pp = ppL;
                out.devc = L;
                break;
            }
            tied01 = tied01 || bad01;
            tied12 = tied12 || bad12;
        }
    }

    // --- Mohr-Coulomb check (manual order: AFTER creep) — single source: softsoil::ss_step
    // with pp=∞ (cap off, deps=0 → the elastic predictor is the identity, only the MC
    // edge-cascaded return).
    {
        const softsoil::StepResult r =
            softsoil::ss_step(S, out.sig, 1.0e300, Eigen::Vector3d::Zero());
        if (r.mc_active) {
            out.sig = r.sig;
            out.mc_active = true;
        }
    }
    return out;
}

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
inline StepResult ssc_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                           const Eigen::Vector3d& deps, double dt, int nsub_fixed = 0) {
    const softsoil::Params S = P.ss();
    const int nsub = nsub_fixed > 0 ? nsub_fixed : ssc_nsub(P, S, sig_c, pp_c, deps, dt);
    const Eigen::Vector3d dsub = deps / nsub;
    const double dtsub = dt / nsub;
    StepResult out{sig_c, std::max(pp_c, softsoil::kPmin), 0.0, false, nsub};
    for (int i = 0; i < nsub; ++i) {
        const StepResult r = ssc_substep(P, S, out.sig, out.pp, dsub, dtsub);
        out.sig = r.sig;
        out.pp = r.pp;
        out.devc += r.devc;
        out.mc_active = out.mc_active || r.mc_active;
    }
    return out;
}

}  // namespace katai::core::softsoilcreep
