#include <cstdlib>
#include <exception>
#include <string>
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/interrupt.hpp>
#include <katai/materials/mohr_coulomb.hpp>   // apply_rankine_cap: the tension cut-off's own return

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <utility>
#include <vector>

namespace katai::core {

HsShearStep hs_shear_step(const HardeningSoilParams& p,
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

HsPrincipalReturn hs_shear_correct(const HardeningSoilParams& p,
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

HsCapStep hs_cap_step(const HardeningSoilParams& p,
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

HsState hs_return_principal(const HardeningSoilParams& p,
                            const Eigen::Vector3d& sigma_n, double gamma_p_n,
                            double pp_n, const Eigen::Vector3d& dstrain,
                            double sigma3_stiff) {
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

HsIntegrated hs_integrate(const HardeningSoilParams& p_in,
                          const Eigen::Vector3d& sigma_n, double gamma_p_n,
                          double pp_n, const Eigen::Vector3d& dstrain,
                          double stol,
                          HsSubstepPlan* plan_out,
                          const HsSubstepPlan* plan_in) {
    // This call's own copy, with the cap's c.cot(phi) worked out once (HardeningSoilParams::cap_cc).
    HardeningSoilParams p = p_in;
    p.cap_cc = p_in.cap_shift();
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
        double s3, s3raw, Eur, Ei, qa, qf;
        Eigen::Matrix3d De;
    };
    // The stiffness and strength laws read the MINOR principal stress (Benz 2007, Eqn 7.19), and
    // the integrator's vector is not assumed sorted: see Order below.
    // The laws of HardeningSoilParams (stiffness_factor, E50/Eur/Ei, q_failure/q_asymptote)
    // written out with their trigonometry evaluated once per call of the integrator rather than
    // six times per evaluation -- the same operations in the same order, so the same bits.
    const double law_cc = p.cohesion * std::cos(p.friction);
    const double law_s = std::sin(p.friction);
    const double law_den = law_cc + p.p_ref * law_s;
    // THE FLOOR IS A STIFFNESS FLOOR. p_limit exists so that E_i and E_ur do not fall to zero
    // (and the hyperbola's 2/E_i blow up) at a stress-free point; it is not part of the failure
    // criterion, which is Mohr-Coulomb in the stress itself. Until 2026-09-29 the strength read the
    // floored stress too, so a cohesionless point shallower than p_limit carried q_f(p_limit) --
    // on a sand of phi = 30 degrees a deviator of 20 kPa at zero confinement, an apparent cohesion
    // of 5.8 kPa, on the unsafe side and largest exactly where a footing's bearing capacity is
    // decided. q_f, q_a = q_f/R_f and the mobilised friction angle now read the minor stress
    // itself. (Reading q_a at the floored stress instead -- a hyperbola far flatter than the
    // strength it ends at -- put the shallow points of a column loaded from zero on the failure
    // plateau after almost no hardening, where they flipped between plateau and hardening and held
    // the equilibrium iteration at 1e-3.) Only the stress-free point itself needs guarding, where
    // q/q_a is 0/0: q_a is kept above its value at a thousandth of the floor.
    //
    // AND ONLY WHERE IT IS 0/0: without cohesion. With c > 0 the strength at zero stress is
    // 2 c cos(phi)/(1 - sin(phi)) and nothing divides by zero until the apex, -c cot(phi); the
    // guard, read at sigma3 = p_limit/1000 regardless, was a kink in the flow there (q_a and its
    // slope switched at sigma3 = 0.01 kPa) that no substep of a walk across it could resolve --
    // measured on a point pulled to the tension cut-off with c = 10 kPa, the one part of its walk
    // that crossed 0.01 kPa held R at 1.2e-5 at every refinement, and the walk to the ceiling.
    // There it stands a hair above the apex, where the flow is on the failure branch anyway.
    const double qf_min = law_cc > 0.0 ? 2.0e-9 * law_cc / (1.0 - law_s)
                                       : 2.0 * (1e-3 * plim * law_s) / (1.0 - law_s);
    auto stiff_compute = [&](const Eigen::Vector3d& s) {
        Stiff k;
        k.s3raw = s.minCoeff();
        k.s3 = std::max(k.s3raw, plim);
        const double g = std::pow((law_cc + k.s3 * law_s) / law_den, p.m);
        k.Eur = p.Eur_ref * g;
        k.Ei = 2.0 * (p.E50_ref * g) / (2.0 - p.Rf) * p.ei_scale;
        k.qf = 2.0 * (law_cc + k.s3raw * law_s) / (1.0 - law_s);
        k.qa = std::max(k.qf, qf_min) / p.Rf;
        // The elastic matrix carries the small-strain multiplier; E_ur in the surface does not
        // (see HardeningSoilParams::elastic_scale).
        k.De = detail::hs_elastic(k.Eur * p.elastic_scale, nu);
        return k;
    };
    // The same stress is asked for its stiffness several times in a row: one increment() reads the
    // substep's start three times and its elastic trial twice (the onset test, the shear trial).
    // Two remembered states, matched by EXACT equality -- the same input gives the same bits, so
    // nothing the integration computes can change; only the repeated pow() and De are saved. Local
    // to this call, so each stress point has its own.
    Eigen::Vector3d memo_s[2];
    Stiff memo_k[2];
    int memo_n = 0, memo_next = 0;
    auto stiff_at = [&](const Eigen::Vector3d& s) -> Stiff {
        for (int i = 0; i < memo_n; ++i)
            if (memo_s[i] == s) return memo_k[i];
        const Stiff k = stiff_compute(s);
        memo_s[memo_next] = s;
        memo_k[memo_next] = k;
        memo_next = 1 - memo_next;
        if (memo_n < 2) ++memo_n;
        return k;
    };
    const Stiff k_n = stiff_at(sigma_n);   // the committed state, for the elastic default
    const bool cap_on = p.cap_beta > 0.0;
    // THE TENSION CUT-OFF (HardeningSoilParams::sigma_t), in this compression-positive frame
    // s_i >= -ten: three planes more for the same return -- the Rankine planes sigma_i - sigma_t <= 0
    // with an associated flow rule and sharp corners per Koiter, as the Mohr-Coulomb return of this
    // tree takes them (mohr_coulomb.hpp). Until 2026-10-29 the wrapper applied it after this return, and
    // with c = 0 the capped state lay OUTSIDE the shear cone -- the minor stress cut to zero, the
    // major not -- so the next increment started from an inadmissible state: measured on a staged
    // excavation, 147 such points near the apex took 466 million substep evaluations of 537.
    const bool ten_on = p.sigma_t < 1e299;
    const double ten = ten_on ? p.sigma_t : 0.0;
    double ev_n = cap_on ? p.cap_ev_from_pc(pp_n) : 0.0;

    // WHICH COMPONENT IS WHICH PRINCIPAL STRESS is decided by the stress being evaluated, every
    // time. The vector this integrator walks keeps the positions its caller gave it -- one per
    // principal DIRECTION -- and the order of the VALUES in it can change inside one increment:
    // in triaxial extension the axial stress passes from major to minor. Until 2026-09 position
    // 0 was taken to be sigma1 and position 2 sigma3 for the whole increment, so once the values
    // had swapped the flow was applied to the wrong directions. Measured on an undrained
    // extension of a normally consolidated clay: the two lateral stresses leap-frogged each
    // other, the mean effective stress DOUBLED on a psi = 0 path and the deviator came out 2 to
    // 2.6 times the reference.
    struct Order { int a, b, c; };   // s(a) >= s(b) >= s(c): major, intermediate, minor
    auto order_of = [](const Eigen::Vector3d& s) {
        int i0 = 0, i1 = 1, i2 = 2;
        if (s(i0) < s(i1)) std::swap(i0, i1);
        if (s(i1) < s(i2)) std::swap(i1, i2);
        if (s(i0) < s(i1)) std::swap(i0, i1);
        return Order{i0, i1, i2};
    };

    // THE CAP MEASURE (Schanz 1998, Eqns 4.71-4.74; Benz 2007, Eqn 7.21): q~ = sigma1 + (delta - 1) sigma2 - delta sigma3,
    // delta = (3 + sin phi)/(3 - sin phi). On the triaxial compression corner q~ = sigma1 - sigma3,
    // which is the von Mises q this cap used to be written with -- so every oedometer and
    // compression result, and the (alpha, beta) calibration done on that path, is unchanged by
    // construction. In extension q~ = delta (sigma1 - sigma3): the cap is met earlier there,
    // which the von Mises form did not do at all.
    const double sphi_f = std::sin(p.friction);
    const double delta = (3.0 + sphi_f) / (3.0 - sphi_f);
    auto mean = [](const Eigen::Vector3d& s) { return (s(0) + s(1) + s(2)) / 3.0; };
    auto qtilde = [&](const Eigen::Vector3d& s) {
        const Order o = order_of(s);
        return s(o.a) + (delta - 1.0) * s(o.b) - delta * s(o.c);
    };
    // THE CAP CLOSES THE COMPRESSIVE SIDE. It is the second yield surface of the model because the
    // shear surfaces leave the elastic region open along the p axis for compaction (Schanz, Vermeer &
    // Bonnier 1999); on the tensile side the region is closed by the Mohr-Coulomb cone and, where it
    // is switched on, the tension cut-off. The ellipse written out is symmetric in p, and read as a
    // yield surface at p <= 0 it made every stress in tension a cap state as soon as p_p was small:
    // with p_p = 0 -- the ground surface, fresh fill -- a point pulled a ten-thousandth of a kPa into
    // tension flowed on a "cap" there, its pair of evaluations disagreeing by the whole step at any
    // substep size (measured on a strip footing on c = 5 kPa sand: 709 walks near the tip at
    // 7 300 substeps each, 5.2 million of the run's 8.2 million, and a third more equilibrium
    // iterations; 68 s -> 8 s without it). At p <= 0 it is no yield surface: inside, whatever q~.
    auto fcap = [&](const Eigen::Vector3d& s, double ppv) {
        const double qt = qtilde(s), pm = mean(s);
        if (pm <= 0.0) return -(1.0 + ppv * ppv);
        return qt * qt / (alpha * alpha) + pm * pm - ppv * ppv;
    };
    // Where the stress sits relative to the corners of the hexagon. 0: a face (sigma1 > sigma2 >
    // sigma3); 1: the compression corner (sigma2 = sigma3); 2: the extension corner (sigma1 =
    // sigma2). A corner is not a special case of the face but the meeting of two of them, and
    // its gradient is the equal-weight combination of the two faces' -- the value Koiter's rule
    // gives when the loading is symmetric, and the one that keeps a symmetric state symmetric.
    enum Corner { kFace = 0, kCompression = 1, kExtension = 2 };
    // Gradient of the cap function: 2 q~/alpha^2 dq~/dsigma + (2p/3) 1.
    auto cap_grad = [&](const Eigen::Vector3d& s, const Order& o, int corner) {
        const double qt = qtilde(s), pm = mean(s);
        Eigen::Vector3d dq = Eigen::Vector3d::Zero();
        if (corner == kCompression) {
            dq(o.a) = 1.0; dq(o.b) = -0.5; dq(o.c) = -0.5;
        } else if (corner == kExtension) {
            dq(o.a) = 0.5 * delta; dq(o.b) = 0.5 * delta; dq(o.c) = -delta;
        } else {
            dq(o.a) = 1.0; dq(o.b) = delta - 1.0; dq(o.c) = -delta;
        }
        return Eigen::Vector3d(2.0 * qt / (alpha * alpha) * dq +
                               (2.0 * pm / 3.0) * Eigen::Vector3d::Ones());
    };
    auto fbar = [&](double q, const Stiff& k) {
        return (2.0 / k.Ei) * q / (1.0 - q / k.qa) - 2.0 * q / k.Eur;
    };
    auto fbar_p = [&](double q, const Stiff& k) {
        const double r = 1.0 - q / k.qa; return (2.0 / k.Ei) / (r * r) - 2.0 / k.Eur;
    };
    auto spm_of = [&](double q, const Stiff& k) { return dil.from_q(q, k.s3raw); };
    // Beyond the Mohr-Coulomb apex (q_f <= 0: the minor stress in tension past c cot phi) no
    // deviator is admissible and the hyperbola has no meaning; the failure branch takes over.
    auto beyond = [&](double q, const Stiff& k) { return q >= k.qf; };
    // The shear yield function OF A STRESS: f = fbar(q; E_i, q_a, E_ur) - gamma_p, with E_i, q_a and
    // E_ur read at that stress's own minor principal -- they are part of the surface, so the
    // surface is a function of the stress alone. Scaled to be dimensionless, continuous across the
    // failure deviator (past it, the Mohr-Coulomb excess).
    //
    // Until 2026-09-29 the ACTIVATION tests read the surface with the stiffness of the state the
    // substep starts from, while the flow (whose gradient carries dfbar/dsigma3) and the drift
    // correction read it at the stress itself. On any path that lowers sigma3 -- the ground next to
    // a footing -- the deviator falls while the surface shrinks faster: the frozen test said
    // "elastic" at the start of every substep and the Euler end point said "plastic", the two
    // evaluations of the modified Euler pair disagreed at every size, the error estimate never
    // shrank faster than the substep, and 29% of the integrations of a footing on sand ran to the
    // 200-substep ceiling with a max stress error of 1e-3.
    //
    // The two parts are joined by a MAX, not switched between at q = q_f: the hyperbola's value there
    // is generally positive (the surface has not hardened to failure yet) and the excess is zero, so
    // a switch jumps -- and a trial just past the failure deviator, far outside the hardening
    // surface, read as barely on it. Measured on the corpus oedometer, whose column starts stress
    // free: +9.2% on the loading settlement and four times the run time, back to -1.0% with the max.
    auto fshear_at = [&](const Eigen::Vector3d& s, double gpv) {
        const Stiff ks = stiff_at(s);
        const double qq = s.maxCoeff() - s.minCoeff();
        const double excess = (qq - ks.qf) / (1.0 + std::fabs(ks.qf));
        if (!(ks.qf > 0.0)) return excess;       // past the apex: no hyperbola to read
        return std::max((fbar(std::min(qq, ks.qf), ks) - gpv) / (1.0 + std::fabs(gpv)), excess);
    };
    auto pc_of = [&](double evv) { return cap_on ? p.cap_pc_from_ev(evv) : pp_n; };
    // AN INADMISSIBLE START IS SEEDED, NOT RETURNED. A committed state whose hardening variables lag
    // its own stress -- a cluster given this material in a later phase keeps its stress and starts
    // with gamma_p = 0, a stress field entered by hand, a history seeded elsewhere -- lies outside
    // its shear surface or its cap before any strain is applied. Walked from there, the first
    // evaluation of a pair pulled the STRESS back to the stale surface and the second did not, by
    // the same amount at every substep size: the error estimate never fell and every such call ran
    // to the ceiling (measured: every call of the HSsmall tangent check, 2000 random states with
    // gamma_p = 0, R 1e-2 at 10 000 substeps), and the stress it returned had moved with no strain.
    // The hardening variables are raised to the stress instead, as the initial state is seeded
    // (hs_seed_from_history: the shear hardening that puts the present stress on its surface, the
    // cap at least the present cap measure): continuous in the state, and nothing for a state that
    // is already inside.
    {
        const Stiff ks = stiff_at(sigma_n);
        const double qn = sigma_n.maxCoeff() - sigma_n.minCoeff();
        if (ks.qf > 0.0 && qn > 0.0) gamma_p_n = std::max(gamma_p_n, fbar(std::min(qn, ks.qf), ks));
        const double pm = mean(sigma_n);
        if (cap_on && pm > 0.0) {
            const double qt = qtilde(sigma_n);
            const double pc_req = std::sqrt(qt * qt / (alpha * alpha) + pm * pm);
            if (pc_req > pp_n) { pp_n = pc_req; ev_n = p.cap_ev_from_pc(pp_n); }
        }
    }
    // The Mohr-Coulomb failure function of the plane (i, j), i the larger: f = s_i - s_j - q_f(s_j),
    // q_f(s) = (c cot phi + s) 2 sin phi/(1 - sin phi). Its gradient carries the dependence of the
    // strength on the smaller stress; the (1, 0, -1) this used to be written with does not, so
    // the multiplier it gave left the stress off the failure surface and a clamp had to finish
    // the job.
    const double dqf_ds = 2.0 * sphi_f / (1.0 - sphi_f);

    // A SHEAR PLANE (i, j) of the hexagon (Benz 2007, Eqns 7.15-7.16, 7.22):
    //   yield      f_ij = fbar(s_i - s_j) - gamma_p          (hardening region)
    //              f_ij = s_i - s_j - q_f                     (on the failure plateau)
    //   potential  g_ij = (s_i - s_j)/2 - (s_i + s_j)/2 sin psi_m   (compression positive)
    //   hardening  d gamma_p = d lambda (dg/ds1 - dg/ds2 - dg/ds3 = 1 on every plane)
    // so that d eps_v^p = sin psi_m d gamma_p (Schanz 1998, Eqn 4.40) holds on a face and at either corner
    // alike. The flow used to be (1, R, R) with R = -(1 + sin psi_m)/(2 - sin psi_m): the shape of
    // the compression corner only, with a dilatancy of -3 sin psi_m/(4 + sin psi_m) per unit
    // gamma_p instead of -sin psi_m -- 27% short at psi_m = 5.5 degrees -- and a plastic strain
    // along the intermediate stress on every face, where the potential has none.
    struct Plane { Eigen::Vector3d n, m; };
    // `floor_side` (16 above, 32 below, 0 as the stress has it): which side of the stiffness floor the
    // gradient is read on, for a stress ON the floor (see the floor event in substep()).
    auto plane = [&](const Eigen::Vector3d& s, int i, int j, double spm, bool at_fail,
                     const Stiff& k, int floor_side = 0) {
        Plane pl;
        pl.n = Eigen::Vector3d::Zero();
        pl.n(i) = 0.5 * (1.0 - spm);
        pl.n(j) = -0.5 * (1.0 + spm);
        pl.m = Eigen::Vector3d::Zero();
        if (at_fail) {
            // The failure strength reads the minor stress itself (see Stiff): d q_f/d sigma_j at
            // every depth, shallow or not.
            pl.m(i) = 1.0;
            pl.m(j) = -1.0 - dqf_ds;
        } else {
            const double qij = s(i) - s(j);
            const double fp = fbar_p(qij, k);
            pl.m(i) = fp;
            pl.m(j) = -fp;
            // fbar also depends on the minor principal stress, through E_i, q_a and E_ur. Leaving
            // that out made every hardening step inconsistent to FIRST order: gamma_p ran ahead
            // of fbar, the stress was left inside the plane it had just flowed on, and the answer
            // moved with the outer step (gamma_p by 4% between 4 and 400 steps of the same
            // oedometer). Where two stresses share the minimum the derivative is split between
            // them, so a corner stays symmetric.
            const double s3 = k.s3;
            {
                // In closed form. E_i and E_ur scale with the same stiffness factor g(sigma3),
                // dg/dsigma3 = g m sin(phi)/(c cos(phi) + sigma3 sin(phi)), and q_a = q_f/R_f is
                // linear in sigma3; fbar = A - 2q/E_ur with A = (2/E_i) q/(1 - q/q_a). Below the
                // floor the stiffness no longer moves with sigma3; q_a does, down to its guard.
                const double sphi = std::sin(p.friction);
                const bool above = floor_side == 16 || (floor_side != 32 && s3 > plim);
                const double r = above ? p.m * sphi / (p.cohesion * std::cos(p.friction) + s3 * sphi)
                                       : 0.0;
                const double dqa = k.qf > qf_min ? 2.0 * sphi / ((1.0 - sphi) * p.Rf) : 0.0;
                const double w = 1.0 - qij / k.qa;
                const double A = (2.0 / k.Ei) * qij / w;
                const double dfd3 = -A * r - (2.0 / k.Ei) * qij * qij / (k.qa * k.qa * w * w) * dqa +
                                    2.0 * qij / k.Eur * r;
                const double smin = s.minCoeff();
                const double tie = 1e-9 * (1.0 + std::fabs(smin));
                int nmin = 0;
                for (int x = 0; x < 3; ++x) nmin += (s(x) - smin <= tie) ? 1 : 0;
                for (int x = 0; x < 3; ++x)
                    if (s(x) - smin <= tie) pl.m(x) += dfd3 / nmin;
            }
        }
        return pl;
    };

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
        bool as = false;       // shear active (after Koiter drops)
        bool ac = false;       // cap active
        bool at_fail = false;  // on the perfectly-plastic MC plateau
        bool at = false;       // a tension cut-off plane active
        int corner = kFace;    // which part of the hexagon the step flowed on
        int floor_side = 0;    // the side of the stiffness floor the step read (allow bits 16/32)
        // At a corner, the two principal DIRECTIONS whose stresses are equal (-1: on a face). The
        // corner TYPE only means something relative to the ordering it was found in; the pair is
        // what carries over to a state whose ordering has changed.
        int tie_i = -1, tie_j = -1;
    };
    // `allow` (bit 0 the shear planes, bit 1 the cap, bit 2 a corner of them, bit 3 the tension
    // cut-off) names what this step may activate; a step held to fewer is how the walk finds where
    // the flow meets the rest. Bits 4 and 5 (16, 32) read the shear planes' gradient above or below
    // the stiffness floor, for a step that starts on it.
    auto increment_core = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                              const Eigen::Vector3d& de_in, int allow = 15) -> Inc {
        Inc r;
        r.floor_side = allow & 48;
        const Stiff k = stiff_at(s_in);   // the stiffness and strength THIS substep starts from
        const Eigen::Matrix3d& De = k.De;
        r.dsig = De * de_in;
        r.tangent = De;
        const double pp = pc_of(ev_in);
        const Eigen::Vector3d sig_tr = s_in + De * de_in;
        // A trial at or past the failure deviator is a shear trial, whatever the hyperbola says:
        // past q_a it changes sign, and a trial that far out used to be read as inside.
        const bool as0 = (allow & 1) && fshear_at(sig_tr, gp_in) > 1e-12;
        const bool ac0 = (allow & 2) && cap_on && fcap(sig_tr, pp) > 1e-10 * (1.0 + pp * pp);
        const bool at0 = (allow & 8) && ten_on && -sig_tr.minCoeff() - ten > 0.0;
        if (!as0 && !ac0 && !at0) return r;
        r.plastic = true;

        // Equal values carry no direction, so a tie is broken by the elastic trial: the stress an
        // isotropic state is being loaded towards is the one that becomes the major.
        const Order o = [&] {
            int i0 = 0, i1 = 1, i2 = 2;
            auto above = [&](int i, int j) {
                return s_in(i) > s_in(j) || (s_in(i) == s_in(j) && sig_tr(i) > sig_tr(j));
            };
            if (above(i1, i0)) std::swap(i0, i1);
            if (above(i2, i1)) std::swap(i1, i2);
            if (above(i1, i0)) std::swap(i0, i1);
            return Order{i0, i1, i2};
        }();
        const double q = s_in(o.a) - s_in(o.c);
        const double spm = spm_of(q, k);
        // On the plateau: a positive deviator at the Mohr-Coulomb line. A state with no deviator
        // at all -- a stress-free point, whose q_f is zero too -- is at the apex, not on the
        // plateau, and a compression takes it into the cone by hardening.
        const bool at_fail = q >= k.qf - 1e-9 * (1.0 + std::fabs(k.qf)) &&
                             q > 1e-12 * (1.0 + s_in.cwiseAbs().maxCoeff());
        // d gamma_p per unit multiplier of a plane, and the cap's modulus: both carry the
        // small-strain factor h_i (1 for plain Hardening Soil).
        const double hi = p.hardening_scale;
        const double H = at_fail ? 0.0 : hi;
        const double pmean = mean(s_in);
        const double Hcap = cap_on ? hi * p.cap_hardening_modulus(pp) : 0.0;
        r.at_fail = at_fail;

        // Solve the active set for one corner hypothesis; returns false if nothing is active.
        struct Sol {
            Eigen::Vector3d dsig; double dgp = 0.0, dev = 0.0;
            Eigen::Matrix3d tangent; bool as = false, ac = false, at = false;
            bool face = false;   // a corner hypothesis solved as its face (the cut-off holds the edge)
        };
        auto solve = [&](int corner, Sol& sol) {
            bool as = as0, ac = ac0;
            // The cut-off planes the trial violates, one per principal direction. All three: the
            // stress is the cut-off's hydrostatic apex, which holds no deviator, so the shear planes
            // and the cap have nothing to add (and through that point they would make the system
            // singular).
            std::array<char, 3> ton{};
            for (int i = 0; i < 3; ++i) ton[i] = at0 && -sig_tr(i) - ten > 0.0;
            // A COMPRESSION CORNER WHOSE TWO SMALLER STRESSES ARE BOTH ON THE CUT-OFF is held there by
            // the cut-off: the second shear face and the second cap face only restate it, and with
            // them the system is five surfaces in a three-dimensional stress space -- singular. The
            // fall-back for a singular system then set the shear planes and the cap aside, and a
            // major stress still hardening on both responded elastically instead: the first
            // evaluation of a pair hardened, the second (from the same edge, a hair further) did not,
            // by the whole of the step at every step size (measured: a strip load on c = 1 kPa sand,
            // walks from rest at the ceiling of 10 000 substeps, R 1e-5 on 0.2 kPa states). It is
            // solved as the face, which with the two cut-off planes keeps the edge, and it is the face
            // the step reports having flowed on (Sol::face): the pair's second evaluation is held to
            // what the first one did.
            const bool pinned = corner == kCompression && ton[o.b] && ton[o.c];
            sol.face = pinned;
            // The planes a corner hypothesis brings: the face (a, c), and at a corner the second
            // face meeting it.
            // At most two planes and the cap: fixed-capacity storage on the stack. This runs
            // several times per substep at every plastic stress point, and heap vectors here cost
            // more than the arithmetic (a Berlin dewatering step: 27 000 plastic points).
            std::array<Plane, 2> planes;
            std::array<double, 2> fp{};   // each plane's yield function where the substep starts
            int nplanes = 0;
            auto add = [&](int i, int j) {
                planes[nplanes] = plane(s_in, i, j, spm, at_fail, k, r.floor_side);
                const double qij = s_in(i) - s_in(j);
                // The same failure strength as the plateau test, the drift correction and the
                // Mohr-Coulomb bound use (Stiff::qf, the minor stress itself). They must agree: when
                // this line read the raw stress while the others read a floored one, a shallow
                // point carried two failure surfaces at once and a strip footing on c = 5 kPa sand
                // stopped converging at 99% of its service load.
                fp[nplanes] = at_fail ? qij - p.q_failure(s_in(j)) : fbar(qij, k) - gp_in;
                ++nplanes;
            };
            add(o.a, o.c);
            if (corner == kCompression) add(o.a, o.b);
            if (corner == kExtension) add(o.b, o.c);
            std::array<char, 2> on{1, (char)(pinned ? 0 : 1)};
            // THE CAP AT A CORNER IS TWO FACES TOO. q~ is piecewise linear in the ordered
            // principals, so the cap has an edge wherever two of them are equal, exactly as the
            // shear hexagon does, and it is treated the same way: each face with its own gradient
            // and its own value where the substep starts, both hardening the one pp (Koiter).
            // With both active the linearised consistency puts the end state on both faces, i.e.
            // back on the edge, and a face whose multiplier comes out negative is dropped -- so
            // the response is continuous where the flow joins or leaves the edge.
            //
            // Until 2026-09-29 a corner carried ONE cap surface with the averaged gradient. The
            // face flow of this cap drives the two smaller stresses together (it is how a K0
            // state loaded in compression behaves), and once they met the averaged surface did
            // not hold them there: substep after substep alternated between face and corner, and
            // the stress returned jumped with the strain increment. Measured on a strip footing on
            // normally consolidated sand (every point starts ON that edge, sigma_xx = sigma_zz),
            // 0.2% of the plastic points jumped at every Newton iteration, deep ones included, and
            // the iteration stalled with the line search finding no descent.
            const Order o_b = corner == kCompression ? Order{o.a, o.c, o.b}
                                                     : Order{o.b, o.a, o.c};   // the other face
            const Eigen::Vector3d nc = cap_on ? cap_grad(s_in, o, kFace) : Eigen::Vector3d::Zero();
            const Eigen::Vector3d nc_b = (cap_on && corner != kFace) ? cap_grad(s_in, o_b, kFace)
                                                                     : Eigen::Vector3d::Zero();
            auto fcap_o = [&](const Order& oo) {
                const double qt = s_in(oo.a) + (delta - 1.0) * s_in(oo.b) - delta * s_in(oo.c);
                const double pm = mean(s_in);
                return qt * qt / (alpha * alpha) + pm * pm - pp * pp;
            };
            const Eigen::Vector3d nc_avg =
                (cap_on && corner != kFace) ? cap_grad(s_in, o, corner) : nc;
            const double fc_in = cap_on ? fcap_o(o) : 0.0;
            const double fc_b = (cap_on && corner != kFace) ? fcap_o(o_b) : 0.0;
            bool ac_b = ac && corner != kFace && !pinned;
            if (ton[0] && ton[1] && ton[2]) { as = false; ac = ac_b = false; }
            // Koiter: a surface with a negative multiplier is dropped and the rest re-solved.
            for (int pass = 0; pass < 12; ++pass) {
                std::array<Eigen::Vector3d, 7> ns, ms;
                std::array<double, 7> f0s{};
                std::array<int, 7> kind{};   // 0 = shear plane, 1 = cap, 2 = tension cut-off
                int na = 0;
                if (as)
                    for (int i = 0; i < nplanes; ++i)
                        if (on[i]) { ns[na] = planes[i].n; ms[na] = planes[i].m; f0s[na] = fp[i]; kind[na] = 0; ++na; }
                // Where BOTH shear planes of the corner are active they already hold the stress on
                // the edge (equal values of the two planes mean equal stresses), and the cap flows
                // there with the one gradient the two faces share by symmetry. Two cap faces on
                // top of them would be four surfaces in a three-dimensional stress space, and
                // that system is singular -- measured: a multiplier of exactly zero and the edge
                // split by 2e-4 of the mean stress in undrained extension.
                const bool shear_edge = as && nplanes == 2 && on[0] && on[1];
                const bool cap_avg = shear_edge && (ac || ac_b);
                if (cap_avg) {
                    ns[na] = nc_avg; ms[na] = nc_avg; f0s[na] = std::max(fc_in, fc_b); kind[na] = 1; ++na;
                } else {
                    if (ac) { ns[na] = nc; ms[na] = nc; f0s[na] = fc_in; kind[na] = 1; ++na; }
                    if (ac_b) { ns[na] = nc_b; ms[na] = nc_b; f0s[na] = fc_b; kind[na] = 1; ++na; }
                }
                for (int i = 0; i < 3; ++i)
                    if (ton[i]) {
                        ns[na] = -Eigen::Vector3d::Unit(i); ms[na] = ns[na];
                        f0s[na] = -s_in(i) - ten; kind[na] = 2; ++na;
                    }
                if (na == 0) { sol.dsig = De * de_in; sol.tangent = De; sol.as = sol.ac = sol.at = false; return false; }
                using SmallMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 7, 7>;
                using SmallVec = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 7, 1>;
                SmallMat A(na, na);
                SmallVec b(na);
                for (int i = 0; i < na; ++i) {
                    // The linearised consistency f + m.dsigma - h dlambda = 0 of each surface,
                    // from the value it STARTS at. Without f the multiplier was the one a state ON
                    // the surface would take: a plane the state sat inside (while the cap was
                    // active) hardened past the stress, went inactive at the next evaluation,
                    // and the substeps alternated between the two -- a response that jumped by
                    // 1e-3 when a lateral strain of 1e-8 flipped the first activation, and a
                    // Newton iteration on an oedometer that stalled at a relative 1e-4.
                    b(i) = f0s[i] + ms[i].dot(De * de_in);
                    for (int j = 0; j < na; ++j) {
                        A(i, j) = ms[i].dot(De * ns[j]);
                        // The planes share gamma_p: every shear multiplier hardens every plane.
                        if (kind[i] == 0 && kind[j] == 0) A(i, j) += H;
                        if (kind[i] == 1 && kind[j] == 1) A(i, j) += 4.0 * pp * pmean * Hcap;
                    }
                }
                const Eigen::FullPivLU<SmallMat> lu(A);
                // More surfaces than the stress space has room for, through one point -- a cut-off
                // plane through the shear cone's own apex (c = 0, sigma_t = 0): the cut-off is the
                // tighter set there, so it is kept and the cone and cap are set aside.
                if (lu.rank() < na && (ton[0] || ton[1] || ton[2]) && (as || ac || ac_b)) {
                    as = false; ac = ac_b = false;
                    continue;
                }
                const SmallVec dl = lu.solve(b);
                bool dropped = false;
                int idx = 0;
                if (as)
                    for (int i = 0; i < nplanes; ++i)
                        if (on[i]) { if (dl(idx) < 0.0) { on[i] = 0; dropped = true; } ++idx; }
                if (cap_avg) {
                    if (dl(idx) < 0.0) { ac = ac_b = false; dropped = true; }
                    ++idx;
                } else {
                    if (ac) { if (dl(idx) < 0.0) { ac = false; dropped = true; } ++idx; }
                    if (ac_b) { if (dl(idx) < 0.0) { ac_b = false; dropped = true; } ++idx; }
                }
                for (int i = 0; i < 3; ++i)
                    if (ton[i]) { if (dl(idx) < 0.0) { ton[i] = 0; dropped = true; } ++idx; }
                if (dropped) {
                    bool any = false;
                    for (int i = 0; i < nplanes; ++i) any = any || on[i];
                    if (!any) as = false;
                    continue;
                }
                Eigen::Vector3d ep = Eigen::Vector3d::Zero();
                double dgp = 0.0, dev = 0.0;
                Eigen::Matrix<double, 3, Eigen::Dynamic, 0, 3, 7> N(3, na), M(3, na);
                for (int i = 0; i < na; ++i) {
                    ep += dl(i) * ns[i];
                    N.col(i) = ns[i];
                    M.col(i) = ms[i];
                    if (kind[i] == 0) dgp += H * dl(i);
                    else if (kind[i] == 1) dev += hi * dl(i) * 2.0 * pmean;   // the cut-off hardens nothing
                }
                sol.dsig = De * (de_in - ep);
                sol.dgp = dgp;
                sol.dev = dev;
                sol.tangent = De - De * N * lu.solve(M.transpose() * De);
                sol.as = as;
                sol.ac = ac || ac_b;
                sol.at = ton[0] || ton[1] || ton[2];
                return true;
            }
            sol.dsig = De * de_in; sol.tangent = De; sol.as = sol.ac = sol.at = false;
            return false;
        };

        // The corner is decided by the RATE, as a Mohr-Coulomb return decides its edges: flow on
        // the face first, and if that flow would carry the intermediate stress past the minor
        // (or the major past the intermediate), the state is on the corner and both faces flow.
        //
        // A step can break BOTH orderings at once: from a state whose two larger values are
        // equal, the face flow lowers one of them and raises the minor past the other. Taking the
        // compression corner first there -- as this did until the fix -- paired the minor with one
        // of the two equal stresses, split them, and let round-off decide which one; the answer
        // then jumped between two values as the strain increment moved by 1e-7, and a Newton
        // iteration on a weightless oedometer stalled at a relative residual of 1e-6. Both corners
        // are solved instead and the one whose own result keeps the ordering it assumed is taken;
        // if both or neither do, the pair that was closer to equal is the corner.
        Sol sol;
        int corner = kFace;
        solve(kFace, sol);
        if (allow & 4) {
            const Eigen::Vector3d s1 = s_in + sol.dsig;
            const bool past_c = s1(o.c) > s1(o.b), past_b = s1(o.b) > s1(o.a);
            if (past_c && !past_b) {
                corner = kCompression;
                solve(corner, sol);
            } else if (past_b && !past_c) {
                corner = kExtension;
                solve(corner, sol);
            } else if (past_c && past_b) {
                Sol sc, se;
                solve(kCompression, sc);
                solve(kExtension, se);
                const Eigen::Vector3d tc = s_in + sc.dsig, te = s_in + se.dsig;
                const bool keeps_c = tc(o.a) >= std::max(tc(o.b), tc(o.c));
                const bool keeps_e = te(o.c) <= std::min(te(o.a), te(o.b));
                if (keeps_c != keeps_e)
                    corner = keeps_c ? kCompression : kExtension;
                else
                    corner = (s_in(o.a) - s_in(o.b) < s_in(o.b) - s_in(o.c)) ? kExtension
                                                                             : kCompression;
                sol = corner == kCompression ? sc : se;
            }
        }
        r.dsig = sol.dsig;
        r.dgp = sol.dgp;
        r.dev = sol.dev;
        r.tangent = sol.tangent;
        r.as = sol.as;
        r.ac = sol.ac;
        r.at = sol.at;
        if (sol.face) corner = kFace;
        r.corner = corner;
        if (corner == kCompression) { r.tie_i = o.b; r.tie_j = o.c; }
        if (corner == kExtension) { r.tie_i = o.a; r.tie_j = o.b; }
        r.plastic = sol.as || sol.ac || sol.at;
        return r;
    };

    // WHERE THE ELASTIC TRIAL OF A STEP MEETS THE SURFACES: the fraction alpha of the strain
    // increment that is elastic, for a step that STARTS INSIDE both surfaces and whose elastic trial
    // ends outside one of them (Sloan, Abbo & Sheng 2001, Sec. 3.1; Illinois regula falsi, of the
    // same family as their Pegasus scheme, on the larger of the two surfaces each in its own units
    // -- continuous in the stress, zero where the path first meets either one). -1 when the step
    // does not cross.
    // `met`, when given, is told which surface the trial meets first (1 the shear planes, 2 the cap,
    // 3 the tension cut-off).
    auto elastic_fraction = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                                const Eigen::Vector3d& de_in, int* met = nullptr) -> double {
        const Stiff k = stiff_at(s_in);
        const double pp = pc_of(ev_in);
        auto f_at = [&](const Eigen::Vector3d& s) {
            const double fs = fshear_at(s, gp_in);
            const double fc = cap_on ? fcap(s, pp) / (1.0 + pp * pp) : -1.0;
            const double ft = ten_on ? (-s.minCoeff() - ten) / pr : -1.0;
            return std::max(std::max(fs, fc), ft);
        };
        const Eigen::Vector3d dse = k.De * de_in;
        const double f0 = f_at(s_in), f1 = f_at(s_in + dse);
        if (!(f0 < 0.0 && f1 > 0.0)) return -1.0;
        double a0 = 0.0, a1 = 1.0, g0 = f0, g1 = f1, alpha = 0.0;
        for (int it = 0; it < 60; ++it) {
            alpha = a1 - g1 * (a1 - a0) / (g1 - g0);
            const double ga = f_at(s_in + alpha * dse);
            if (std::fabs(ga) < 1e-14 || a1 - a0 < 1e-14) break;
            if ((ga > 0.0) == (g1 > 0.0)) { a1 = alpha; g1 = ga; g0 *= 0.5; }
            else { a0 = alpha; g0 = ga; g1 *= 0.5; }
        }
        if (met) {
            const Eigen::Vector3d sa = s_in + alpha * dse;
            const double fc = cap_on ? fcap(sa, pp) / (1.0 + pp * pp) : -1.0;
            const double ft = ten_on ? (-sa.minCoeff() - ten) / pr : -1.0;
            const double fs = fshear_at(sa, gp_in);
            *met = ft >= fs && ft >= fc ? 3 : (fs >= fc ? 1 : 2);
        }
        return alpha;
    };
    // One evaluation of a step that may cross: the part before the surface is elastic, only the
    // rest flows. Without the split the multiplier is computed from the whole strain increment as
    // if the start were on the surface, so the response jumps -- from elastic to a finite plastic
    // step -- where the trial crosses the surface; a Newton iteration whose points sit there (every
    // point of a column that the cap correction left just inside its cap) then flips points
    // between the two answers and stalls. The walks below split their SUBSTEPS there as well (see
    // substep()), so that this inner split only meets what is left of a crossing.
    auto increment = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                         const Eigen::Vector3d& de_in, int allow = 15) -> Inc {
        const double alpha = elastic_fraction(s_in, gp_in, ev_in, de_in);
        if (alpha < 0.0) return increment_core(s_in, gp_in, ev_in, de_in, allow);
        const Eigen::Vector3d dse = stiff_at(s_in).De * de_in;
        Inc r = increment_core(s_in + alpha * dse, gp_in, ev_in, (1.0 - alpha) * de_in, allow);
        r.dsig += alpha * dse;
        return r;
    };
    // THE SECOND EVALUATION OF A PAIR FLOWS ON WHAT THE FIRST DID. In Sloan's scheme it is the
    // elastoplastic matrix AT the Euler end point times the substep's strain (his algorithm step 8)
    // -- a rate, not a second finite step with a crossing of its own. Taken as a finite step from
    // the end point, it looked a whole substep further ahead and switched on a surface, or the
    // corner, that the substep never reached: a substep ending just inside the cap was elastic in
    // its first evaluation and plastic in its second (measured: R 3.7e-5 at 400 substeps, rising
    // as they shrank). Every crossing inside a substep is located before its pair (substep()), so
    // the second evaluation is held to the surfaces, and the corner, the first one found active;
    // it may still drop one (Koiter), as unloading must.
    auto activated = [&](const Inc& k) {
        return (k.as ? 1 : 0) | (k.ac ? 2 : 0) | (k.plastic && k.corner != kFace ? 4 : 0) | (k.at ? 8 : 0) |
               k.floor_side;
    };

    // --- DRIFT CORRECTION, applied to an ACCEPTED substep ------------------------------------
    // A stress violating a surface is pulled back along the CONSISTENT elastoplastic direction
    // De.n (Potts & Gens 1985): dsigma = -(f/(a.De.n + h)) De.n, a = df/dsigma. At a corner n is
    // the sum of the two faces' flows and h the gamma_p that sum hardens by, so a symmetric state
    // is corrected symmetrically. The cap is associated and projected on its own gradient. The
    // surfaces corrected are the ones the substep was yielding on.
    //
    // Which corner the corrected state is on is decided HERE, from the pair of directions the
    // substep found equal and the ordering of the state being corrected. Handing over the
    // substep's corner type instead was wrong whenever the ordering had changed within the step:
    // a step that started with the two lateral stresses equal and ABOVE the axial one (an
    // extension corner) and ended with them below it (a compression corner) was corrected as an
    // extension corner of the new ordering -- which pairs the axial stress with one lateral, and
    // split two stresses that had been equal.
    auto corner_of = [](const Order& o, int ti, int tj) {
        auto is = [&](int x, int y) { return (x == ti && y == tj) || (x == tj && y == ti); };
        if (ti < 0) return (int)kFace;
        if (is(o.b, o.c)) return (int)kCompression;
        if (is(o.a, o.b)) return (int)kExtension;
        return (int)kFace;
    };
    auto correct_drift = [&](Eigen::Vector3d& s, double& gp_io, double& ev_io,
                             bool as_act, bool ac_act, bool at_fail, int tie_i, int tie_j,
                             bool as_both = false, bool ac_both = false) -> bool {
        // BOTH WAYS for a surface the substep flowed on throughout (Sloan, Abbo & Sheng 2001,
        // Sec. 3.3: the correction is made when |f| > FTOL). The pair's mean leaves the stress off the
        // surface by O(R) either way; left INSIDE, the next substep read the surface as one not yet
        // met and located it as an event every substep of a path that stays on it (a normally
        // consolidated load: 2.3 million such events in a staged excavation). Only where both
        // evaluations flowed on it -- a substep that turned to unloading inside is left inside.
        for (int it = 0; it < 5; ++it) {
            const Stiff k = stiff_at(s);   // the correction moves sigma3, so it moves these too
            const Order o = order_of(s);
            const int corner = corner_of(o, tie_i, tie_j);
            const double ppc = pc_of(ev_io);
            const double qd = s(o.a) - s(o.c);
            const bool fail_here = at_fail || beyond(qd, k);
            const double fs = fail_here ? (qd - k.qf) : (fbar(qd, k) - gp_io);
            const double fcp = cap_on ? fcap(s, ppc) : -1.0;
            bool corr = false;
            const double tol_s = 1e-9 * (1.0 + std::fabs(gp_io) + std::fabs(k.qf));
            if (as_act && (fs > tol_s || (as_both && fs < -tol_s))) {
                const double spm = spm_of(qd, k);
                const Plane main = plane(s, o.a, o.c, spm, fail_here, k);
                Eigen::Vector3d n = main.n;
                double h = fail_here ? 0.0 : p.hardening_scale;
                if (corner == kCompression) { n += plane(s, o.a, o.b, spm, fail_here, k).n; h *= 2.0; }
                if (corner == kExtension) { n += plane(s, o.b, o.c, spm, fail_here, k).n; h *= 2.0; }
                const Eigen::Vector3d Den = k.De * n;
                const double dlam_d = fs / (main.m.dot(Den) + h);
                s -= dlam_d * Den;
                gp_io += h * dlam_d;
                corr = true;
            }
            // The cap is corrected the same consistent way, with its hardening: along De.g, and
            // pp moving with the plastic volume that step takes. Only a state OUTSIDE is
            // corrected. The bare projection this replaces, s -= f g/|g|^2 on |f|, also pulled a
            // state that had ended inside the cap back OUT onto it -- at constant pp, along a
            // gradient whose volumetric part is always compressive -- and that raised the mean
            // effective stress of an undrained psi = 0 path, where it cannot rise.
            const double tol_c = 1e-8 * (1.0 + ppc * ppc);
            if (ac_act && (fcp > tol_c || (ac_both && fcp < -tol_c))) {
                const Eigen::Vector3d g = cap_grad(s, o, corner);
                const Eigen::Vector3d Deg = k.De * g;
                const double pm = mean(s);
                const double Hc_ = p.hardening_scale * p.cap_hardening_modulus(ppc);
                const double dlam_c = fcp / (g.dot(Deg) + 4.0 * ppc * pm * Hc_);
                s -= dlam_c * Deg;
                ev_io += p.hardening_scale * dlam_c * 2.0 * pm;
                corr = true;
            }
            if (!corr) break;
        }
        // THE TENSION CUT-OFF, by its own associated return (the Rankine planes, with the elastic
        // constants at this state): whatever the integration left past sigma_t is brought back to
        // it, and the Mohr-Coulomb bound below then reads the minor stress the cut-off allows.
        if (ten_on && -s.minCoeff() - ten > 1e-12 * (1.0 + std::fabs(ten) + s.cwiseAbs().maxCoeff())) {
            const double E = stiff_at(s).Eur * p.elastic_scale;
            const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), mu = E / (2.0 * (1.0 + nu));
            int perm[3] = {0, 1, 2};
            std::sort(perm, perm + 3, [&](int a, int b) { return s(a) < s(b); });   // most tensile first
            double r[3] = {-s(perm[0]), -s(perm[1]), -s(perm[2])};                 // tension positive, descending
            apply_rankine_cap(r, ten, lam, mu, nullptr);
            for (int i = 0; i < 3; ++i) s(perm[i]) = -r[i];
        }
        // MC failure bound (shear): q <= q_f, the last line. With the failure gradient above the
        // drift loop lands on the surface itself; this only removes what is left of the
        // integration error. At the extension corner both major stresses are lowered, so the
        // corner stays one.
        {
            const Stiff k = stiff_at(s);
            const Order o = order_of(s);
            const int corner = corner_of(o, tie_i, tie_j);
            const double over = (s(o.a) - s(o.c)) - k.qf;
            if (k.qf <= 0.0) {
                // Past the apex no deviator is admissible: the stress IS the apex, -c cot(phi)
                // (the tension cut-off applied after this return then takes what it allows).
                if (over > 0.0) {
                    s.setConstant(law_s > 1e-12 ? -law_cc / law_s : 0.0);
                    return true;
                }
            } else if (over > 1e-12 * (1.0 + k.qf)) {
                s(o.a) -= over;
                if (corner == kExtension) s(o.b) -= over;
            }
            // A state ON the failure line is on the hardening surface too: the hyperbola ends at
            // q_f, so gamma_p there is at least fbar(q_f). A state brought to the line by this
            // bound (or by the cap, or by a path that reached it before it hardened -- a column
            // loaded from zero, whose first elastic trial lies past the line) carried the gamma_p
            // of before, far inside what the stress had reached: an admissible stress on an
            // inadmissible surface. The next increment then chose between the plateau and a
            // hardening return that pulled the whole stress back to the stale surface, and a
            // perturbation of 1e-6 of the Newton step flipped it -- measured: 1.2 kPa jumps in a
            // column at 1.25 kPa, holding the equilibrium iteration at 2e-3.
            const double qn = s.maxCoeff() - s.minCoeff();
            const Stiff kn = stiff_at(s);
            if (kn.qf > 0.0 && qn >= kn.qf - 1e-9 * (1.0 + kn.qf))
                gp_io = std::max(gp_io, fbar(kn.qf, kn));
        }
        return false;
    };

    bool any_plastic = false;
    // ON THE TIP AND NEAR IT. A stress-free point that has never been loaded -- freshly placed
    // fill -- stands at the tip of its cap: p = p_p = 0, where the cap's gradient is zero. An
    // explicit step from there reads no cap at all and is elastic; the same step evaluated at its
    // own end is plastic. The pair disagree by the same proportion at every step size, so the
    // measured subdivision ran to its ceiling on every call and reported an integration that
    // had not met its tolerance. Measured on the first time step of an embankment placed on
    // 15-noded elements: every stress point of the fill at the ceiling of 200 substeps, the
    // balance held at 9e-3 of its scale for sixty iterations of seven seconds each, the step
    // halved, and the same again. Three cases are told apart below.
    // NEAR THE TIP: a point that has left the tip by half a kPa and is asked for an increment of
    // ten is in the same position: its surfaces grow many times over within the increment. Both
    // are walked in graded substeps (below) whenever the elastic trial is more than four times the
    // size of the state (the larger of its stress and its cap pressure).
    Eigen::Vector3d s0 = sigma_n, dstrain_w = dstrain;
    double gp0 = gamma_p_n, ev0 = ev_n;

    // The environment override is for a RUN's tolerance. The calibration integrates at its own
    // (kHsCalibrationTol), which the override must not reach: alpha and beta are properties of the
    // material, and a sweep of the run tolerance that silently re-calibrated the cap measured two
    // things at once.
    const double env_tol = stol == kHsCalibrationTol ? 0.0 : hs_env_substep_tol();
    const double tol = env_tol > 0.0 ? env_tol
                                     : (stol > 0.0 ? stol : hs_default_substep_tol());
    // Ceiling on the subdivision. A GUARD, the place of Sloan's minimum step dT_min, not a policy:
    // the tolerance decides the count, and a walk the ceiling stopped short of it is reported.
    const int kMaxSubsteps = hs_max_substeps();
    // The relative local error of one substep (Sloan, Abbo & Sheng 2001, Eqn 24, the stress part;
    // see WHAT THE TOLERANCE MEANS below). The denominator has a FLOOR: at a point that carries next
    // to nothing -- fresh fill, the ground surface -- a relative tolerance is an absolute one of
    // nothing, and resolving stresses of a thousandth of a kPa asked for the ceiling. p_ref / 20
    // ties it to the stiffness law's own stress level (measured: a fill on a soft base 141 s ->
    // 97 s at p_ref / 200 -> p_ref / 20, same answer).
    auto rel_error = [&](const Eigen::Vector3d& d1, const Eigen::Vector3d& d2, const Eigen::Vector3d& s_hat) {
        return 0.5 * (d2 - d1).norm() / std::max(s_hat.norm(), pr / 20.0);
    };
    // A walk is refined until its largest R is under 0.81 STOL -- Sloan's 0.9 on the step (his
    // Eqn 27), squared onto R -- by at most this many walks.
    constexpr double kAim = 0.81;
    constexpr int kMaxWalks = 8;
    // ONE SUBSTEP OF A WALK, split where its elastic trial meets a surface -- Sloan's order
    // (Sloan, Abbo & Sheng 2001, algorithm steps 3-5): the intersection first, the elastic part by
    // the elastic law alone, and the modified-Euler pair only over the plastic rest, from the
    // surface. With the split inside ONE evaluation instead, the pair straddled the kink: the
    // first evaluation was elastic up to the surface and the second, from its end on the surface,
    // plastic throughout, so the two differed by (D_ep - D_e) times the elastic strain to the
    // surface -- an amount that does not shrink with the substep -- and their mean left the stress
    // inside by half that distance, which the next substep halved again. Measured on a staged
    // excavation in sand: the largest R of a walk stayed at 2.05e-5 from 50 to 400 substeps,
    // halving substep by substep, and walks were refined to thousands of substeps for it.
    // The elastic part uses the stress-dependent E_ur at both ends (the elastic law is nonlinear
    // here) and its own R is measured with the rest. `whole`, when given, is the pair already
    // evaluated over this substep from this state, used if no split is needed.
    auto substep = [&](Eigen::Vector3d& sg, double& g, double& e, const Eigen::Vector3d& de,
                       Eigen::Matrix3d& tang, bool& plastic, double& r_max, bool elastic_tangent_at_end,
                       const Inc* whole1, const Inc* whole2) {
        poll_interrupt();   // a run's Cancel reaches inside one point's walk (interrupt.hpp)
        Eigen::Vector3d de_p = de;
        int met = 0;   // the surface the elastic part ends on, when there is one
        const double alpha = elastic_fraction(sg, g, e, de, &met);
        if (alpha >= 0.0) {
            const Eigen::Vector3d de_e = alpha * de;
            const Eigen::Vector3d d1 = stiff_at(sg).De * de_e;
            const Eigen::Vector3d d2 = stiff_at(sg + d1).De * de_e;
            sg += 0.5 * (d1 + d2);
            r_max = std::max(r_max, rel_error(d1, d2, sg));
            de_p = (1.0 - alpha) * de;
            whole1 = whole2 = nullptr;
        }
        // THE EVENTS OF THE PLASTIC PART. On the surfaces it has met, the flow can still change what
        // it flows on within the substep: meet the OTHER surface -- the cap and the shear planes
        // meet near the normally consolidated state, where much of a loaded ground lies -- or close
        // two principal stresses onto an EDGE of the hexagon, from where it flows on both faces (a
        // K0 state loaded in compression closes on its compression corner). Taken inside one
        // evaluation, each of these is the kink again: the first evaluation of the pair on the old
        // set, the second from its end on the new, the two apart by an amount that does not shrink
        // with the substep, and their mean leaving the stress short of the new set by half of it.
        // Measured on the staged excavation: walks refined to thousands of substeps at R constant
        // from 3 to 412 of them -- inside both surfaces by 1e-6 .. 1e-3 of their scale, or two
        // stresses 0.018 kPa from equal.
        //
        // So the plastic part is walked from event to event, as the elastic part is ended at the
        // surface: from the set active where the piece starts (the surfaces it is on, the edge if
        // it is on one), the flow HELD to that set is followed over the rest of the substep; where
        // it would meet another surface, or make two stresses equal, the fraction beta of the strain
        // at which it does is found (Illinois regula falsi, as for the surface), the earliest one
        // is taken, the pair is taken over that piece held to the set, and the rest starts again
        // from the new state. This is Sloan's intersection (Sloan, Abbo & Sheng 2001, Sec. 3.1)
        // carried to the other events of a model with two surfaces and corners; his models have one
        // smooth surface.
        // The held step over the rest of the substep, when the loop ends on it without an event: it
        // IS the first evaluation of the final pair, and the better one -- the unrestricted step would
        // switch on, through its linearised consistency, a surface the flow does not reach.
        Inc k1_held;
        bool k1_ready = false;
        for (int piece = 0; piece < 4 && !whole1; ++piece) {
            const double pp0 = pc_of(e);
            const double fs0 = fshear_at(sg, g);
            const double fc0 = cap_on ? fcap(sg, pp0) / (1.0 + pp0 * pp0) : -1.0;
            int held = 0;
            // ON a surface, or an edge, is a start close enough that the kink it can open is under a
            // tenth of the tolerance: a start at stress distance d from where a set changes opens a
            // kink of at most R = d / (2 |sigma|) (the elastic strain to it, through D_e - D_ep), so
            // d <= 0.2 STOL |sigma| is on it. The distance is read along the substep's elastic trial,
            // by the secant of the surface's value. Read stricter (to the drift correction's own
            // 1e-9), the walk on both surfaces of an oedometer -- which the pair leaves inside one or
            // the other by O(R) every substep, and the correction pulls back only from outside --
            // opened an event piece in half of 1.4 million substeps, at 2 million root evaluations.
            const double scale = std::max(sg.norm(), pr / 20.0);
            const double near = 0.2 * tol * scale;
            const Eigen::Vector3d ds_e = stiff_at(sg).De * de_p;
            auto distance = [&](double v0, double v1) {
                if (!(v0 < 0.0)) return 0.0;
                if (!(v1 > v0)) return 1e300;   // moving away from it
                return -v0 / (v1 - v0) * ds_e.norm();
            };
            const double fs1 = fshear_at(sg + ds_e, g);
            const double fc1 = cap_on ? fcap(sg + ds_e, pp0) / (1.0 + pp0 * pp0) : -1.0;
            const double ft0 = ten_on ? -sg.minCoeff() - ten : -1e300;
            const double ft1 = ten_on ? -(sg + ds_e).minCoeff() - ten : -1e300;
            if (distance(fs0, fs1) <= near || met == 1) held |= 1;
            if ((cap_on && distance(fc0, fc1) <= near) || met == 2) held |= 2;
            if ((ten_on && distance(ft0, ft1) <= near) || met == 3) held |= 8;
            // The cut-off is THREE planes, one per principal direction, and the ones the stress is
            // on are its smallest values: nt_on of them. The next one to meet it is an event of its
            // own (below) -- held to "the cut-off", the step would switch it on unannounced.
            int nt_on = 0;
            if (held & 8) {
                for (int i = 0; i < 3; ++i)
                    if (distance(-sg(i) - ten, -(sg + ds_e)(i) - ten) <= near) ++nt_on;
                nt_on = std::max(nt_on, 1);
            }
            met = 0;
            if (!held) break;   // inside both: nothing to flow on, and no crossing (see above)
            const Order oc = order_of(sg);
            auto gaps = [&](const Eigen::Vector3d& x) {
                return Eigen::Vector2d(x(oc.b) - x(oc.c), x(oc.a) - x(oc.b));
            };
            const Eigen::Vector2d gap0 = gaps(sg);
            if (gap0(0) <= near || gap0(1) <= near) held |= 4;
            // ON THE STIFFNESS FLOOR (within the same distance): the gradient is read on the side the
            // flow is going to, in both evaluations, as a start on a surface flows on it.
            const double s3_0 = sg.minCoeff();
            if ((held & 1) && std::fabs(s3_0 - plim) <= near) {
                const Inc probe = increment_core(sg, g, e, de_p, held);
                held |= (sg + probe.dsig).minCoeff() > s3_0 ? 16 : 32;
            }
            const Inc whole_held = increment_core(sg, g, e, de_p, held);
            if (!whole_held.plastic) { k1_held = whole_held; k1_ready = true; break; }
            // The (k+1)-th smallest principal stress: the next one to meet the cut-off when k are on it.
            auto nth = [](const Eigen::Vector3d& x, int k) {
                double v[3] = {x(0), x(1), x(2)};
                std::sort(v, v + 3);
                return v[k];
            };
            // THE STIFFNESS FLOOR is where the shear planes' gradient jumps: below p_limit E_i, E_ur
            // no longer move with sigma3 and the dfbar/dsigma3 term of the flow (plane()) is zero, at
            // it that term is finite. The stress the laws read is continuous there, the direction the
            // flow takes is not -- the same first-order kink as a surface met inside a substep
            // (measured: a walk from rest refined to the ceiling of 10 000 substeps by the one
            // substep that crossed sigma3 = p_limit, R 9.2e-5 -> 3.2e-5 for substeps 3.3 times
            // shorter). Which way the minor stress crosses it is read from where it starts.
            const double floor_sign = s3_0 < plim ? 1.0 : -1.0;
            // The event functions: negative before the event, positive after it.
            auto event = [&](int kind, const Inc& k) {
                const Eigen::Vector3d x = sg + k.dsig;
                if (kind == 0) return fshear_at(x, g + k.dgp);
                if (kind == 1) {
                    const double ppk = pc_of(e + k.dev);
                    return fcap(x, ppk) / (1.0 + ppk * ppk);
                }
                if (kind == 4) return (x(oc.a) - x(oc.c)) - stiff_at(x).qf;   // the failure plateau
                if (kind == 5) return -nth(x, nt_on) - ten;                  // the next cut-off plane
                if (kind == 6) return floor_sign * (x.minCoeff() - plim);    // the stiffness floor
                return -gaps(x)(kind - 2);   // kind 2: the compression edge, 3: the extension edge
            };
            double beta = 2.0;
            // THE PLATEAU is an event too: hardening on the shear planes (H > 0) ends where q reaches
            // q_f and the flow is perfectly plastic from there (H = 0, the failure gradient).
            const double qf0 = stiff_at(sg).qf;
            const double plateau0 = (sg(oc.a) - sg(oc.c)) - qf0;
            const bool hardening = (held & 1) && plateau0 < -1e-9 * (1.0 + std::fabs(qf0));
            for (int kind = 0; kind < 7; ++kind) {
                if (kind == 0 && (held & 1)) continue;
                if (kind == 1 && (!cap_on || (held & 2))) continue;
                if ((kind == 2 || kind == 3) && (held & 4)) continue;
                if (kind == 4 && !hardening) continue;
                if (kind == 5 && (!ten_on || nt_on >= 3)) continue;
                if (kind == 6 && (!hardening || (held & 48))) continue;
                const double q_start = kind == 0 ? fs0 : kind == 1 ? fc0 : kind == 4 ? plateau0
                                     : kind == 5 ? -nth(sg, nt_on) - ten
                                     : kind == 6 ? floor_sign * (s3_0 - plim) : -gap0(kind - 2);
                const double q_end = event(kind, whole_held);
                if (!(q_start < 0.0 && q_end > 0.0)) continue;
                double b0 = 0.0, b1 = 1.0, q0 = q_start, q1 = q_end, b = 0.0;
                for (int it = 0; it < 4; ++it) {
                    b = b1 - q1 * (b1 - b0) / (q1 - q0);
                    const double qb = event(kind, increment_core(sg, g, e, b * de_p, held));
                    // To 1e-3 of the event function's range, in a few evaluations: what an error in
                    // beta leaves is a kink of that fraction of the jump, at the start of the next
                    // piece -- whose own R is measured and refined like any other, so the location
                    // need not be exact, only close. NOT CLOSER: the event function is read through
                    // the held step, whose active set can change inside the bracket, and iterated to
                    // round-off (tried 2026-10-29) the location converged onto that switch, where a
                    // perturbation of the strain moves it discontinuously -- a step of 1e-6 kPa in
                    // the stress over 1e-11 of strain, read by a central difference as a slope of
                    // 1e4 (the tension cut-off tangent check, 6e-9 -> 0.43).
                    if (std::fabs(qb) <= 1e-3 * std::max(-q_start, q_end) || b1 - b0 < 1e-6) break;
                    if ((qb > 0.0) == (q1 > 0.0)) { b1 = b; q1 = qb; q0 *= 0.5; }
                    else { b0 = b; q0 = qb; q1 *= 0.5; }
                }
                beta = std::min(beta, b);
            }
            if (beta > 1.0) { k1_held = whole_held; k1_ready = true; break; }   // no event in the substep
            // The piece up to the event is as long as the event is far, whatever the substep: Sloan's
            // own local rule refines it (his Sec. 3.4), into m equal parts of a real count, m grown by
            // (R / 0.81 STOL)^(1/2) while its largest R is above that -- continuous, as the walks are.
            const Eigen::Vector3d de_b = beta * de_p;
            const Eigen::Vector3d sg_b = sg;
            const double g_b = g, e_b = e;
            double m = 1.0, r_piece = 0.0;
            for (int tries = 0; tries < 4; ++tries) {
                sg = sg_b; g = g_b; e = e_b;
                r_piece = 0.0;
                const double dm = 1.0 / m;
                int parts = 0;
                for (double T = 0.0; T < 1.0 - 1e-12 && parts < (int)std::ceil(m - 1e-9) + 1; ++parts) {
                    const Eigen::Vector3d dep = std::min(dm, 1.0 - T) * de_b;
                    const Inc h_1 = increment_core(sg, g, e, dep, held);
                    const Inc h_2 = increment_core(sg + h_1.dsig, g + h_1.dgp, e + h_1.dev, dep, held & activated(h_1));
                    sg += 0.5 * (h_1.dsig + h_2.dsig);
                    g += 0.5 * (h_1.dgp + h_2.dgp);
                    e += 0.5 * (h_1.dev + h_2.dev);
                    r_piece = std::max(r_piece, rel_error(h_1.dsig, h_2.dsig, sg));
                    if (h_1.plastic || h_2.plastic) {
                        plastic = true;
                        tang = h_2.plastic ? h_2.tangent : h_1.tangent;
                        if (correct_drift(sg, g, e, h_1.as || h_2.as, h_1.ac || h_2.ac,
                                          h_1.at_fail || h_2.at_fail, h_2.plastic ? h_2.tie_i : h_1.tie_i,
                                          h_2.plastic ? h_2.tie_j : h_1.tie_j,
                                          h_1.as && h_2.as, h_1.ac && h_2.ac))
                            tang.setZero();
                    }
                    T += std::min(dm, 1.0 - T);
                }
                if (!(r_piece > kAim * tol) || m >= 256.0) break;
                m = std::min(m * std::sqrt(r_piece / (kAim * tol)), 256.0);
            }
            r_max = std::max(r_max, r_piece);
            de_p = (1.0 - beta) * de_p;
        }
        const Inc k1 = whole1 ? *whole1 : (k1_ready ? k1_held : increment(sg, g, e, de_p));
        const Inc k2 = whole2 ? *whole2 : increment_core(sg + k1.dsig, g + k1.dgp, e + k1.dev, de_p, activated(k1));
        sg += 0.5 * (k1.dsig + k2.dsig);
        g += 0.5 * (k1.dgp + k2.dgp);
        e += 0.5 * (k1.dev + k2.dev);
        // Measured before the drift correction, as Sloan's steps 10 and 13 are ordered.
        r_max = std::max(r_max, rel_error(k1.dsig, k2.dsig, sg));
        if (k1.plastic || k2.plastic) {
            plastic = true;
            tang = k2.plastic ? k2.tangent : k1.tangent;
            // A state returned to the apex does not depend on the strain at all: its derivative is
            // zero (the solver's vertex floor keeps such a point from making K_T singular).
            if (correct_drift(sg, g, e, k1.as || k2.as, k1.ac || k2.ac, k1.at_fail || k2.at_fail,
                              k2.plastic ? k2.tie_i : k1.tie_i, k2.plastic ? k2.tie_j : k1.tie_j,
                              k1.as && k2.as, k1.ac && k2.ac))
                tang.setZero();
        } else if (elastic_tangent_at_end) {
            tang = stiff_at(sg).De;
        }
    };
    // A POINT WITHOUT COHESION THAT ENDS AT THE APEX IS RETURNED THERE. With c = 0 and the cut-off
    // through the apex (sigma_t = 0), a point whose minor stress is in the stiffness floor has an
    // elastic law that no longer moves with the stress, and the apex is the answer exactly when the
    // strain that would be plastic there, x = De^-1 sigma + d eps, is a flow the apex has: the
    // cut-off planes give any extension, and a shear face compression in its major direction with
    // (1 + sin psi) / (1 - sin psi) as much extension in its minor -- so the apex holds the increment
    // when (1 + sin psi) sum(x+) <= (1 - sin psi) sum(x-). On the failure plateau (gamma_p at least
    // the value the stress fails at, which falls to zero with the stress) the way there hardens
    // nothing and compacts nothing, as at the apex itself (the rule below). Walked instead, the point
    // slid along the cut-off to the apex with the last substep's elastic part left over, its error
    // falling as 1/n: such points -- the ground heaved and torn beside a wall pushed into sand, the
    // floor of an excavation -- ran their walks to the ceiling of 10 000 substeps, one to five
    // seconds each, and the first phase of a wall pushed 50 mm into Hardening Soil sand took 33
    // times as long as the same wall in Mohr-Coulomb sand. It is tried where the walk starts and at
    // every substep of it, so a point that enters the floor on its way is returned from there.
    const bool apex_rule = law_cc <= 0.0 && ten_on && std::fabs(ten) <= 1e-12 * pr;
    const double sin_psi_apex = std::clamp(dil(dil.sin_phi), -0.99, 0.99);
    auto ends_at_apex = [&](const Eigen::Vector3d& s, double g, const Eigen::Vector3d& de) -> bool {
        if (!apex_rule) return false;
        const Order o = order_of(s);
        if (!(s(o.c) <= plim)) return false;
        const Stiff k = stiff_at(s);
        if (!(g >= fbar(std::max(k.qf, 0.0), k))) return false;
        const Eigen::Vector3d x = k.De.inverse() * s + de;
        double pos = 0.0, neg = 0.0;
        for (int i = 0; i < 3; ++i) (x(i) > 0.0 ? pos : neg) += std::fabs(x(i));
        return pos + neg > 0.0 && (1.0 + sin_psi_apex) * pos <= (1.0 - sin_psi_apex) * neg;
    };
    {
        const double tiny = 1e-9 * pr;
        const Eigen::Vector3d ds_tr = k_n.De * dstrain;
        const double size0 = std::max(pp_n, sigma_n.cwiseAbs().maxCoeff());
        const double step = ds_tr.cwiseAbs().maxCoeff();
        // On the tip: no stress and no cap pressure (whatever the shear hardening stands at).
        const bool on_tip = size0 <= tiny;
        // AT THE APEX OF A CONE WITHOUT COHESION, PULLED, IT STAYS THERE. With c = 0 a stress-free
        // point has no strength (q_f = 0) whatever its cap pressure, so an increment whose elastic
        // trial is not compressive leaves it at the apex: its stress does not depend on the strain,
        // it is on the failure plateau (where gamma_p does not grow, H = 0) and it compacts nothing
        // (so p_p does not either) -- the cap tip's own rule below, for a cap that is not at its tip.
        // Walked instead, each substep returned to the apex but the last one's elastic part stayed,
        // so the stress went as 1/n (6.63 / 0.663 / 0.0663 kPa at 10 / 100 / 1000 substeps), the error
        // measure never settled, and such points -- the floor of an excavation -- ran every walk to the
        // ceiling: 203 of the staged excavation's 11 600 costliest calls, and nearly all of its time.
        if (step > 0.0 && ((law_cc <= 0.0 && sigma_n.cwiseAbs().maxCoeff() <= tiny && mean(ds_tr) <= 0.0) ||
                           ends_at_apex(sigma_n, gamma_p_n, dstrain))) {
            HsIntegrated apex;
            apex.stress = sigma_n.cwiseAbs().maxCoeff() <= tiny ? sigma_n : Eigen::Vector3d::Zero();
            apex.gamma_p = gamma_p_n;
            apex.pp = pp_n;
            apex.tangent.setZero();
            apex.plastic = true;
            apex.nsub = 0;
            apex.saturated = 0;
            if (plan_out) plan_out->clear();
            return apex;
        }
        // IN THE STIFFNESS FLOOR, PULLED ACROSS THE CUT-OFF, THE INCREMENT IS RETURNED IMPLICITLY.
        // Below p_limit the elastic stiffness no longer moves with the stress, so the elastic trial is
        // exact, and the fully implicit (backward Euler) return solves the increment in one step: the
        // stress and the multipliers of the surfaces active at its end (the closest-point return,
        // Simo & Hughes 1998, Sec. 3.3, with several surfaces by Koiter's rule, Koiter 1960) -- the
        // three cut-off planes, the shear face of the major and minor stresses and, at the compression
        // corner, the face of the major and intermediate (both harden the one gamma_p; on the plateau
        // the Mohr-Coulomb line), and the cap. The set is found by Koiter's rule from the cut-off planes
        // the trial crosses: a negative multiplier leaves it, a surface the answer crosses joins it.
        //
        // Walked instead, such a point -- fill with a cohesion, pulled into tension at the surface --
        // slides along the cut-off with the shear surface and the cap meeting it, the two evaluations
        // of each pair read the planes a hair apart, and the error estimate falls only as 1/n: walks at
        // the ceiling of 10 000 substeps, a tenth to four tenths of a second per call, and an embankment
        // on soft soil whose first phase ran an hour where it had run three minutes (measured: 517 of
        // the 551 thread-seconds of its integrations in 13 298 such calls).
        //
        // ITS ERROR IS THAT OF ONE IMPLICIT STEP: first order in the increment, the order of the
        // implicit time stepping that hands the increment over, and it vanishes with it (the return is
        // consistent). Measured on 532 calls sampled from that phase, against the limit of the walk
        // (100 000 substeps): the stress within 1e-2 of the 5 kPa the tolerance is read against, 5e-5
        // at the median, the hardening variables within 2%; 15 times less integration time. Where the
        // return does not converge to a set that holds -- multipliers, crossings, ordering, the floor,
        // the plateau reached by hardening -- the walk integrates the increment as before.
        //
        // A TEST OF PATH INDEPENDENCE WAS TRIED AND TAKEN OUT: the answer only where one step, two
        // and four agreed to the walk's tolerance. It misjudged both ways -- an increment whose
        // cut-off plane is met and left in its first 0.3% agrees at every count up to sixteen and is
        // 9e-4 off, and a point pulled from rest, 7e-5 off, did not agree -- and it made the
        // integrator a jump between two schemes wherever it changed its mind.
        struct BeOut {
            Eigen::Vector3d stress;
            double gp = 0.0, ev = 0.0;
            Eigen::Matrix3d tangent;
            bool plastic = false;
        };
        // One backward-Euler return over de_in from (s_in, g_in, e_in), the state in the floor.
        auto be_return = [&](const Eigen::Vector3d& s_in, double g_in, double e_in, const Eigen::Vector3d& de_in,
                             BeOut& res) -> bool {
            const double pp_in = pc_of(e_in);
            const Eigen::Vector3d s_tr = s_in + k_n.De * de_in;
            {
                const Order ot = order_of(s_tr);
                int ia = ot.a, ib = ot.b, ic = ot.c;   // the ordering the shear faces and the cap read
                const Eigen::Matrix3d& De = k_n.De;
                const double hi = p.hardening_scale;
                // A cap at the tip is seeded with a vanishing size, as the walk from the tip seeds it.
                const double pp_seed = 1e-6 * (1.0 + s_tr.cwiseAbs().maxCoeff());
                const double ev0 = (cap_on && pp_in < pp_seed) ? p.cap_ev_from_pc(pp_seed) : e_in;
                const double dqf = 2.0 * law_s / (1.0 - law_s);
                const double dqa = 2.0 * law_s / ((1.0 - law_s) * p.Rf);
                auto qf_of = [&](double sc) { return 2.0 * (law_cc + sc * law_s) / (1.0 - law_s); };
                auto qa_of = [&](double sc) { return std::max(qf_of(sc), qf_min) / p.Rf; };
                const double Ei = k_n.Ei, Eur = k_n.Eur;
                auto fb = [&](double q, double qa) { return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur; };
                Eigen::Vector3d dqt;
                auto set_order = [&](int a, int b, int c) {
                    ia = a; ib = b; ic = c;
                    dqt.setZero();
                    dqt(ia) = 1.0; dqt(ib) = delta - 1.0; dqt(ic) = -delta;
                };
                set_order(ia, ib, ic);
                int reorders = 0;
                int corner_adds = 0;
                using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 9, 9>;
                using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 9, 1>;
                // surfaces: 0..2 the cut-off plane of each direction, 3 the shear face (a, c), 4 the cap,
                // 5 the shear face (a, b) -- the compression corner
                constexpr int NS = 6;
                std::array<char, NS> act{};
                bool plateau = false;
                Eigen::Vector3d sg = s_tr;
                std::array<double, NS> dl{};
                double gp_end = g_in, ev_end = ev0, cap_dfdl = 0.0;
                // one shear face (i the larger, j the smaller): its value, gradient and flow
                auto face = [&](const Eigen::Vector3d& x, int i, int j, double gsum, double& fv, Eigen::Vector3d& nv,
                                Eigen::Vector3d& mv, double& gpv) {
                    const double q = x(i) - x(j), sj = x(j);
                    const double qf = qf_of(sj), qa = qa_of(sj);
                    Stiff kk = k_n;
                    kk.s3raw = sj;
                    const double spm = spm_of(q, kk);
                    mv.setZero();
                    mv(i) = 0.5 * (1.0 - spm);
                    mv(j) = -0.5 * (1.0 + spm);
                    nv.setZero();
                    if (plateau || !(qf > 0.0)) {
                        fv = q - qf;
                        nv(i) = 1.0;
                        nv(j) = -1.0 - dqf;
                        gpv = qf > 0.0 ? std::max(g_in, fb(qf * (1.0 - 1e-12), qa)) : g_in;
                    } else {
                        gpv = g_in + hi * gsum;
                        const double qq = std::min(q, qf * (1.0 - 1e-9));
                        const double w = 1.0 - qq / qa;
                        fv = fb(qq, qa) - gpv;
                        const double fp = (2.0 / Ei) / (w * w) - 2.0 / Eur;
                        nv(i) = fp;
                        nv(j) = -fp - (2.0 / Ei) * qq * qq / (qa * qa * w * w) * dqa;
                    }
                };
                auto eval = [&](const Eigen::Vector3d& x, const std::array<double, NS>& l, std::array<double, NS>& f,
                                std::array<Eigen::Vector3d, NS>& n, std::array<Eigen::Vector3d, NS>& m, double& gpv,
                                double& evv, Eigen::Matrix3d& dm_cap) {
                    for (int t = 0; t < 3; ++t) {
                        f[t] = -x(t) - ten;
                        n[t] = -Eigen::Vector3d::Unit(t);
                        m[t] = n[t];
                    }
                    double g5 = g_in;
                    face(x, ia, ic, l[3] + l[5], f[3], n[3], m[3], gpv);
                    face(x, ia, ib, l[3] + l[5], f[5], n[5], m[5], g5);
                    // the cap, written out for any mean stress (a yield surface only where it is positive)
                    const double pm = x.sum() / 3.0;
                    // at the compression corner the cap's two faces averaged (as increment_core does)
                    Eigen::Vector3d dqc = dqt;
                    if (act[5]) { dqc.setZero(); dqc(ia) = 1.0; dqc(ib) = -0.5; dqc(ic) = -0.5; }
                    const double qt = x.dot(dqc);
                    evv = ev0 + hi * 2.0 * pm * l[4];
                    const double pcv = pc_of(evv);
                    if (cap_on) {
                        f[4] = qt * qt / (alpha * alpha) + pm * pm - pcv * pcv;
                        const double Hc = p.cap_hardening_modulus(pcv);
                        m[4] = 2.0 * qt / (alpha * alpha) * dqc + (2.0 * pm / 3.0) * Eigen::Vector3d::Ones();
                        n[4] = m[4] - 2.0 * pcv * Hc * hi * 2.0 * l[4] / 3.0 * Eigen::Vector3d::Ones();
                        cap_dfdl = -2.0 * pcv * Hc * hi * 2.0 * pm;
                        dm_cap = (2.0 / (alpha * alpha)) * dqc * dqc.transpose() + (2.0 / 9.0) * Eigen::Matrix3d::Ones();
                    } else {
                        f[4] = -1.0;
                        m[4] = n[4] = Eigen::Vector3d::Zero();
                        cap_dfdl = 0.0;
                        dm_cap.setZero();
                    }
                };
                // d f_r / d dl_c: both shear faces harden the one gamma_p; the cap hardens p_p
                auto dfdl = [&](int r, int c) {
                    if ((r == 3 || r == 5) && (c == 3 || c == 5)) return plateau ? 0.0 : -hi;
                    if (r == 4 && c == 4) return cap_dfdl;
                    return 0.0;
                };
                // crossed: the cap only where the mean stress is positive; the face (a, b) only where the
                // intermediate stress is not held by its own cut-off plane (there the edge is the face)
                auto crossed = [&](int t, const Eigen::Vector3d& x, double fv) {
                    if (t == 4) return x.sum() > 0.0 && fv > 1e-10 * (1.0 + x.cwiseAbs().maxCoeff());
                    if (t == 3) return fv > 1e-12 || x(ia) - x(ic) > qf_of(x(ic)) * (1.0 + 1e-9);
                    // the corner face only where the ordering holds: an intermediate below the minor is a
                    // swapped ordering (re-read below), not a corner -- taken as one, the set cycled
                    if (t == 5) return !act[ib] && act[3] && x(ib) >= x(ic) - 1e-9 * (1.0 + x.cwiseAbs().maxCoeff()) &&
                                       (fv > 1e-12 || x(ia) - x(ib) > qf_of(x(ib)) * (1.0 + 1e-9));
                    return fv > 1e-12 * (1.0 + x.cwiseAbs().maxCoeff());
                };
                std::array<double, NS> f{};
                std::array<Eigen::Vector3d, NS> n, m;
                Eigen::Matrix3d dm_cap;
                double gpv = g_in, evv = e_in;
                eval(s_tr, dl, f, n, m, gpv, evv, dm_cap);
                const double fscale = 1.0 + s_tr.cwiseAbs().maxCoeff();
                for (int t = 0; t < 3; ++t) act[t] = crossed(t, s_tr, f[t]) ? 1 : 0;
                bool solved = false;
                int retries = 0;
                bool line_tried = false;
                std::array<int, 16> dropped{};
                int ndropped = 0;
                std::array<char, NS> readded{};
                // the cap's multiplier where p_p passes through x (an initial guess: the solution is the set's)
                auto seed_cap = [&](const Eigen::Vector3d& x) {
                    if (!cap_on || !act[4]) return;
                    const double pm = x.sum() / 3.0;
                    if (!(pm > 0.0)) return;
                    Eigen::Vector3d dqc = dqt;
                    if (act[5]) { dqc.setZero(); dqc(ia) = 1.0; dqc(ib) = -0.5; dqc(ic) = -0.5; }
                    const double qt = x.dot(dqc);
                    const double target = std::sqrt(qt * qt / (alpha * alpha) + pm * pm);
                    if (!(target > pc_of(ev0 + hi * 2.0 * pm * dl[4]))) return;
                    dl[4] = (p.cap_ev_from_pc(target) - ev0) / (hi * 2.0 * pm);
                };
                if (!act[0] && !act[1] && !act[2]) {
                    if (crossed(3, s_tr, f[3])) {
                        act[3] = 1;
                    } else if (crossed(4, s_tr, f[4])) {
                        act[4] = 1;
                        seed_cap(s_tr);
                    } else {
                        res.stress = s_tr;
                        res.gp = g_in;
                        res.ev = e_in;
                        res.plastic = false;
                        res.tangent = De;
                        return true;
                    }
                }
                Mat Jlast;
                std::array<int, NS> idx{};
                int nact = 0;
                for (int pass = 0; pass < 24 && !solved; ++pass) {
                    nact = 0;
                    for (int t = 0; t < NS; ++t) if (act[t]) idx[nact++] = t;
                    if (nact == 0) break;
                    if (act[3] && !plateau && sg(ia) - sg(ic) >= qf_of(sg(ic)) * (1.0 - 1e-9) &&
                        s_in(ia) - s_in(ic) < qf_of(s_in(ic)) * (1.0 - 1e-9)) {
                        sg = s_in;
                        dl.fill(0.0);
                        seed_cap(sg);
                    }
                    bool conv = false;
                    const Eigen::Vector3d sg_pass = sg;
                    const std::array<double, NS> dl_pass = dl;
                    const bool plateau_pass = plateau;
                    Eigen::Vector3d sg_ok = sg;
                    for (int it = 0; it < 60; ++it) {
                        eval(sg, dl, f, n, m, gpv, evv, dm_cap);
                        Eigen::Vector3d flow = Eigen::Vector3d::Zero();
                        for (int a = 0; a < nact; ++a) flow += dl[idx[a]] * m[idx[a]];
                        Vec R(3 + nact);
                        R.head<3>() = sg - s_tr + De * flow;
                        for (int a = 0; a < nact; ++a) R(3 + a) = f[idx[a]];
                        Mat J = Mat::Zero(3 + nact, 3 + nact);
                        J.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity() + (act[4] ? dl[4] : 0.0) * De * dm_cap;
                        for (int a = 0; a < nact; ++a) {
                            J.block<3, 1>(0, 3 + a) = De * m[idx[a]];
                            J.block<1, 3>(3 + a, 0) = n[idx[a]].transpose();
                            for (int b = 0; b < nact; ++b) J(3 + a, 3 + b) = dfdl(idx[a], idx[b]);
                        }
                        Jlast = J;
                        const double rn = R.head<3>().norm() / fscale + R.tail(nact).cwiseAbs().maxCoeff() / fscale;
                        if (rn < 1e-13) { conv = true; break; }
                        const Eigen::FullPivLU<Mat> lu(J);
                        if (lu.rank() < 3 + nact) break;
                        const Vec dz = lu.solve(-R);
                        sg_ok = sg;
                        const std::array<double, NS> dl_old = dl;
                        sg += dz.head<3>();
                        for (int a = 0; a < nact; ++a) dl[idx[a]] += dz(3 + a);
                        if (!sg.allFinite()) break;
                        if (act[3] && !plateau) {
                            // a hyperbola iterate past q_f reads no surface: the step is cut back until it is not
                            double th = 1.0;
                            for (int h = 0; h < 40 && sg(ia) - sg(ic) >= qf_of(sg(ic)) * (1.0 - 1e-9); ++h) {
                                th *= 0.5;
                                sg = sg_ok + th * dz.head<3>();
                                dl = dl_old;
                                for (int a = 0; a < nact; ++a) dl[idx[a]] += th * dz(3 + a);
                            }
                            // gamma_p at its failure value: no root of the hyperbola below q_f, the line is the surface
                            const double scf = sg(ic);
                            const double gpf = fb(qf_of(scf) * (1.0 - 1e-12), qa_of(scf));
                            if (qf_of(scf) > 0.0 && g_in + hi * (dl[3] + dl[5]) >= gpf) plateau = true;
                        }
                    }
                    if (!conv) {
                        // the surfaces the last finite iterate crosses join the set, solved again from the pass's start
                        const Eigen::Vector3d xr = sg.allFinite() ? sg : sg_ok;
                        int addr = -1;
                        if (retries < 6 && xr.allFinite()) {
                            std::array<double, NS> fr{};
                            std::array<Eigen::Vector3d, NS> nr, mr;
                            Eigen::Matrix3d dmr;
                            double gr = g_in, er = e_in;
                            eval(xr, dl, fr, nr, mr, gr, er, dmr);
                            for (int t = 0; t < NS && addr < 0; ++t)
                                if (!act[t] && t != 5 && crossed(t, xr, fr[t])) addr = t;
                        }
                        if (addr < 0 && retries < 6)
                            for (int k = ndropped - 1; k >= 0 && addr < 0; --k)
                                if (!act[dropped[k]] && !readded[dropped[k]]) addr = dropped[k];
                        if (addr >= 0) readded[addr] = 1;
                        if (addr >= 0) {
                            sg = sg_pass;
                            dl = dl_pass;
                            plateau = plateau_pass;
                            act[addr] = 1;
                            dl[addr] = 0.0;
                            seed_cap(sg);
                            ++retries;
                            continue;
                        }
                        if (act[3] && !plateau && !line_tried) {
                            line_tried = true;
                            sg = sg_pass;
                            dl = dl_pass;
                            plateau = true;
                            continue;
                        }
                        break;
                    }
                    // Koiter: a negative multiplier leaves the set (and a cap at a tensile mean stress)...
                    // (negative beyond the rounding of the set's own multipliers: at the tip they are 1e-9)
                    // ...but on the edge where both lateral stresses are held by their planes, the flow of the
                    // shear face and of the cap is not one direction: any mix of the faces (a, c) and (a, b) and
                    // of the cap's two orderings is a gradient there (the edge's subdifferential). It moves
                    // multiplier between the two planes and leaves the stress and the hardening as they are:
                    // what the minor's plane is short, the intermediate's gives. (Read with the face (a, c) and
                    // the ordering as solved, the minor's plane is at its least; dropped for it, the set cycled
                    // between the orderings -- measured: soft clay at the surface, pulled sideways under fill.)
                    if (act[ib] && act[ic] && act[3] && !act[5] && dl[ic] < 0.0) {
                        double give = std::max(0.0, -m[3](ic) * dl[3]);
                        const double qt = sg.dot(dqt);
                        if (act[4] && qt > 0.0)
                            give += std::max(0.0, 2.0 * qt / (alpha * alpha) * (2.0 * delta - 1.0) * dl[4]);
                        if (dl[ic] + give >= 0.0 && dl[ib] + dl[ic] >= 0.0) {
                            dl[ib] += dl[ic];
                            dl[ic] = 0.0;
                        }
                    }
                    double lsum = 0.0;
                    for (int a = 0; a < nact; ++a) lsum += std::fabs(dl[idx[a]]);
                    int worst = -1;
                    double wv = -1e-9 * lsum;
                    // At the compression corner the two faces are one surface -- sigma_b = sigma_c and one
                    // gamma_p -- and where the set came to it because the face (a, c) alone turned the
                    // ordering over, only their sum has to be positive: the split between them carries the
                    // cap's corner gradient, averaged here, against the cap's face gradient the point flowed
                    // along before it reached the corner (measured: a split of -21 % of the face's multiplier,
                    // with the stress and gamma_p those of the walk to 1e-4).
                    const bool one_corner = act[3] && act[5] && corner_adds > 0 && dl[3] + dl[5] > 0.0;
                    for (int a = 0; a < nact; ++a)
                        if (!(one_corner && (idx[a] == 3 || idx[a] == 5)) && dl[idx[a]] < wv) {
                            wv = dl[idx[a]];
                            worst = idx[a];
                        }
                    if (worst < 0 && act[4] && !(sg.sum() > 0.0)) worst = 4;
                    if (worst >= 0) {
                        act[worst] = 0; dl[worst] = 0.0;
                        if (worst != 5 && ndropped < 16) dropped[ndropped++] = worst;
                        if (worst == 3) { plateau = false; act[5] = 0; dl[5] = 0.0; }
                        continue;
                    }
                    // ...the intermediate fallen below the minor while the face (a, c) holds alone is the
                    // compression corner, not a new ordering: the face (a, b) joins, and the two faces, sharing
                    // gamma_p, hold the two stresses equal. Read as a new ordering (or as the intermediate's
                    // plane, when it had gone below it too), the set cycled between the two orderings --
                    // measured: fill at the toe of an embankment on 15-noded elements, the walk ending on
                    // the corner with sigma_b = sigma_c = 0.05 kPa, off both planes.
                    if (!act[5] && act[3] && !act[ib] && !act[ic] && corner_adds < 2 &&
                        sg(ib) < sg(ic) - 1e-9 * fscale) {
                        act[5] = 1;
                        dl[5] = 0.0;
                        ++corner_adds;
                        continue;
                    }
                    // ...a surface the answer crosses joins it (a held intermediate takes the corner face away)...
                    std::array<double, NS> f2{};
                    std::array<Eigen::Vector3d, NS> n2, m2;
                    Eigen::Matrix3d dm2;
                    double g2 = g_in, e2 = e_in;
                    eval(sg, dl, f2, n2, m2, g2, e2, dm2);
                    int add = -1;
                    for (int t = 0; t < NS && add < 0; ++t)
                        if (!act[t] && crossed(t, sg, f2[t])) add = t;
                    if (add >= 0) {
                        act[add] = 1; dl[add] = 0.0;
                        if (add == 4) seed_cap(sg);
                        if (add == ib && act[5]) { act[5] = 0; dl[5] = 0.0; }
                        continue;
                    }
                    // ...then the answer's own ordering: the trial's most-pulled stress need not end the smallest
                    const Order o2 = order_of(sg);
                    const double tie2 = 1e-9 * fscale;
                    if (reorders < 4 && (sg(o2.a) > sg(ia) + tie2 || sg(o2.c) < sg(ic) - tie2)) {
                        set_order(o2.a, o2.b, o2.c);
                        ++reorders;
                        continue;
                    }
                    gp_end = g2;
                    ev_end = e2;
                    solved = true;
                }
                // the floor and the ordering assumed must hold at the end; else the walk decides
                if (solved) {
                    const double tie = 1e-9 * fscale;
                    if (!(sg(ia) >= sg(ib) - tie && sg(ib) >= sg(ic) - tie) || sg(ib) > plim) solved = false;
                    // an answer on the line has hardened gamma_p to its failure value on the way
                    if (solved && act[3] && plateau) {
                        const double scf = sg(ic);
                        const double gpf = fb(qf_of(scf) * (1.0 - 1e-12), qa_of(scf));
                        if (g_in < gpf && g_in + hi * (dl[3] + dl[5]) < gpf * (1.0 - 1e-9)) solved = false;
                    }
                }
                if (solved) {
                    res.stress = sg;
                    res.gp = gp_end;
                    res.ev = act[4] ? ev_end : e_in;
                    res.plastic = true;
                    // the consistent tangent: d sigma = (J^-1)_ss De d eps
                    const Eigen::FullPivLU<Mat> lu(Jlast);
                    res.tangent = lu.inverse().topLeftCorner<3, 3>() * De;
                    if (!res.tangent.allFinite()) res.tangent.setZero();
                    return true;
                }
            }
            return false;
        };
        if (ten_on && law_cc > 0.0 && step > 0.0 && !(plan_in && plan_in->n > 0 && plan_in->dT > 0.0)) {
            const Order on = order_of(sigma_n);
            // (a replay of an implicit answer is returned the same way, whether or not its perturbed
            // trial still crosses the cut-off)
            if (sigma_n(on.b) <= plim && ((sigma_n + ds_tr).minCoeff() < -ten || (plan_in && plan_in->implicit))) {
                BeOut r1;
                if (be_return(sigma_n, gamma_p_n, ev_n, dstrain, r1)) {
                    HsIntegrated out;
                    out.stress = r1.stress;
                    out.gamma_p = r1.gp;
                    out.pp = r1.ev == ev_n ? pp_n : pc_of(r1.ev);
                    out.plastic = r1.plastic;
                    out.nsub = 1;
                    out.saturated = 0;
                    out.tangent = r1.tangent;
                    if (plan_out) { plan_out->clear(); plan_out->implicit = true; }
                    return out;
                }
            }
        }
        // Only a COMPRESSIVE trial walks from the tip in graded substeps: the cap closes the
        // compressive side only, so a trial with a tensile mean stress has no cap to grow and no
        // singular first step -- with c = 0 it stays at the apex (the rule above), with c > 0 it
        // is carried elastically up to the cone, walked the usual way. (Until 2026-10-29 such a
        // trial from the tip stayed on it, read as "a cap of size zero admits nothing", and with
        // c > 0 a stress-free point could not carry the tension its cohesion allows.)
        if (cap_on && step > 0.0 && (on_tip || size0 < 0.25 * step) && mean(ds_tr) > 0.0) {
            // PUSHED FROM THE TIP, OR FAR PAST WHAT THE STATE AMOUNTS TO: A GRADED WALK. Each
            // substep is a fraction g of the distance covered, the state's own size counted as
            // covered (never less than kLead of the increment, or from nothing it would never
            // start), so the substeps grow geometrically from a very small first one.
            //
            // THE FIRST SUBSTEP DECIDES THE PATH. With c = 0 the strength is zero at the tip while
            // the stiffness keeps its floor, and the elastic oedometric ratio nu/(1-nu) lies below
            // the active ratio (1-sin phi)/(1+sin phi): the exact path from rest is ON the failure
            // line from its first infinitesimal step. A coarse first substep reads the tip's zero
            // cap gradient, steps elastically, never meets the line, and walks a different path to
            // the end. Measured at the stress point (the corpus oedometer's sand, one increment from
            // rest): 200 equal substeps -- the walk used here until 2026-10-29 -- gave
            // sigma_3/sigma_1 = 0.483 and p_p 8-10% low wherever the increment exceeded 1.6e-4 of
            // strain (the count needed grows with the increment), against 0.271 (= the active
            // ratio) and the 25 600-substep answer to 1e-7 for this walk. On the oedometer from
            // rest (KV-NUM-009) the 200-400 kPa range sat 0.35 points off the other two; with this
            // walk the three agree to every printed digit.
            //
            // Its error is measured and refined like the uniform walk's below: the largest R of its
            // substeps, and g shrunk by (0.81 STOL / R)^(1/p) until R is under 0.81 STOL -- the
            // factor is 1 at that threshold, so the walk is a continuous function of the strain.
            constexpr double kLead = 1.0e-3;
            // FROM THE TIP ITSELF the walk starts at an isotropic, normally consolidated state of
            // vanishing size -- p = p_p = delta, q = 0, on both surfaces -- of which the stress-free
            // state is the limit: delta a millionth of the elastic trial, never below 1e-9 p_ref (at
            // 1e-8 of a trial of 0.03 kPa the seed fell under the round-off of the walk's own tests and
            // the failure line was taken again). Measured on one oedometric increment from rest: the
            // seeds 1e-4, 1e-6 and 1e-8 of the trial give the same branch at every increment size and
            // tolerance, and the last two the same stress to 4e-6.
            constexpr double kSeed = 1.0e-6;
            const bool seeded = on_tip;
            const double delta = seeded ? std::max(kSeed * step, 1e-9 * pr) : 0.0;
            const Eigen::Vector3d s_start = seeded ? Eigen::Vector3d::Constant(delta) : sigma_n;
            const double ev_start = seeded ? p.cap_ev_from_pc(delta) : ev_n;
            // The substeps a growth g takes to cover the increment: done = covered ((1+g)^k - 1).
            auto span_of = [](double cov) { return std::log((1.0 + cov) / cov); };
            struct Lead {
                Eigen::Vector3d sig;
                double gp = 0.0, ev = 0.0;
                Eigen::Matrix3d tangent;
                bool plastic = false, ok = false;
                int walked = 0;
                double r_max = 0.0;
                double r_first = 0.0, r_rest = 0.0;   // R of the first substep, and of the others
            };
            auto lead_walk = [&](double g, double covered) {
                Lead l;
                l.sig = s_start; l.gp = gamma_p_n; l.ev = ev_start; l.tangent = k_n.De;
                double done = 0.0;
                const int budget = (int)std::ceil(span_of(covered) / std::log1p(g)) + 2;
                bool fine = l.sig.allFinite() && std::isfinite(l.gp) && std::isfinite(l.ev);
                while (fine && done < 1.0 - 1e-12 && l.walked < budget) {
                    const double dT = std::min(g * (covered + done), 1.0 - done);
                    const Eigen::Vector3d de = dT * dstrain;
                    double r_one = 0.0;
                    substep(l.sig, l.gp, l.ev, de, l.tangent, l.plastic, r_one, true, nullptr, nullptr);
                    l.r_max = std::max(l.r_max, r_one);
                    (l.walked == 0 ? l.r_first : l.r_rest) = std::max(l.walked == 0 ? l.r_first : l.r_rest, r_one);
                    done += dT;
                    ++l.walked;
                    fine = l.sig.allFinite() && std::isfinite(l.gp) && std::isfinite(l.ev);
                }
                l.ok = fine && done >= 1.0 - 1e-12;
                return l;
            };
            // TWO THINGS ARE REFINED, EACH BY ITS OWN ERROR. The first substep leaves the tip itself,
            // where nothing is self-similar -- the cap's gradient is zero there and the pair disagree
            // by the whole of what the substep does -- so its R is FIRST order in its own size,
            // g covered, and falls only as that does. Every later substep starts from a state the
            // size of what it covers and is the usual second order in g. Refining g alone for the
            // first one made every substep of the walk smaller for the sake of one (measured: 8 341
            // substeps where 1 100 met the tolerance once the first was shrunk instead); shrinking
            // the lead costs only ln(1/f)/g substeps more. So the lead is shrunk by the first
            // substep's R (first order), and g by the rest's (the order measured, as below).
            // The lead starts at the state's own size (never less than kLead of the increment, or
            // from nothing the walk would never start); a lead shrunk below the state's size only
            // makes the first substeps smaller.
            const double cov_min = 1e-12;
            double covered = seeded ? delta / step : std::max(size0 / step, kLead);
            // A fiftieth of the distance covered to start with (at a tenth the walk was 0.7% off
            // on the oedometer from rest), and never finer than the ceiling allows.
            auto g_min_of = [&](double cov) { return std::expm1(span_of(cov) / kMaxSubsteps); };
            double g = std::max(0.02, g_min_of(covered));
            Lead l = lead_walk(g, covered);
            bool lead_met = false;
            double g_prev = 0.0, r_prev = 0.0;
            for (int it = 0; it < kMaxWalks && l.ok && l.r_max > kAim * tol; ++it) {
                bool refined = false;
                if (l.r_first > kAim * tol && covered > cov_min) {
                    covered = std::max(covered * (kAim * tol / l.r_first), cov_min);
                    refined = true;
                }
                const double g_min = g_min_of(covered);
                if (l.r_rest > kAim * tol && g > g_min) {
                    double order = 2.0;   // Euler's local error, R ~ (g covered)^2
                    if (g_prev > 0.0 && r_prev > l.r_rest)
                        order = std::clamp(std::log(r_prev / l.r_rest) / std::log(g_prev / g), 1.0, 2.0);
                    else if (g_prev > 0.0)
                        order = 1.0;
                    g_prev = g;
                    r_prev = l.r_rest;
                    g = std::max(g * std::pow(kAim * tol / l.r_rest, 1.0 / order), g_min);
                    refined = true;
                }
                g = std::max(g, g_min);
                if (!refined) break;
                const Lead l_prev = l;
                l = lead_walk(g, covered);
                // Two successive walks that agree to the tolerance have met it (see the uniform
                // walk below, THE ANSWER IS WHAT THE TOLERANCE IS FOR).
                if (l.ok && l_prev.ok && l.r_max > kAim * tol) {
                    const double scale = std::max(l.sig.norm(), pr / 20.0);
                    if ((l.sig - l_prev.sig).norm() <= kAim * tol * scale &&
                        std::fabs(l.gp - l_prev.gp) <= kAim * tol * std::max(std::fabs(l.gp), 1e-12) &&
                        std::fabs(l.ev - l_prev.ev) <= kAim * tol * std::max(std::fabs(l.ev), 1e-12)) {
                        lead_met = true;
                        break;
                    }
                }
            }
            if (l.ok) {
                HsIntegrated led;
                led.stress = l.sig;
                led.gamma_p = l.gp;
                led.pp = p.cap_pc_from_ev(l.ev);
                led.tangent = l.tangent;
                led.plastic = l.plastic;
                led.nsub = l.walked;
                led.saturated = l.r_max > tol && !lead_met ? 1 : 0;
                if (plan_out) plan_out->clear();
                return led;
            }
            // The graded walk failed (a non-finite state): the increment is walked the usual way.
        }
    }
    // The committed state's elastic matrix, until a plastic substep replaces it with its own
    // continuum tangent (an elastic increment returns this one).
    Eigen::Matrix3d tangent = k_n.De;

    if (plan_out) plan_out->clear();

    const bool replaying = plan_in && plan_in->n > 0 && plan_in->dT > 0.0;

    // --- WHAT THE TOLERANCE MEANS, AND HOW MANY SUBSTEPS MEET IT ------------------------------
    //
    // THE TOLERANCE IS SLOAN'S (Sloan, Abbo & Sheng 2001, Engineering Computations 18:121-154,
    // Sec. 3.4): every substep's relative local error, half the difference of its Euler and
    // modified-Euler increments over the stress it ends at (their Eqn 24, the stress part), is at
    // most STOL. That is what a substep is accepted on there, and what they measure the stress
    // error against (Tables II-V: below STOL in every case). Until 2026-10-28 the tolerance had no
    // one meaning here: the count was estimated as sqrt(STOL dT^3 / e) -- about sqrt(n) more
    // substeps than this law asks, as a margin for substeps that cross a surface, a corner or the
    // plateau -- and the warning said "ran out" whenever the ceiling clipped that estimate, without
    // measuring anything. Measured on a staged excavation in undrained clay over sand (14491
    // clipped calls, 604 replayed against 40 000 substeps): the modified-Euler result met the
    // tolerance at 200 substeps in 590, while almost every one failed the per-substep test.
    //
    // Two departures from the paper are kept, both forced by measurement.
    // (1) The substeps are EQUAL and their count is chosen for the whole increment, not grown and
    //     shrunk substep by substep with rejections. Accept/reject makes the subdivision a
    //     discontinuous function of the strain increment: neighbouring Newton iterates integrate
    //     along different subdivisions and the global iteration cannot converge below that noise
    //     (measured on the corpus oedometer, 2026-08-20: the residual stalled at 1.2e-3, the line
    //     search halving to 2.4e-4 without descent).
    // (2) The count is a REAL number, with a partial last substep. With ceil() the answer jumps
    //     where the count does (the same oedometer stalled at 5.7e-6); with a real count it is a
    //     CONTINUOUS function of the strain increment.
    // Sloan's own rule is therefore applied to the walk as a whole: walk, measure the largest R
    // of the substeps taken, and if it is above 0.81 STOL -- the 0.9 of his Eqn 27 on the step,
    // squared onto R -- walk again with n (R / 0.81 STOL)^(1/p). The order p is his 2 (Euler's
    // local error ~ dT^2) for the first refinement and the one the two walks measured after it,
    // held to [1, 2]: a substep that meets a corner or the plateau is first order, which is what
    // Sloan's smoothed surfaces never meet and what the old margin was guessing at. At the
    // threshold the factor is 1, so neither the count nor the answer jumps there.
    //
    // A FOURTH-ORDER STEP WAS TRIED HERE AND TAKEN OUT (2026-10). The classical four-stage step
    // with a step-doubling estimate of its own error needs a tenth of the substeps and ran the
    // boundary-value benches twice as fast with the same answers, but on the oedometer path from
    // rest (KV-NUM-009) the deviation from the idealised law, which the power law's
    // self-similarity makes the same at every stress level, varied by 0.74 points with it against
    // 1e-4 with the pair: the stages of one step fall on different branches wherever the active
    // set changes inside it, and the order is lost there.
    double nreal = 1.0;
    // The estimate's first pass over the whole increment is the walk's own step when the
    // increment needs no subdivision -- the same two evaluations at the same state, strain and
    // stiffness. Kept and taken as that step instead of being made again: the integration of a
    // point that stays elastic, or nearly so, cost twice what it does.
    Inc whole_k1, whole_k2;
    bool have_whole = false;
    if (replaying) {
        // Replay: the subdivision is uniform, so its one step size reproduces it exactly.
        nreal = 1.0 / plan_in->dT;
    } else if (dstrain_w.squaredNorm() > 0.0) {
        // The first estimate, from the first substep only: Sloan's step law, twice -- over the
        // whole increment, then over the substep that proposed, where R ~ dT^2 holds better.
        // No plasticity gate: E_ur moves with sigma3, so even an elastic step has an error.
        for (int pass = 0; pass < 2; ++pass) {
            const double dT = 1.0 / nreal;
            const Eigen::Vector3d de = dT * dstrain_w;
            double r = 0.0;
            if (elastic_fraction(s0, gp0, ev0, de) < 0.0) {
                const Inc f1 = increment(s0, gp0, ev0, de);
                const Inc f2 = increment_core(s0 + f1.dsig, gp0 + f1.dgp, ev0 + f1.dev, de, activated(f1));
                // Only an ELASTIC pair is kept: a plastic one was evaluated without the walk's own
                // treatment of the events inside the substep (substep()), so the call would be
                // integrated one way and its replays -- the consistent tangent's perturbed runs,
                // which have no estimate -- another, and their difference divided by the
                // perturbation.
                if (pass == 0 && dT == 1.0 && !f1.plastic && !f2.plastic) {
                    whole_k1 = f1; whole_k2 = f2; have_whole = true;
                }
                r = rel_error(f1.dsig, f2.dsig, s0 + 0.5 * (f1.dsig + f2.dsig));
            } else {   // the first substep crosses a surface: measured as the walk will take it
                Eigen::Vector3d sc = s0;
                double gc = gp0, ec = ev0;
                Eigen::Matrix3d tc = tangent;
                bool pc = false;
                substep(sc, gc, ec, de, tc, pc, r, false, nullptr, nullptr);
            }
            if (!(r > tol)) break;
            nreal = std::min(nreal * std::sqrt(r / tol), (double)kMaxSubsteps);
        }
    }

    struct Walk {
        Eigen::Vector3d sig;
        double gp = 0.0, ev = 0.0;
        Eigen::Matrix3d tangent;
        bool plastic = false;
        int taken = 0;
        double r_max = 0.0;   // the largest local error of the substeps taken
    };
    auto walk = [&](double n) {
        Walk w;
        w.sig = s0; w.gp = gp0; w.ev = ev0; w.tangent = tangent;
        const double dT_full = 1.0 / n;
        if (plan_out) { plan_out->clear(); plan_out->dT = dT_full; }
        // The walk's own count bounds the loop, never the ceiling: a walk planned past the ceiling
        // and stopped at it would return the integral of part of the strain, silently.
        const int planned = (int)std::ceil(n - 1e-9) + 1;
        for (double T = 0.0; T < 1.0 - 1e-12 && w.taken < planned; ++w.taken) {
            if (T > 0.0 && ends_at_apex(w.sig, w.gp, (1.0 - T) * dstrain_w)) {   // (above)
                w.sig.setZero();
                w.tangent.setZero();
                w.plastic = true;
                break;
            }
            const double dT = std::min(dT_full, 1.0 - T);
            const Eigen::Vector3d de = dT * dstrain_w;
            const bool whole = have_whole && w.taken == 0 && dT == 1.0;
            substep(w.sig, w.gp, w.ev, de, w.tangent, w.plastic, w.r_max, false,
                    whole ? &whole_k1 : nullptr, whole ? &whole_k2 : nullptr);
            if (plan_out) ++plan_out->n;
            T += dT;
        }
        return w;
    };

    Walk w = walk(nreal);
    bool answer_met = false;   // two successive walks agreed (see below)
    if (!replaying) {
        double n_prev = 0.0, r_prev = 0.0;
        for (int it = 0; it < kMaxWalks && w.r_max > kAim * tol && nreal < (double)kMaxSubsteps; ++it) {
            double order = 2.0;
            if (n_prev > 0.0 && r_prev > w.r_max && nreal > n_prev)
                order = std::clamp(std::log(r_prev / w.r_max) / std::log(nreal / n_prev), 1.0, 2.0);
            else if (n_prev > 0.0)
                order = 1.0;   // the error did not fall with the count: refine as for first order
            n_prev = nreal;
            r_prev = w.r_max;
            nreal = std::min(nreal * std::pow(w.r_max / (kAim * tol), 1.0 / order), (double)kMaxSubsteps);
            const Walk w_prev = w;
            w = walk(nreal);
            // THE ANSWER IS WHAT THE TOLERANCE IS FOR. R bounds the local error of each substep's
            // Euler step (Sloan's test, conservative by design); the walk's own result, refined, is
            // the error that is handed on. Where the flow slides along a surface it is neutral to --
            // the shear planes held at a zero multiplier while the cap does the work -- the two
            // evaluations of a pair read the two sides of that zero, R falls only as 1/n while the
            // mean is already the sliding rate: measured on the Berlin Sand triaxial at STOL 1e-8,
            // the answer at 13 substeps within 2.5e-10 of the one at 1600, R 1.9e-5 there and still
            // 5e-8 at the ceiling of 10 000 (the test 33 s -> 50 min). So two successive walks
            // whose answers -- stress and both hardening variables -- agree to the tolerance in the
            // same measure R uses have met it, R notwithstanding, and the finer one is kept.
            if (w.r_max > kAim * tol) {
                const double scale = std::max(w.sig.norm(), pr / 20.0);
                const bool same =
                    (w.sig - w_prev.sig).norm() <= kAim * tol * scale &&
                    std::fabs(w.gp - w_prev.gp) <= kAim * tol * std::max(std::fabs(w.gp), 1e-12) &&
                    std::fabs(w.ev - w_prev.ev) <= kAim * tol * std::max(std::fabs(w.ev), 1e-12);
                if (same) { answer_met = true; break; }
            }
        }
    }
    any_plastic = any_plastic || w.plastic;
    // Saturation is REPORTED, never absorbed: a walk that left a substep above the tolerance --
    // because the ceiling or the walk budget stopped the refinement -- did not meet it, and the
    // caller is told so rather than handed a number that looks like every other one. A replay
    // measures nothing and reports nothing.
    const int clipped = (!replaying && w.r_max > tol && !answer_met) ? 1 : 0;
    HsIntegrated out;
    out.stress = w.sig;
    out.gamma_p = w.gp;
    out.pp = cap_on ? p.cap_pc_from_ev(w.ev) : pp_n;
    out.tangent = w.tangent;
    out.plastic = any_plastic;
    out.nsub = w.taken;
    out.saturated = clipped;
    return out;
}

HsSeed hs_seed_from_history(const HardeningSoilParams& p, double sxx, double syy,
                            double sxy, double szz, int mode, double OCR, double POP) {
    auto principal = [](double xx, double yy, double xy, double zz) {
        const double m = 0.5 * (xx + yy), r = std::sqrt(0.25 * (xx - yy) * (xx - yy) + xy * xy);
        return Eigen::Vector3d(m + r, m - r, zz);
    };
    const double k0nc = p.K0nc > 0.0 ? p.K0nc : 1.0 - std::sin(p.friction);
    double syy_c = syy;
    if (mode == 1) syy_c = std::max(1.0, OCR) * syy;
    else if (mode == 2) syy_c = syy + std::max(0.0, POP);
    const Eigen::Vector3d now = principal(sxx, syy, sxy, szz);
    const Eigen::Vector3d past = principal(k0nc * syy_c, syy_c, sxy, k0nc * syy_c);
    // THE PRE-CONSOLIDATION STATE IS WHAT BOTH SURFACES REMEMBER: the ground has been at that stress,
    // on the K0nc line, and both the cap and the shear surface were pushed out to it. An
    // overconsolidated state lies inside them and answers elastically up to where it left them.
    // Measured on the undrained element tests of four clays against the published curves (rms of the
    // deviators and mean stresses 3.5% -> 3.1%, first yield of a clay of POP 32 at 107 kPa against the
    // published 112 and 103 with the shear surface on the present state), and on a strutted excavation
    // in clay with a pre-overburden of 5 kPa (wall movement after each stage, published 6.1 / 24.1 / 57.1
    // / 141.6 mm): 6.04 / 23.11 / 57.06 / 156.34 against 6.11 / 25.01 / 59.73 / 162.65. This same
    // excavation had once chosen the other way (22.61 / 55.50 against 23.68 / 57.25), but on a joint that
    // read the wrong water.
    return {std::max(hs_initial_pp(p, now), hs_initial_pp(p, past)),
            std::max(hs_initial_gamma_p(p, now), hs_initial_gamma_p(p, past))};
}

void hs_oedometer_probe(const HardeningSoilParams& p, double& Eoed_pref,
                        double& K0) {
    // The walk STARTS ON THE NORMAL-CONSOLIDATION LINE it is asked about -- K0nc, on the cap,
    // with the gamma_p of that state -- and above the stiffness floor. It used to start isotropic
    // at 2% of p_ref and assume the path had forgotten that by p_ref; it had not: at alpha = 2.4
    // the K0 read at p_ref was 0.405 from a 2 kPa start and 0.453 from a 20 kPa one, so every
    // calibrated cap carried an arbitrary starting point. From the line itself the calibration's
    // root is exact whatever the walk's length: the right alpha keeps the path at K0nc, a wrong
    // one drifts towards its own K0, and the sign of the drift is all the bisection reads.
    const double pr = p.p_ref;
    const double k0s = p.K0nc > 0.0 ? p.K0nc : 1.0 - std::sin(p.friction);
    const double s10 = std::max(0.3 * pr, 0.15 * pr / k0s);
    Eigen::Vector3d sig(s10, k0s * s10, k0s * s10);
    double gp = hs_initial_gamma_p(p, sig), pp = hs_initial_pp(p, sig), s1b = s10;
    // The walk up to the reference pressure is taken in large steps (the integrator's own error
    // control sets the accuracy, not the outer step), and the READING is taken AT p_ref: the step
    // that crosses it is repeated to land on p_ref, and Eoed is the response to a strain
    // increment a hundredth of that. It used to be the secant over the crossing step, which
    // spans about 6% of p_ref: with Eoed rising as sigma^m that secant is the tangent a few
    // percent above p_ref, and every calibrated cap came out that much too soft -- measured as a
    // +1.7% oedometric settlement against the closed form on the corpus sand.
    const double de_coarse = 2.0e-3, de_fine = 2.0e-4;
    Eoed_pref = 0.0; K0 = 0.0;
    for (int i = 0; i < 3000; ++i) {
        const double de1 = (sig(0) < 0.9 * pr) ? de_coarse : de_fine;
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, 0, 0),
                                            kHsCalibrationTol);
        if (s1b < pr && r.stress(0) >= pr) {
            // Land on p_ref: the crossing step, shortened to the fraction that reaches it, twice
            // (the second pass corrects the curvature the first one assumed away).
            double de_to = de1 * (pr - s1b) / (r.stress(0) - s1b);
            HsIntegrated at = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de_to, 0, 0),
                                           kHsCalibrationTol);
            for (int pass = 0; pass < 2; ++pass) {
                const double slope = (at.stress(0) - s1b) / de_to;
                de_to += (pr - at.stress(0)) / slope;
                at = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de_to, 0, 0), kHsCalibrationTol);
            }
            const double de_t = 1e-2 * de_fine;   // secant curvature ~1e-4, far above the drift residual
            const HsIntegrated tip = hs_integrate(p, at.stress, at.gamma_p, at.pp,
                                                  Eigen::Vector3d(de_t, 0, 0), kHsCalibrationTol);
            Eoed_pref = (tip.stress(0) - at.stress(0)) / de_t;
            K0 = at.stress(2) / at.stress(0);
            return;
        }
        s1b = r.stress(0); sig = r.stress; gp = r.gamma_p; pp = r.pp;
        if (sig(0) > 1.25 * pr) break;
    }
}

// E50 as a drained triaxial test reads it: from an isotropic stress at p_ref on the cap, the minor stress held, the
// secant at half the failure deviator.
static double hs_triaxial_e50(const HardeningSoilParams& p) {
    const double s3 = p.p_ref, qf = p.q_failure(s3);
    if (!(qf > 0.0)) return p.E50(s3);
    Eigen::Vector3d sig(s3, s3, s3);
    double gp = 0.0, pp = hs_initial_pp(p, sig);
    double e1 = 0.0, q_prev = 0.0, e_prev = 0.0;
    const double de1 = 0.04 * 0.5 * qf / p.E50(s3);   // 25 steps to the strain of the hyperbola alone
    double lateral = -0.2;                            // lateral over axial strain, carried from step to step
    for (int i = 0; i < 200; ++i) {
        double de3 = lateral * de1;
        const auto walk = [&](double d3) {
            return hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, d3, d3), kHsCalibrationTol);
        };
        HsIntegrated r = walk(de3);
        double err = r.stress(1) - s3;
        if (std::fabs(err) >= 1e-7 * s3) {
            const double slope = (walk(de3 + 1e-3 * de1).stress(1) - r.stress(1)) / (1e-3 * de1);
            for (int it = 0; it < 12 && std::fabs(slope) > 0.0 && std::fabs(err) >= 1e-7 * s3; ++it) {
                de3 -= err / slope;
                r = walk(de3);
                err = r.stress(1) - s3;
            }
        }
        lateral = de3 / de1;
        sig = r.stress; gp = r.gamma_p; pp = r.pp;
        e1 += de1;
        const double q = sig(0) - sig(1);
        if (q >= 0.5 * qf) return 0.5 * qf / (e_prev + (0.5 * qf - q_prev) / (q - q_prev) * (e1 - e_prev));
        q_prev = q; e_prev = e1;
    }
    return p.E50(s3);
}

void hs_calibrate_cap(HardeningSoilParams& p, double K0_NC) {
    p.K0nc = K0_NC;   // the line the oedometer probe starts on
    // ONE CALIBRATION PER SOIL, NOT PER PHASE. The two nested bisections below walk an oedometer
    // path some three thousand times, each walk at the ceiling of substeps, and the material
    // models are built afresh for every phase of a job: measured on a tunnel analysis of six
    // phases with two such soils, 37 calibrations of the same two soils took 45% of all the
    // constitutive substeps of the run (37% on a single strip load). The answer is a function
    // of the parameters alone, so it is kept. The key is every member the walk can read; the
    // size check is there so that a member added to the struct is added here as well.
    // cap_cc is not a parameter but the integrator's copy of c.cot(phi), derived from two that are
    // keyed (and NaN on every set of parameters that reaches this function), so it is not keyed.
    // sigma_t is not keyed either: the material wrapper sets it on its own copy for each call, and
    // the calibration walks a compressive oedometer path the cut-off never reaches.
    static_assert(sizeof(HardeningSoilParams) == 21 * sizeof(double),
                  "a member was added to HardeningSoilParams: add it to the calibration key");
    using Key = std::array<double, 16>;
    const Key key{p.E50_ref, p.Eur_ref, p.Eoed_ref, p.m, p.p_ref, p.cohesion, p.friction,
                  p.dilatancy, p.dilatancy_cut ? 1.0 : 0.0, p.Rf, p.nu_ur, p.G0_ref, p.gamma07,
                  p.hardening_scale, p.elastic_scale, K0_NC};
    static std::mutex cache_mutex;
    static std::vector<std::pair<Key, std::array<double, 3>>> cache;
    {
        const std::lock_guard<std::mutex> lock(cache_mutex);
        for (const auto& e : cache)
            if (e.first == key) { p.cap_alpha = e.second[0]; p.cap_beta = e.second[1]; p.ei_scale = e.second[2]; return; }
    }
    // Stored on every way out of the search below but an exception: a calibration stopped half way
    // -- a run cancelled inside it (katai/materials/interrupt.hpp) -- holds an alpha and a beta that
    // are not this material's, and every later build in the process would have read them back.
    struct Keep {
        HardeningSoilParams& p; const Key& key;
        const int unwinding = std::uncaught_exceptions();
        ~Keep() {
            if (std::uncaught_exceptions() > unwinding) return;
            const std::lock_guard<std::mutex> lock(cache_mutex);
            for (const auto& e : cache) if (e.first == key) return;
            cache.push_back({key, {p.cap_alpha, p.cap_beta, p.ei_scale}});
        }
    } keep{p, key};
    const double Kp = hs_cap_Kp(p, K0_NC);
    // THE TWO ROOTS ARE FOUND BY REGULA FALSI (the Illinois variant), not by bisection. Each
    // evaluation is an oedometer walk at the calibration's own integration tolerance, and the
    // bisections took about 290 of them -- 4 of the 6 s the initial phase of a staged excavation
    // with two Hardening Soils spent before its first equilibrium step. Both functions are smooth
    // and monotone in the bracket (Eoed rises with k; K0 falls with alpha, read here in log alpha
    // where it is closer to straight), so the false position, with the Illinois halving that keeps
    // it from stalling on one end, needs a fraction of the walks to reach a root far inside the
    // brackets the bisections stopped at (1e-4 on k, 2e-3 on alpha). The roots are the same
    // roots; they are now found to 1e-6 on k and 1e-5 on log alpha rather than to the midpoint of
    // the last bracket.
    const auto illinois = [](auto&& f, double a, double fa, double b, double fb, double xtol,
                             int maxit) {
        // f(a) and f(b) of opposite sign; returns the root's estimate.
        double c = a;
        int side = 0;
        for (int it = 0; it < maxit; ++it) {
            const double c_prev = c;
            c = (a * fb - b * fa) / (fb - fa);
            if (!std::isfinite(c)) c = 0.5 * (a + b);
            const double fc = f(c);
            if (fc == 0.0 || (it > 0 && std::fabs(c - c_prev) < xtol)) return c;
            if ((fc > 0.0) == (fb > 0.0)) {
                b = c; fb = fc;
                if (side == -1) fa *= 0.5;
                side = -1;
            } else {
                a = c; fa = fc;
                if (side == +1) fb *= 0.5;
                side = +1;
            }
            if (std::fabs(b - a) < xtol) return c;
        }
        return c;
    };
    // Each root starts next to the last one where it brackets there, and that is most of the walks: k moves little with
    // alpha, and alpha little with the initial stiffness.
    double k_last = 0.0, alpha_last = 0.0;
    // Inner: at a given alpha, the k (beta = p_ref / (k Kp)) that gives Eoed_ref; returns K0.
    auto solve_k_return_K0 = [&](double alpha) {
        p.cap_alpha = alpha;
        const auto eoed_gap = [&](double k) {
            p.cap_beta = p.p_ref / (k * Kp);
            double e, k0; hs_oedometer_probe(p, e, k0);
            return e - p.Eoed_ref;   // rises with k (a stiffer cap)
        };
        double klo = 0.2, khi = 5.0, glo = 0.0, ghi = 0.0;
        bool near = false;
        if (k_last > 0.0) {
            const double a = std::max(klo, 0.8 * k_last), b = std::min(khi, 1.25 * k_last);
            glo = eoed_gap(a);
            ghi = eoed_gap(b);
            near = (glo > 0.0) != (ghi > 0.0);
            if (near) { klo = a; khi = b; }
        }
        if (!near) { glo = eoed_gap(klo); ghi = eoed_gap(khi); }
        const double k = (glo > 0.0) == (ghi > 0.0)
                             ? (std::fabs(glo) < std::fabs(ghi) ? klo : khi)
                             : illinois(eoed_gap, klo, glo, khi, ghi, 1e-6, 60);
        k_last = k;
        p.cap_beta = p.p_ref / (k * Kp);
        double e, k0; hs_oedometer_probe(p, e, k0);
        return k0;
    };
    // Outer: alpha for K0 = K0_NC (K0 falls with alpha), in log alpha.
    const auto solve_cap = [&] {
        const double alo = 0.1, ahi = 80.0;
        const auto k0_gap = [&](double la) { return solve_k_return_K0(std::exp(la)) - K0_NC; };
        if (alpha_last > 0.0) {
            const double a = std::max(alo, alpha_last / 1.3), b = std::min(ahi, alpha_last * 1.3);
            const double ga = solve_k_return_K0(a) - K0_NC, gb = solve_k_return_K0(b) - K0_NC;
            if (ga * gb <= 0.0) {
                alpha_last = std::exp(illinois(k0_gap, std::log(a), ga, std::log(b), gb, 1e-5, 60));
                solve_k_return_K0(alpha_last);
                return;
            }
        }
        const double klo = solve_k_return_K0(alo), khi = solve_k_return_K0(ahi);
        if ((K0_NC - klo) * (K0_NC - khi) <= 0.0) {  // bracketed -> reachable
            alpha_last = std::exp(illinois(k0_gap, std::log(alo), klo - K0_NC, std::log(ahi), khi - K0_NC, 1e-5, 60));
        } else {
            alpha_last = std::fabs(klo - K0_NC) < std::fabs(khi - K0_NC) ? alo : ahi;
        }
        solve_k_return_K0(alpha_last);
    };
    solve_cap();
    // THE INITIAL STIFFNESS IS AN INTERNAL PARAMETER LIKE THE OTHER TWO. "In double hardening situations, i.e. both
    // yield loci are hardened simultaneously, analytical back-calculation of internal model parameters is impossible.
    // Therefore, the internal parameters [E_i, H and alpha] are solved for in an iterative scheme so that the HS model
    // simulates the user input E50 in triaxial element tests and both, E_oed, and K0nc in oedometer element tests"
    // (Benz 2007, p. 81). The hyperbola's 2 E50 / (2 - R_f) is the answer only with the cap far away. In a normally
    // consolidated test the cap flows as well, and until this was solved for E50 came out low by what the cap added:
    // read in the test, 0.75 of the value entered on a clay of E_oed / E50 = 0.8, 0.87 on a sand of 1, and on the
    // published points of a drained test of Berlin sand III 389.8 / 502.7 / 576.8 kPa against 400 / 512 / 587 at 0.1 /
    // 0.2 / 0.3% axial strain. The test is from an isotropic stress at p_ref, on the cap, as the definition of E50 has it.
    // The cap is solved again with the stiffness found, since its own walk feels the shear surface, and two rounds settle.
    {
        HardeningSoilParams t = p;
        t.G0_ref = 0.0;
        t.hardening_scale = 1.0;
        t.elastic_scale = 1.0;
        for (int round = 0; round < 2; ++round) {
            t.cap_alpha = p.cap_alpha;
            t.cap_beta = p.cap_beta;
            const auto gap = [&](double f) {
                t.ei_scale = f;
                return hs_triaxial_e50(t) - p.E50_ref;
            };
            const double flo = 1.0, fhi = HardeningSoilParams::kMaxEiScale;
            const double glo = gap(flo), ghi = glo < 0.0 ? gap(fhi) : 0.0;
            p.ei_scale = !(glo < 0.0) ? 1.0 : ghi < 0.0 ? fhi : illinois(gap, flo, glo, fhi, ghi, 1e-4, 30);
            solve_cap();
        }
    }
}

}  // namespace katai::core
