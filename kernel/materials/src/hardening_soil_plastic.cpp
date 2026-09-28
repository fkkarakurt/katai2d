#include <katai/materials/hardening_soil_plastic.hpp>

#include <array>

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

HsIntegrated hs_integrate(const HardeningSoilParams& p,
                          const Eigen::Vector3d& sigma_n, double gamma_p_n,
                          double pp_n, const Eigen::Vector3d& dstrain,
                          double stol,
                          HsSubstepPlan* plan_out,
                          const HsSubstepPlan* plan_in) {
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
    // The stiffness and strength laws read the MINOR principal stress (Benz 2007, Eqn 7.19), and
    // the integrator's vector is not assumed sorted: see Order below.
    // The laws of HardeningSoilParams (stiffness_factor, E50/Eur/Ei, q_failure/q_asymptote)
    // written out with their trigonometry evaluated once per call of the integrator rather than
    // six times per evaluation -- the same operations in the same order, so the same bits.
    const double law_cc = p.cohesion * std::cos(p.friction);
    const double law_s = std::sin(p.friction);
    const double law_den = law_cc + p.p_ref * law_s;
    auto stiff_at = [&](const Eigen::Vector3d& s) {
        Stiff k;
        k.s3 = std::max(s.minCoeff(), plim);
        const double g = std::pow((law_cc + k.s3 * law_s) / law_den, p.m);
        k.Eur = p.Eur_ref * g;
        k.Ei = 2.0 * (p.E50_ref * g) / (2.0 - p.Rf);
        k.qf = 2.0 * (law_cc + k.s3 * law_s) / (1.0 - law_s);
        k.qa = k.qf / p.Rf;
        k.De = detail::hs_elastic(k.Eur, nu);
        return k;
    };
    const Stiff k_n = stiff_at(sigma_n);   // the committed state, for the elastic default
    const bool cap_on = p.cap_beta > 0.0;
    const double ev_n = cap_on ? p.cap_ev_from_pc(pp_n) : 0.0;

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
    auto fcap = [&](const Eigen::Vector3d& s, double ppv) {
        const double qt = qtilde(s), pm = mean(s);
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
    auto spm_of = [&](double q, const Stiff& k) { return dil.from_q(q, k.s3); };
    auto pc_of = [&](double evv) { return cap_on ? p.cap_pc_from_ev(evv) : pp_n; };
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
    auto plane = [&](const Eigen::Vector3d& s, int i, int j, double spm, bool at_fail,
                     const Stiff& k) {
        Plane pl;
        pl.n = Eigen::Vector3d::Zero();
        pl.n(i) = 0.5 * (1.0 - spm);
        pl.n(j) = -0.5 * (1.0 + spm);
        pl.m = Eigen::Vector3d::Zero();
        if (at_fail) {
            // The failure strength reads the minor stress through the same floor as the
            // stiffness (see Stiff), so below it the strength does not move with sigma_j.
            pl.m(i) = 1.0;
            pl.m(j) = -1.0 - (s(j) > plim ? dqf_ds : 0.0);
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
            if (s3 > plim) {
                // In closed form. E_i and E_ur scale with the same stiffness factor g(sigma3),
                // dg/dsigma3 = g m sin(phi)/(c cos(phi) + sigma3 sin(phi)), and q_a = q_f/R_f is
                // linear in sigma3; fbar = A - 2q/E_ur with A = (2/E_i) q/(1 - q/q_a).
                const double sphi = std::sin(p.friction);
                const double r = p.m * sphi / (p.cohesion * std::cos(p.friction) + s3 * sphi);
                const double dqa = 2.0 * sphi / ((1.0 - sphi) * p.Rf);
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
        int corner = kFace;    // which part of the hexagon the step flowed on
        // At a corner, the two principal DIRECTIONS whose stresses are equal (-1: on a face). The
        // corner TYPE only means something relative to the ordering it was found in; the pair is
        // what carries over to a state whose ordering has changed.
        int tie_i = -1, tie_j = -1;
    };
    auto increment_core = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                              const Eigen::Vector3d& de_in) -> Inc {
        Inc r;
        const Stiff k = stiff_at(s_in);   // the stiffness and strength THIS substep starts from
        const Eigen::Matrix3d& De = k.De;
        r.dsig = De * de_in;
        r.tangent = De;
        const double pp = pc_of(ev_in);
        const Eigen::Vector3d sig_tr = s_in + De * de_in;
        const double q_tr = sig_tr.maxCoeff() - sig_tr.minCoeff();
        const bool as0 = fbar(q_tr, k) - gp_in > 1e-12 * (1.0 + std::fabs(gp_in));
        const bool ac0 = cap_on && fcap(sig_tr, pp) > 1e-10 * (1.0 + pp * pp);
        if (!as0 && !ac0) return r;
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
        const bool at_fail = q >= k.qf - 1e-9 * (1.0 + k.qf);
        const double H = at_fail ? 0.0 : 1.0;   // d gamma_p per unit multiplier of a plane
        const double pmean = mean(s_in);
        const double Hcap = cap_on ? p.cap_hardening_modulus(pp) : 0.0;
        r.at_fail = at_fail;

        // Solve the active set for one corner hypothesis; returns false if nothing is active.
        struct Sol {
            Eigen::Vector3d dsig; double dgp = 0.0, dev = 0.0;
            Eigen::Matrix3d tangent; bool as = false, ac = false;
        };
        auto solve = [&](int corner, Sol& sol) {
            bool as = as0, ac = ac0;
            // The planes a corner hypothesis brings: the face (a, c), and at a corner the second
            // face meeting it.
            // At most two planes and the cap: fixed-capacity storage on the stack. This runs
            // several times per substep at every plastic stress point, and heap vectors here cost
            // more than the arithmetic (a Berlin dewatering step: 27 000 plastic points).
            std::array<Plane, 2> planes;
            std::array<double, 2> fp{};   // each plane's yield function where the substep starts
            int nplanes = 0;
            auto add = [&](int i, int j) {
                planes[nplanes] = plane(s_in, i, j, spm, at_fail, k);
                const double qij = s_in(i) - s_in(j);
                // The same failure strength as the plateau test, the drift correction and the
                // Mohr-Coulomb bound use: read at the floored minor stress. Written here with
                // the raw one, a stress point shallower than the floor carried two failure
                // surfaces at once -- this consistency pulled it down to the raw strength while
                // everything else held it at the floored one -- and a strip footing on
                // c = 5 kPa sand stopped converging at 99% of its service load.
                fp[nplanes] = at_fail ? qij - p.q_failure(std::max(s_in(j), plim))
                                      : fbar(qij, k) - gp_in;
                ++nplanes;
            };
            add(o.a, o.c);
            if (corner == kCompression) add(o.a, o.b);
            if (corner == kExtension) add(o.b, o.c);
            std::array<char, 2> on{1, 1};
            const Eigen::Vector3d nc = cap_on ? cap_grad(s_in, o, corner) : Eigen::Vector3d::Zero();
            const double fc_in = cap_on ? fcap(s_in, pp) : 0.0;
            // Koiter: a surface with a negative multiplier is dropped and the rest re-solved.
            for (int pass = 0; pass < 6; ++pass) {
                std::array<Eigen::Vector3d, 3> ns, ms;
                std::array<double, 3> f0s{};
                std::array<int, 3> kind{};   // 0 = shear plane, 1 = cap
                int na = 0;
                if (as)
                    for (int i = 0; i < nplanes; ++i)
                        if (on[i]) { ns[na] = planes[i].n; ms[na] = planes[i].m; f0s[na] = fp[i]; kind[na] = 0; ++na; }
                if (ac) { ns[na] = nc; ms[na] = nc; f0s[na] = fc_in; kind[na] = 1; ++na; }
                if (na == 0) { sol.dsig = De * de_in; sol.tangent = De; sol.as = sol.ac = false; return false; }
                using SmallMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 3, 3>;
                using SmallVec = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3, 1>;
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
                const SmallVec dl = lu.solve(b);
                bool dropped = false;
                int idx = 0;
                if (as)
                    for (int i = 0; i < nplanes; ++i)
                        if (on[i]) { if (dl(idx) < 0.0) { on[i] = 0; dropped = true; } ++idx; }
                if (ac) { if (dl(idx) < 0.0) { ac = false; dropped = true; } ++idx; }
                if (dropped) {
                    bool any = false;
                    for (int i = 0; i < nplanes; ++i) any = any || on[i];
                    if (!any) as = false;
                    continue;
                }
                Eigen::Vector3d ep = Eigen::Vector3d::Zero();
                double dgp = 0.0, dev = 0.0;
                Eigen::Matrix<double, 3, Eigen::Dynamic, 0, 3, 3> N(3, na), M(3, na);
                for (int i = 0; i < na; ++i) {
                    ep += dl(i) * ns[i];
                    N.col(i) = ns[i];
                    M.col(i) = ms[i];
                    if (kind[i] == 0) dgp += H * dl(i);
                    else dev += dl(i) * 2.0 * pmean;
                }
                sol.dsig = De * (de_in - ep);
                sol.dgp = dgp;
                sol.dev = dev;
                sol.tangent = De - De * N * lu.solve(M.transpose() * De);
                sol.as = as;
                sol.ac = ac;
                return true;
            }
            sol.dsig = De * de_in; sol.tangent = De; sol.as = sol.ac = false;
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
        {
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
        r.corner = corner;
        if (corner == kCompression) { r.tie_i = o.b; r.tie_j = o.c; }
        if (corner == kExtension) { r.tie_i = o.a; r.tie_j = o.b; }
        r.plastic = sol.as || sol.ac;
        return r;
    };

    // A substep that STARTS INSIDE the yield surfaces and whose elastic trial ends outside is
    // split where it meets them (Sloan, Abbo & Sheng 2001): the part before the surface is
    // elastic, only the rest flows. Without the split the multiplier is computed from the whole
    // strain increment as if the start were on the surface, so the response jumps -- from
    // elastic to a finite plastic step -- where the trial crosses the surface; a Newton iteration
    // whose points sit there (every point of a column that the cap correction left just inside
    // its cap) then flips points between the two answers and stalls.
    auto increment = [&](const Eigen::Vector3d& s_in, double gp_in, double ev_in,
                         const Eigen::Vector3d& de_in) -> Inc {
        const Stiff k = stiff_at(s_in);
        const double pp = pc_of(ev_in);
        // The larger of the two surfaces, each scaled to its own units: continuous in the stress,
        // and zero where the path first meets either one.
        auto f_at = [&](const Eigen::Vector3d& s) {
            const double qq = s.maxCoeff() - s.minCoeff();
            const double fs = std::max((fbar(std::min(qq, k.qf), k) - gp_in) / (1.0 + std::fabs(gp_in)),
                                       (qq - k.qf) / (1.0 + k.qf));
            const double fc = cap_on ? fcap(s, pp) / (1.0 + pp * pp) : -1.0;
            return std::max(fs, fc);
        };
        const Eigen::Vector3d dse = k.De * de_in;
        const double f0 = f_at(s_in), f1 = f_at(s_in + dse);
        if (!(f0 < 0.0 && f1 > 0.0)) return increment_core(s_in, gp_in, ev_in, de_in);
        // Illinois regula falsi on the elastic fraction.
        double a0 = 0.0, a1 = 1.0, g0 = f0, g1 = f1, alpha = 0.0;
        for (int it = 0; it < 60; ++it) {
            alpha = a1 - g1 * (a1 - a0) / (g1 - g0);
            const double ga = f_at(s_in + alpha * dse);
            if (std::fabs(ga) < 1e-14 || a1 - a0 < 1e-14) break;
            if ((ga > 0.0) == (g1 > 0.0)) { a1 = alpha; g1 = ga; g0 *= 0.5; }
            else { a0 = alpha; g0 = ga; g1 *= 0.5; }
        }
        Inc r = increment_core(s_in + alpha * dse, gp_in, ev_in, (1.0 - alpha) * de_in);
        r.dsig += alpha * dse;
        return r;
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
                             bool as_act, bool ac_act, bool at_fail, int tie_i, int tie_j) {
        for (int it = 0; it < 5; ++it) {
            const Stiff k = stiff_at(s);   // the correction moves sigma3, so it moves these too
            const Order o = order_of(s);
            const int corner = corner_of(o, tie_i, tie_j);
            const double ppc = pc_of(ev_io);
            const double qd = s(o.a) - s(o.c);
            const double fs = at_fail ? (qd - k.qf) : (fbar(qd, k) - gp_io);
            const double fcp = cap_on ? fcap(s, ppc) : -1.0;
            bool corr = false;
            if (as_act && fs > 1e-9 * (1.0 + std::fabs(gp_io) + k.qf)) {
                const double spm = spm_of(qd, k);
                const Plane main = plane(s, o.a, o.c, spm, at_fail, k);
                Eigen::Vector3d n = main.n;
                double h = at_fail ? 0.0 : 1.0;
                if (corner == kCompression) { n += plane(s, o.a, o.b, spm, at_fail, k).n; h *= 2.0; }
                if (corner == kExtension) { n += plane(s, o.b, o.c, spm, at_fail, k).n; h *= 2.0; }
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
            if (ac_act && fcp > 1e-8 * (1.0 + ppc * ppc)) {
                const Eigen::Vector3d g = cap_grad(s, o, corner);
                const Eigen::Vector3d Deg = k.De * g;
                const double pm = mean(s);
                const double Hc_ = p.cap_hardening_modulus(ppc);
                const double dlam_c = fcp / (g.dot(Deg) + 4.0 * ppc * pm * Hc_);
                s -= dlam_c * Deg;
                ev_io += dlam_c * 2.0 * pm;
                corr = true;
            }
            if (!corr) break;
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
            if (over > 1e-12 * (1.0 + k.qf)) {
                s(o.a) -= over;
                if (corner == kExtension) s(o.b) -= over;
            }
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
                          k1.at_fail || k2.at_fail, k2.plastic ? k2.tie_i : k1.tie_i,
                          k2.plastic ? k2.tie_j : k1.tie_j);
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

void hs_calibrate_cap(HardeningSoilParams& p, double K0_NC) {
    p.K0nc = K0_NC;   // the line the oedometer probe starts on
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
