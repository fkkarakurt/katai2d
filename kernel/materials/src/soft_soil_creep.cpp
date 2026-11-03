#include <katai/materials/soft_soil_creep.hpp>

namespace katai::core::softsoilcreep {

namespace {
// WHERE ON THE ELLIPSE A CREEP STEP ENDS. The creep goes along the normal of the ellipse, m = n / tr(n), and the
// deviatoric part of a normal is along w, the weights of q~, which the order of the principal stresses fixes. So a
// step moves the stress along w, and all the direction asks of the step's end is X = q~ / (M (p' + c cot phi)), the
// position on the ellipse: 1 is its top. With A the trial's q~, B = M (p' + c cot phi) at the step's end and
// C = 4 G L (w . w) / M, backward Euler on the direction reads
//     A - B X = C X / (1 - X^2),   that is   h(X) = (A - B X)(1 - X^2) - C X = 0.
// There is one root below the top, in (0, min(1, A/B)), and one above it, beyond max(1, A/B), each by a monotone
// function, whatever the size of the step.
double position_on_ellipse(double A, double B, double C, bool above_top) {
    const auto h = [&](double X) { return (A - B * X) * (1.0 - X * X) - C * X; };
    double lo, hi;
    if (!above_top) {
        lo = 0.0;
        hi = std::min(1.0, A / B);
        if (!(hi > 0.0)) return 0.0;
        if (C <= 0.0) return hi;
    } else {
        lo = std::max(1.0, A / B);
        if (C <= 0.0) return lo;
        hi = 2.0 * lo;
        for (int i = 0; i < 200 && h(hi) < 0.0; ++i) hi *= 2.0;
    }
    double X = 0.5 * (lo + hi);
    for (int it = 0; it < 100; ++it) {
        const double f = h(X);
        if ((f > 0.0) != above_top) lo = X; else hi = X;   // h falls through its root below the top, rises above it
        const double df = -B * (1.0 - X * X) - 2.0 * X * (A - B * X) - C;
        double next = df != 0.0 ? X - f / df : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (std::fabs(next - X) <= 1e-14 * std::max(X, 1e-300)) return next;
        X = next;
    }
    return X;
}
}  // namespace

StepResult ssc_substep(const Params& P, const softsoil::Params& S,
                       const Eigen::Vector3d& sig_c, double pp_c,
                       const Eigen::Vector3d& deps, double dt) {
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi);
    const double ccot = P.phi > 1e-12 ? P.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = softsoil::M_from_K0nc(S);
    const double beta = (P.lam_star - P.kap_star) / P.mu_star;
    const double ppfloor = std::max(ccot, softsoil::kPmin);

    // --- The same exponentially exact elastic predictor as SS ---------------------------------
    const double p_c = sig_c.mean();
    const double dv = deps.sum();
    const double p_tr = softsoil::mean_after(p_c, dv, P.kap_star, ccot);
    const double K_tr = softsoil::bulk_at(p_tr, P.kap_star, ccot);
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const Eigen::Vector3d de = deps - Eigen::Vector3d::Constant(dv / 3.0);
    const Eigen::Vector3d s_dev_c = sig_c - Eigen::Vector3d::Constant(sig_c.mean());
    const Eigen::Vector3d sig_tr =
        Eigen::Vector3d::Constant(p_tr) + s_dev_c + 2.0 * G * de;

    StepResult out{sig_tr, std::max(pp_c, ppfloor), 0.0, false, 1};

    // --- Creep correction: backward Euler on L = d eps_v^c AND on the direction ------------------
    // The direction used to be the trial's, kept through the step. Near the top of the ellipse that is not a
    // matter of accuracy: tr(n) goes to 0 there, and a step that ends elsewhere than the trial moved by the trial's
    // enormous m. Measured on an undrained triaxial test of a clay with M = 1.296 against a friction line of 1.287
    // (K0nc 0.61): q 121 or 23 kPa for 297 at 200 steps, 1.15e7 kPa from one large step, each reported as solved, and
    // a strip footing stopping at 15 to 36% of its load; the coefficient K0 of an oedometer +2% at the substeps of
    // the rule. Read at the step's end the same test is within 0.01% of an independent integration of the rate
    // equations at every step size, and K0 within 0.03%.
    //
    // A tie in the ordering (Koiter): if the creep step breaks the trial's order, the pair's roles are averaged and
    // the pair equalized. With w2 = w3 the equivalent pressure depends only on s2 + s3, so the equalization changes
    // neither it nor the hardening.
    const double pe_tr = p_eq(S, sig_tr);
    if (pe_tr > 1e-12 && dt > 0.0) {
        int idx[3] = {0, 1, 2};
        std::sort(idx, idx + 3, [&](int a, int b) { return sig_tr(a) > sig_tr(b); });
        // The side of the top the step starts on. The rate is normalised on its volumetric part, so at the top the
        // deviatoric creep grows without bound and drives the stress away from it, from either side: a stress on
        // one side stays there, and under shear a clay that reaches the top from below stays at it.
        const bool above_c = softsoil::detail::q_tilde(sig_c, delta) > M * std::max(p_c + ccot, 1e-9);
        bool tied01 = false, tied12 = false;
        for (int attempt = 0; attempt < 3; ++attempt) {
            Eigen::Vector3d w(1.0, delta - 1.0, -delta);
            if (tied01 && tied12) {
                w.setConstant(0.0);   // full tie: no q~ contribution (purely volumetric)
            } else if (tied01) {
                const double mw = 0.5 * (w(0) + w(1)); w(0) = mw; w(1) = mw;
            } else if (tied12) {
                const double mw = 0.5 * (w(1) + w(2)); w(1) = mw; w(2) = mw;
            }
            const double W = w.squaredNorm();
            const double q0 = w(0) * sig_tr(idx[0]) + w(1) * sig_tr(idx[1]) + w(2) * sig_tr(idx[2]);
            const Eigen::Vector3d s_tr = sig_tr - Eigen::Vector3d::Constant(sig_tr.mean());

            // A creep too small to matter does not hold the stress at the top. To hold it the deviatoric creep would
            // have to be over a million times the volumetric one (|tr n| below 1e-6, the floor this direction has
            // always had), and then the stress crosses with its trial: an overconsolidated clay is elastic there,
            // and passes. Decided once for the step, from the creep the trial asks for.
            bool above = above_c;
            {
                const double pb = std::max(p_tr + ccot, 1e-9);
                const double B = M * pb;
                const double X0 = position_on_ellipse(q0, B, 0.0, above_c);
                const double pe0 = p_tr + X0 * X0 * pb;
                const double L0 = dt * (P.mu_star / P.tau_day) * std::pow(std::max(pe0, 1e-12) / out.pp, beta);
                if (above_c ? q0 < B : q0 > B) {
                    const double X = position_on_ellipse(q0, B, 4.0 * G * L0 * W / M, above_c);
                    if (std::fabs(1.0 - X * X) < 1e-6) above = !above_c;
                }
            }

            auto state_at = [&](double L, Eigen::Vector3d& sL, double& ppL) {
                const double pL = softsoil::mean_after(p_tr, -L, P.kap_star, ccot);
                const double B = M * std::max(pL + ccot, 1e-9);
                const double X = W > 0.0 ? position_on_ellipse(q0, B, 4.0 * G * L * W / M, above) : 0.0;
                const double lam = W > 0.0 ? (q0 - B * X) / W : 0.0;
                sL = Eigen::Vector3d::Constant(pL) + s_tr;
                for (int k = 0; k < 3; ++k) sL(idx[k]) -= lam * w(k);
                // The tied pair takes its common value HERE, so that the creep rate is read at
                // the state that is kept (the Soft Soil lesson, soft_soil.cpp).
                if (tied01 && tied12) {
                    sL.setConstant(sL.mean());
                } else if (tied01) {
                    const double ms = 0.5 * (sL(idx[0]) + sL(idx[1]));
                    sL(idx[0]) = ms; sL(idx[1]) = ms;
                } else if (tied12) {
                    const double ms = 0.5 * (sL(idx[1]) + sL(idx[2]));
                    sL(idx[1]) = ms; sL(idx[2]) = ms;
                }
                ppL = std::max(out.pp * std::exp(L / (P.lam_star - P.kap_star)), ppfloor);
            };
            auto resid = [&](double L, Eigen::Vector3d& sL, double& ppL) {
                state_at(L, sL, ppL);
                const double pe = std::max(p_eq(S, sL), 1e-12);
                return L - dt * (P.mu_star / P.tau_day) * std::pow(pe / ppL, beta);
            };
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

    // --- Mohr-Coulomb check (order: AFTER creep) — single source: softsoil::ss_step
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

StepResult ssc_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                    const Eigen::Vector3d& deps, double dt, int nsub_fixed) {
    const softsoil::Params S = P.ss();
    const int nsub = nsub_fixed > 0 ? nsub_fixed : ssc_nsub(P, S, sig_c, pp_c, deps, dt);
    // Full substeps of the size the rule allows and a last one for the remainder, as in the
    // Soft Soil step (soft_soil.cpp, ss_step): the returned stress is then continuous in the
    // increment where the count changes, instead of jumping by the integration error. The size
    // is read from the same two measures ssc_nsub counts with -- the strain and the creep of
    // the step -- and strain and time are cut alike. At the ceiling of substeps the rule no
    // longer holds the size, and the parts are equal as before.
    double full = nsub > 0 ? 1.0 / nsub : 1.0, last = full;
    if (nsub > 1 && nsub < 500) {
        const double beta = (P.lam_star - P.kap_star) / P.mu_star;
        const double pe = std::max(p_eq(S, sig_c), 1e-12);
        const double ratio = pe / std::max(pp_c, softsoil::kPmin);
        const double Lexp = dt * (P.mu_star / P.tau_day) * std::pow(std::min(ratio, 10.0), beta);
        const double em = deps.sum() / 3.0;
        const double dev = (deps - Eigen::Vector3d::Constant(em)).cwiseAbs().maxCoeff();
        const double measure = std::max(Lexp, std::fabs(deps.sum()) + 2.0 * dev);
        full = std::min(0.05 * P.kap_star / std::max(measure, 1e-300), 1.0 / (nsub - 1));
        last = 1.0 - (nsub - 1) * full;
    }
    StepResult out{sig_c, std::max(pp_c, softsoil::kPmin), 0.0, false, nsub};
    for (int i = 0; i < nsub; ++i) {
        const double part = i + 1 < nsub ? full : last;
        const StepResult r = ssc_substep(P, S, out.sig, out.pp, deps * part, dt * part);
        out.sig = r.sig;
        out.pp = r.pp;
        out.devc += r.devc;
        out.mc_active = out.mc_active || r.mc_active;
    }
    return out;
}

}  // namespace katai::core::softsoilcreep
