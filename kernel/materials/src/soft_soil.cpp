#include <cstdio>
#include <cstdlib>
#include <katai/materials/soft_soil.hpp>

namespace katai::core::softsoil {

// MODIFIED CAM-CLAY, ONE SUBSTEP. The predictor is Soft Soil's -- the exact exponential mean and
// the deviatoric part with G at the trial p. The return is BACKWARD EULER: the plastic strain
// goes along the normal of the ellipse at the RETURNED state, taken back through the ln law for
// its volumetric part, with p_p hardened exponentially. In the von Mises q the normal is radial
// in the deviatoric plane, so the deviatoric stress stays along the trial's and the return is a
// question of two scalars, p and q: s = s_tr (1 - 3 G a L / q_tr), a = df/dq. The normal is found
// by fixed point -- the return along a given normal, the normal at its result, again -- starting
// from the normal at the trial. A return along the trial's normal alone (what Soft Soil does)
// carries a bias that does not vanish on a steady path: on one-dimensional compression the stress
// ratio it settles at depends on the size of the step, measured +2.66% in K0 at the substep size
// the finite elements use. With the normal of the returned state it is -0.02% there (what is left
// is the elastic shear modulus held at the trial's p over the step), and the undrained path,
// which the trial's normal stopped short of the critical state, reaches it.
// The multiplier along a given normal is BRACKETED before it is refined. On the dry side of the
// ellipse the plastic volume change is a dilation, p_p falls with it and the surface shrinks as
// the stress moves; a Newton step from zero can land beyond the root. Where there is no root at
// all -- the softening outruns the elastic unloading, a material instability -- the state is left
// at the closest the search came and flagged.
StepResult cc_substep(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                      const Eigen::Vector3d& deps) {
    const double M = M_of(P), M2 = M * M;
    const double p_c = sig_c.mean(), dev_v = deps.sum();
    const double p_tr = mean_after(p_c, dev_v, P.kap_star);
    const double K_tr = bulk_at(p_tr, P.kap_star);
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const Eigen::Vector3d de = deps - Eigen::Vector3d::Constant(dev_v / 3.0);
    const Eigen::Vector3d s_tr = (sig_c - Eigen::Vector3d::Constant(p_c)) + 2.0 * G * de;
    const double q_tr = std::sqrt(1.5 * s_tr.squaredNorm());
    const double pp0 = std::max(pp_c, kPmin);
    StepResult out{Eigen::Vector3d::Constant(p_tr) + s_tr, pp0, false, false};
    const auto f = [&](double p, double q, double pp) { return q * q / (M2 * std::max(p, 1e-9)) + p - pp; };
    const double f0 = f(p_tr, q_tr, pp0);
    const double ftol = 1e-11 * std::max(1.0, pp0);
    if (f0 <= ftol) return out;

    const double hard = P.lam_star - P.kap_star;
    // The return along a normal (nv = df/dp, the plastic volume per unit multiplier; a = df/dq):
    // the multiplier that puts the state on the hardened ellipse. False where there is none.
    const auto along = [&](double nv, double a, double& p, double& q, double& pp) -> bool {
        const auto state = [&](double L) {
            p = mean_after(p_tr, -L * nv, P.kap_star);
            q = std::fabs(q_tr - 3.0 * G * a * L);
            pp = std::max(pp0 * std::exp(L * nv / hard), kPmin);
            return f(p, q, pp);
        };
        const double dfdl0 = -(K_tr * nv * nv + 3.0 * G * a * a + pp0 * nv / hard);
        const double L_q = a > 0.0 ? q_tr / (3.0 * G * a) : 1e300;   // beyond: through the axis
        double lo = 0.0, flo = f0, hi = -1.0, fhi = 0.0;
        double L = dfdl0 < 0.0 ? std::min(-f0 / dfdl0, L_q) : 1e-6 * std::min(L_q, 1.0);
        double best_L = 0.0, best_f = f0;
        for (int k = 0; k < 200 && hi < 0.0; ++k) {
            const double fl = state(L);
            if (std::fabs(fl) < std::fabs(best_f)) { best_f = fl; best_L = L; }
            if (fl <= 0.0) { hi = L; fhi = fl; break; }
            lo = L; flo = fl;
            if (L >= L_q) break;
            L = std::min(2.0 * L, L_q);
        }
        if (hi < 0.0) { state(best_L); return false; }
        // Illinois: regula falsi that halves the stale end, to |f| <= ftol.
        int side = 0;
        L = hi;
        double fl = fhi;
        for (int it = 0; it < 100 && std::fabs(fl) > ftol && hi - lo > 1e-15 * hi; ++it) {
            L = (lo * fhi - hi * flo) / (fhi - flo);
            fl = state(L);
            if (fl > 0.0) { lo = L; flo = fl; if (side == 1) fhi *= 0.5; side = 1; }
            else { hi = L; fhi = fl; if (side == -1) flo *= 0.5; side = -1; }
        }
        state(L);
        return true;
    };
    const auto normal_at = [&](double p, double q, double& nv, double& a) {
        const double pb = std::max(p, 1e-9);
        nv = 1.0 - q * q / (M2 * pb * pb);
        a = 2.0 * q / (M2 * pb);
    };
    double nv, a, p, q, pp;
    normal_at(p_tr, q_tr, nv, a);
    bool ok = along(nv, a, p, q, pp);
    // The normal of the returned state, by fixed point. It contracts at the rate of the step's
    // size; a fixed point that does not settle keeps the last return that did.
    for (int it = 0; ok && it < 50; ++it) {
        double nv2, a2;
        normal_at(p, q, nv2, a2);
        if (std::fabs(nv2 - nv) <= 1e-13 * (1.0 + std::fabs(nv)) &&
            std::fabs(a2 - a) <= 1e-13 * (1.0 + std::fabs(a)))
            break;
        double p2, q2, pp2;
        if (!along(nv2, a2, p2, q2, pp2)) break;
        nv = nv2; a = a2; p = p2; q = q2; pp = pp2;
    }
    if (!ok) out.unstable = true;
    const double scale = q_tr > 0.0 ? q / q_tr : 0.0;
    out.sig = Eigen::Vector3d::Constant(p) + scale * s_tr;
    out.pp = pp;
    out.cap_active = true;
    return out;
}

StepResult ss_substep(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                      const Eigen::Vector3d& deps) {
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi), spsi = std::sin(P.psi);
    const double ccot = P.phi > 1e-12 ? P.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = M_from_K0nc(P);

    // --- Exponentially EXACT elastic mean + deviatoric predictor with G at the trial p --------
    const double p_c = sig_c.mean();
    const double dev_v = deps.sum();
    const double p_tr = mean_after(p_c, dev_v, P.kap_star, ccot);
    const double K_tr = bulk_at(p_tr, P.kap_star, ccot);
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const Eigen::Vector3d de = deps - Eigen::Vector3d::Constant(dev_v / 3.0);
    const Eigen::Vector3d s_dev_c = sig_c - Eigen::Vector3d::Constant(sig_c.mean());
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p_tr) + s_dev_c + 2.0 * G * de;
    double pp = std::max(pp_c, std::max(ccot, kPmin));

    // Frozen elastic matrix (principal): D = K 1⊗1 + 2G (I − 1⊗1/3).
    auto De_mul = [&](const Eigen::Vector3d& n) {
        const double nv = n.sum();
        return Eigen::Vector3d(Eigen::Vector3d::Constant(K_tr * nv) +
                               2.0 * G * (n - Eigen::Vector3d::Constant(nv / 3.0)));
    };
    auto f_cap = [&](const Eigen::Vector3d& s, double ppv) {
        const double p = s.mean();
        const double qt = detail::q_tilde(s, delta);
        return qt * qt / (M * M * std::max(p + ccot, 1e-9)) + p - ppv;
    };
    auto f_mc = [&](const Eigen::Vector3d& s) {
        Eigen::Vector3d o = s;
        std::sort(o.data(), o.data() + 3, std::greater<double>());
        return (o(0) - o(2)) - (o(0) + o(2)) * sphi - 2.0 * P.c * cphi;
    };

    // --- ANALYTIC surface gradients (sorted-role) + VALIDITY-based edge return ----------------
    // q̃ and MC are defined through the principal ORDERING → cone edges on the σ2=σ3 (TC) /
    // σ1=σ2 (TE) meridians. Single-sided surface flow branches discontinuously at an edge
    // (lateral K0NC split in the oedometer — test (g1) caught it); tolerance-based role
    // averaging also leaves the map DISCONTINUOUS at the tolerance boundary and breaks the
    // global Newton (test (g4) caught it). The correct treatment (Koiter active set; the SS
    // counterpart of the region-consistent logic in our MC model): first the single-sided
    // SURFACE return; if the return breaks the trial ordering, the corresponding EDGE
    // return — the tied pair's roles are averaged as (n_A+n_B)/2 and the pair is equalized
    // to a COMMON value. With w2=w3, f depends only on σ2+σ3, so the equalization (a pure
    // 2G deviatoric transfer) changes NEITHER f NOR the volumetric hardening → this is
    // exactly the two-multiplier edge solution, and the surface/edge composition COINCIDES
    // at the selection boundary (continuous map; derivation in soft-soil-formulation.md).
    // No tolerance → independent of scale and ulp noise.
    auto sort_idx3 = [](const Eigen::Vector3d& s, int idx[3]) {
        idx[0] = 0; idx[1] = 1; idx[2] = 2;
        std::sort(idx, idx + 3, [&](int a, int b) { return s(a) > s(b); });
    };
    auto n_from_w = [](const int idx[3], const Eigen::Vector3d& w) {
        Eigen::Vector3d n;
        for (int k = 0; k < 3; ++k) n(idx[k]) = w(k);
        return n;
    };
    // ∂f_cap/∂σ = a·w + b·1 (sorted-role): a = 2q̃/(M²p̄), b = (1 − q̃²/(M²p̄²))/3.
    auto cap_normal = [&](const Eigen::Vector3d& s, const int idx[3], const Eigen::Vector3d& w) {
        const double pb = std::max(s.mean() + ccot, 1e-9);
        const double qt = s(idx[0]) + (delta - 1.0) * s(idx[1]) - delta * s(idx[2]);
        const double a = 2.0 * qt / (M * M * pb);
        const double b = (1.0 - qt * qt / (M * M * pb * pb)) / 3.0;
        return Eigen::Vector3d(a * n_from_w(idx, w) + Eigen::Vector3d::Constant(b));
    };

    // EXPONENTIAL-MEAN return: the volumetric part of the plastic flow is taken back
    // EXPONENTIALLY through the elastic ln-law (dp = (p/κ*)dε_v^e ⇒
    // p(Δλ) = p₀·exp(−Δλ·n_v/κ*)) — instead of the linear −Δλ·K·n_v. This reproduces the
    // virgin isotropic line EXACTLY in closed form at ANY step size (the
    // Δλ = Δε_v(λ*−κ*)/λ* derivation is in the formulation document) and preserves
    // p-consistency in the MC return too. The deviatoric part is linear with 2G (G at the
    // step's trial p; the remaining step-size error sits only in the deviatoric-volumetric
    // cross terms, which the V&V bands measure).
    auto exp_return = [&](const Eigen::Vector3d& s0, const Eigen::Vector3d& m, double L) {
        const double mv = m.sum();
        const Eigen::Vector3d mdev = m - Eigen::Vector3d::Constant(mv / 3.0);
        const double p0 = s0.mean();
        const double pL = mean_after(p0, -L * mv, P.kap_star, ccot);
        return Eigen::Vector3d(Eigen::Vector3d::Constant(pL) +
                               (s0 - Eigen::Vector3d::Constant(s0.mean())) - L * 2.0 * G * mdev);
    };
    // Scalar Δλ Newton (secant derivative): f(σ(Δλ), pp(Δλ)) = 0. harden=false → pp constant (MC).
    auto scalar_return = [&](const Eigen::Vector3d& s0, double pp0, auto&& ffun,
                             const Eigen::Vector3d& m, bool harden, Eigen::Vector3d& sL,
                             double& ppL, double dl0) {
        double dl = dl0;
        auto eval = [&](double L, Eigen::Vector3d& s, double& ppv) {
            s = exp_return(s0, m, L);
            const double devp = L * m.sum();
            // The hardening law is the model's as written, in p_p itself:
            // p_p = p_p0 exp(-eps_v^p / (lam* - kap*)). The ELASTIC law is in the shifted
            // pressure p' + c cot(phi) (soft_soil.hpp, mean_after); the two are not shifted alike.
            ppv = harden ? std::max(pp0 * std::exp(devp / (P.lam_star - P.kap_star)),
                                    std::max(ccot, kPmin))
                         : pp0;
            return ffun(s, ppv);
        };
        double fL = eval(0.0, sL, ppL);
        // Tolerance scale from the INITIAL violation (|f(0)|): deriving it from pp0 inflated
        // the tolerance to infinity with pp=∞ (the cap-off MC use, soft_soil_creep) and
        // silently skipped the return (measured: MC never capped in the SSC triaxial,
        // q → 9e7).
        const double fscale = std::max(1.0, std::fabs(fL));
        for (int it = 0; it < 60 && std::fabs(fL) > 1e-11 * fscale; ++it) {
            const double h = std::max(1e-13, 1e-7 * (std::fabs(dl) + 1e-7));
            Eigen::Vector3d s2; double pp2;
            const double f2 = eval(dl + h, s2, pp2);
            const double dfdl = (f2 - fL) / h;
            if (std::fabs(dfdl) < 1e-30) break;
            dl -= fL / dfdl;
            if (dl < 0.0) dl = 0.0;
            fL = eval(dl, sL, ppL);
        }
        return dl;
    };

    StepResult out{sig, pp, false, false};

    // Surface→edge cascade (block comment above): attempt 0 is the single-sided surface; if
    // the return breaks the trial ordering, the broken pair(s) are tied and it repeats. A
    // tied pair is equalized to a common value after the return (f/hardening-invariant
    // transfer). w: sorted-role coefficients; is_cap: normal a·w+b·1 (cap) / w (MC flow).
    // Returns true when out.sig (and, if harden, out.pp) was updated.
    auto return_with_edges = [&](auto&& ffun, const Eigen::Vector3d& w, bool harden,
                                 bool is_cap) -> bool {
        int idx[3];
        sort_idx3(out.sig, idx);
        bool tied01 = false, tied12 = false;
        for (int attempt = 0; attempt < 3; ++attempt) {
            Eigen::Vector3d wa = w;
            if (tied01 && tied12) {
                wa.setConstant((w(0) + w(1) + w(2)) / 3.0);
            } else if (tied01) {
                const double mw = 0.5 * (wa(0) + wa(1)); wa(0) = mw; wa(1) = mw;
            } else if (tied12) {
                const double mw = 0.5 * (wa(1) + wa(2)); wa(1) = mw; wa(2) = mw;
            }
            const Eigen::Vector3d m = is_cap ? cap_normal(out.sig, idx, wa) : n_from_w(idx, wa);
            Eigen::Vector3d sL; double ppL;
            // THE SURFACE IS MET BY THE STATE THAT IS KEPT. On an edge the tied pair is brought
            // to a common value after the return, and the multiplier used to be solved with the
            // surface read at the state BEFORE that: q~ weighs the two members of the pair
            // differently, so making them equal moved the stress off the surface it had just
            // been returned to -- measured, 0.11 kPa inside the cap at p' = 40 kPa -- and the
            // update jumped by that much wherever an increment passed from the surface return
            // to the edge return. An equilibrium iteration cannot get under such a jump: it was
            // the out-of-balance force that would not fall in two dimensions, and the reason the
            // numerical tangent never matched the update it was taken from. The surface is now
            // read at the equalized state, so the edge return ends on it and the two returns
            // give one answer where they meet.
            auto tie = [&](Eigen::Vector3d s) {
                if (tied01 && tied12) {
                    s.setConstant(s.mean());
                } else if (tied01) {
                    const double ms = 0.5 * (s(idx[0]) + s(idx[1]));
                    s(idx[0]) = ms; s(idx[1]) = ms;
                } else if (tied12) {
                    const double ms = 0.5 * (s(idx[1]) + s(idx[2]));
                    s(idx[1]) = ms; s(idx[2]) = ms;
                }
                return s;
            };
            auto f_tied = [&](const Eigen::Vector3d& s, double ppv) { return ffun(tie(s), ppv); };
            const double dl = scalar_return(out.sig, out.pp, f_tied, m, harden, sL, ppL, 0.0);
            if (dl <= 0.0) {
                // We were called with f > tol: dl=0 is NOT "inside the surface", it is a
                // secant failure (under single-sided flow the σ2/σ3 crossing kink makes
                // f(Δλ) non-monotone — measured on a large oedometric step: the cap was
                // silently skipped, giving elastic exponential swelling). ADVANCE the
                // cascade: TC edge first (the common case), then the full tie.
                if (!tied12) { tied12 = true; continue; }
                if (!tied01) { tied01 = true; continue; }
                return false;
            }
            if (tied01 && tied12) {
                sL.setConstant(sL.mean());
            } else if (tied01) {
                const double ms = 0.5 * (sL(idx[0]) + sL(idx[1]));
                sL(idx[0]) = ms; sL(idx[1]) = ms;
            } else if (tied12) {
                const double ms = 0.5 * (sL(idx[1]) + sL(idx[2]));
                sL(idx[1]) = ms; sL(idx[2]) = ms;
            }
            // THE CAP'S NORMAL AT THE RETURNED STATE (backward Euler). The return above goes along
            // the normal at the trial, which on a steady path is not the normal of the state the
            // path is on: on one-dimensional compression the stress ratio it settled at depended
            // on the step -- measured +1.2% in K0NC through the finite elements at their substep
            // size, against -0.016% now at the material point's step and -0.0016% at a tenth of
            // it (c = 0, where the formula for M is exact; test_soft_soil (e)). The normal is
            // taken again at the returned state, with the same pairs tied, and the return redone
            // from the trial along it, from the multiplier it ended at; a redone return that
            // fails keeps the last one that did not. The ordering is judged on the state this
            // ends at, so a return whose own normal would cross an edge hands over to the edge
            // return, as the trial's did: where the two meet, the tied state is the untied one,
            // and so is its normal.
            // TWO CORRECTIONS, ALWAYS. The fixed point contracts by about a fifth per pass, so two
            // leave a twenty-fifth of the trial normal's error -- the bias above falls from +1.2%
            // to below the step's own -- and a fixed count keeps the update a smooth function of
            // the strain increment, which the numerical tangent differentiates. Iterated to 1e-13
            // it took sixteen passes on average and made a consolidation of this soil 7.6 times
            // slower, for a K0NC that moved in the fifth digit.
            if (is_cap) {
                double dl_k = dl;
                for (int fp = 0; fp < 2; ++fp) {
                    const Eigen::Vector3d m2 = cap_normal(sL, idx, wa);
                    Eigen::Vector3d s2; double pp2;
                    const double dl2 = scalar_return(out.sig, out.pp, f_tied, m2, harden, s2, pp2, dl_k);
                    if (dl2 <= 0.0) break;
                    dl_k = dl2;
                    if (tied01 && tied12) {
                        s2.setConstant(s2.mean());
                    } else if (tied01) {
                        const double ms = 0.5 * (s2(idx[0]) + s2(idx[1]));
                        s2(idx[0]) = ms; s2(idx[1]) = ms;
                    } else if (tied12) {
                        const double ms = 0.5 * (s2(idx[1]) + s2(idx[2]));
                        s2(idx[1]) = ms; s2(idx[2]) = ms;
                    }
                    sL = s2; ppL = pp2;
                }
            }
            const double slack = 1e-9 * (1.0 + std::fabs(sL(idx[0])));
            const bool bad01 = !tied01 && sL(idx[0]) - sL(idx[1]) < -slack;
            const bool bad12 = !tied12 && sL(idx[1]) - sL(idx[2]) < -slack;
            if (!bad01 && !bad12) {
                out.sig = sL;
                if (harden) out.pp = ppL;
                return true;
            }
            tied01 = tied01 || bad01;
            tied12 = tied12 || bad12;
        }
        return false;
    };

    for (int pass = 0; pass < 6; ++pass) {
        bool changed = false;
        // --- Mohr-Coulomb return (flat surface; non-associated ψ; pp constant) -----------------
        if (f_mc(out.sig) > 1e-9) {
            if (return_with_edges([&](const Eigen::Vector3d& s, double) { return f_mc(s); },
                                  Eigen::Vector3d(1.0 - spsi, 0.0, -(1.0 + spsi)), false, false)) {
                out.mc_active = true;
                changed = true;
            }
        }
        // --- Cap return (associated + exponential p_p hardening) -------------------------------
        if (f_cap(out.sig, out.pp) > 1e-9) {
            if (return_with_edges(
                    [&](const Eigen::Vector3d& s, double ppv) { return f_cap(s, ppv); },
                    Eigen::Vector3d(1.0, delta - 1.0, -delta), true, true)) {
                out.cap_active = true;
                changed = true;
            }
        }
        if (!changed) break;
    }
    return out;
}

StepResult ss_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                   const Eigen::Vector3d& deps, int nsub_fixed) {
    // THE SUBSTEPS HAVE A FIXED SIZE AND THE LAST ONE TAKES WHAT IS LEFT. The increment used to
    // be cut into nsub EQUAL parts, with nsub an integer of the increment's size -- so an
    // increment a hair larger than the last one was integrated with one substep more, and the
    // stress it returned jumped by the difference in integration error. Inside an equilibrium
    // iteration, where each iterate changes the increment a little, that jump is noise in the
    // internal force that no iteration can get under: measured on a footing pushed into this
    // soil, the out-of-balance force stopped at 3-4e-3 of the forces in play, in pairs of
    // opposite sign on neighbouring nodes all over the mesh, and the phase was abandoned at
    // 70% of its settlement at a tolerance of 3e-3 (and took 1 300 iterations at 1e-2).
    // With full substeps of the size the rule allows and a last one for the remainder, the
    // returned stress is a continuous function of the increment: where the count goes up by
    // one the new last substep has no length.
    // A pinned count (the columns of the numerical tangent) keeps the count of the run it
    // follows and the same rule for the sizes, so the two runs differ by the perturbation only.
    const double em = deps.sum() / 3.0;
    const double dev = (deps - Eigen::Vector3d::Constant(em)).cwiseAbs().maxCoeff();
    const double measure = std::fabs(deps.sum()) + 2.0 * dev;
    const double h0 = 0.05 * P.kap_star;   // the size ss_nsub allows
    const int nsub = nsub_fixed > 0 ? nsub_fixed : ss_nsub(P, deps);
    const double full = nsub > 1 ? std::min(h0 / std::max(measure, 1e-300), 1.0 / (nsub - 1)) : 0.0;
    const double last = 1.0 - (nsub - 1) * full;
    StepResult out{sig_c, pp_c, false, false, nsub};
    for (int i = 0; i < nsub; ++i) {
        const Eigen::Vector3d dsub = deps * (i + 1 < nsub ? full : last);
        const StepResult r = P.cam_clay ? cc_substep(P, out.sig, out.pp, dsub)
                                        : ss_substep(P, out.sig, out.pp, dsub);
        out.sig = r.sig;
        out.pp = r.pp;
        out.cap_active = out.cap_active || r.cap_active;
        out.mc_active = out.mc_active || r.mc_active;
        out.unstable = out.unstable || r.unstable;
    }
    return out;
}

}  // namespace katai::core::softsoil
