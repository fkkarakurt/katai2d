#pragma once
// SOFT SOIL model — STAGE 1: the material-point core (PLAXIS MMM §10 verbatim; locked
// formulation docs/references/soft-soil-formulation.md). Cam-Clay type: ln-law compression
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

// Manual: p' never drops below unit stress (1 kPa) — a floor against the ln-law/K=p'/κ*
// singularity. ss_step and the FE wrapper (material_model.hpp ss_return_core) must use the
// SAME floor so the elastic-predictor invertibility (exact reconstruction of the trial)
// is not broken.
inline constexpr double kPmin = 1.0;

struct Params {
    double lam_star = 0.1;   // modified compression index λ* [-]
    double kap_star = 0.02;  // modified swelling index κ* [-]
    double nu_ur = 0.15;     // unloading/reloading Poisson ratio
    double c = 0.0;          // effective cohesion [kPa]
    double phi = 0.0;        // effective friction angle [rad] (0 FORBIDDEN — manual §10.3.3)
    double psi = 0.0;        // dilatancy [rad] (SS default 0)
    double K0nc = 0.5;       // normal-consolidation lateral pressure coefficient (→ M)
};

// M(K0NC) — Brinkgreve (1994) / MMM Eq 10-13. Behaviour pin: oedometer primary loading must
// produce σ'_h/σ'_v → K0NC (test_soft_soil (e) measures this DIRECTLY; a transcription error
// in the formula blows up there).
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
// q̃ = σ1 + (δ−1)σ2 − δσ3, sorted σ1 ≥ σ2 ≥ σ3 (compression-positive) — MMM Eq 10-11 (the HS cap measure).
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
inline StepResult ss_substep(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                             const Eigen::Vector3d& deps) {
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi), spsi = std::sin(P.psi);
    const double ccot = P.phi > 1e-12 ? P.c * cphi / sphi : 0.0;
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double M = M_from_K0nc(P);

    // --- Exponentially EXACT elastic mean + deviatoric predictor with G at the trial p --------
    const double p_c = std::max(sig_c.mean(), kPmin);
    const double dev_v = deps.sum();
    const double p_tr = std::max(p_c * std::exp(dev_v / P.kap_star), kPmin);
    const double K_tr = p_tr / P.kap_star;
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
        const double p0 = std::max(s0.mean(), kPmin);
        const double pL = std::max(p0 * std::exp(-L * mv / P.kap_star), kPmin);
        return Eigen::Vector3d(Eigen::Vector3d::Constant(pL) +
                               (s0 - Eigen::Vector3d::Constant(s0.mean())) - L * 2.0 * G * mdev);
    };
    // Scalar Δλ Newton (secant derivative): f(σ(Δλ), pp(Δλ)) = 0. harden=false → pp constant (MC).
    auto scalar_return = [&](const Eigen::Vector3d& s0, double pp0, auto&& ffun,
                             const Eigen::Vector3d& m, bool harden, Eigen::Vector3d& sL,
                             double& ppL) {
        double dl = 0.0;
        auto eval = [&](double L, Eigen::Vector3d& s, double& ppv) {
            s = exp_return(s0, m, L);
            const double devp = L * m.sum();
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
            const double dl = scalar_return(out.sig, out.pp, ffun, m, harden, sL, ppL);
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

// FE initial preconsolidation (K0 seeding; the parallel of hs_initial_pp): pp = f̄(σ0)·OCR_eq
// — f̄ is the equivalent pressure at which the cap passes through σ0 (the counterpart of MMM
// §2.8 p_eq). NC (OCR_eq=1) ⇒ the state sits EXACTLY on the cap (f=0, admissible; a pp=0
// seed would make the cap yield from the start). OCR_eq: OCR mode is the ratio directly;
// POP mode is converted at the caller to the equivalent ratio (σ'_v0+POP)/σ'_v0 (f̄ is NOT
// first-order homogeneous in σ (the c·cotφ shift), but the ratio scaling is the consistent
// counterpart of PLAXIS's vertical-stress-based definition). Floor: pp ≥ max(c·cotφ, unit
// stress) — the threshold ellipse.
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
inline StepResult ss_step(const Params& P, const Eigen::Vector3d& sig_c, double pp_c,
                          const Eigen::Vector3d& deps, int nsub_fixed = 0) {
    const int nsub = nsub_fixed > 0 ? nsub_fixed : ss_nsub(P, deps);
    const Eigen::Vector3d dsub = deps / nsub;
    StepResult out{sig_c, pp_c, false, false, nsub};
    for (int i = 0; i < nsub; ++i) {
        const StepResult r = ss_substep(P, out.sig, out.pp, dsub);
        out.sig = r.sig;
        out.pp = r.pp;
        out.cap_active = out.cap_active || r.cap_active;
        out.mc_active = out.mc_active || r.mc_active;
    }
    return out;
}

}  // namespace katai::core::softsoil
