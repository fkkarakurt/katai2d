// HSsmall -- small-strain stiffness law verification (Faz A.5). The Hardening Soil model with
// small-strain stiffness adds two parameters (G0, gamma0.7) and a Hardin-Drnevich degradation:
//   G_s/G0 = 1/(1+0.385 |g|/g07)          (Eq 7-3; G_s = 0.722 G0 at g=g07)
//   tau    = G_s g = G0 g/(1+0.385 g/g07) (Eq 7-7)
//   G_t    = dtau/dg = G0/(1+0.385 g/g07)^2, with cut-off G_t >= G_ur (Eq 7-8/7-9)
//   g_cut  = (1/0.385)(sqrt(G0/G_ur)-1) g07  (Eq 7-10)
// plus the stress-level power law (same as HS). All checked against the closed form (round-off).
// (See docs/references/hssmall-formulation.md; PLAXIS 2D Material Models Manual ch.7.)
#include <katai/materials/hardening_soil.hpp>
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/material_model.hpp>

#include <cmath>
#include <cstdio>

#include <Eigen/Dense>

using katai::core::GaussState;
using katai::core::HardeningSoilParams;
using katai::core::MaterialModel;
using katai::core::MaterialType;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(b)); }

void test_hssmall_law() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.cohesion = 0.0; p.friction = 35.0 * std::acos(-1.0) / 180.0; p.nu_ur = 0.2;
    p.G0_ref = 1.2e5; p.gamma07 = 1.5e-4;   // small-strain params
    const double pr = p.p_ref;               // at reference stress -> stiffness_factor = 1

    // (1) Eq 7-3: G_s(g07) = 0.722 G0.
    const double Gs07 = p.g_secant(p.gamma07, pr);
    std::printf("  G_s(g07)/G0 = %.4f (Hardin-Drnevich 0.722)\n", Gs07 / p.G0(pr));
    check(close(Gs07 / p.G0(pr), 0.722, 2e-3), "G_s = 0.722 G0 at gamma = gamma0.7");

    // (2) Eq 7-7 hyperbola + (3) Eq 7-8 tangent = d(tau)/d(gamma) (uncapped, small strain).
    auto tau = [&](double g) { return p.g_secant(g, pr) * g; };
    double max_t = 0.0;
    for (double g = 1e-6; g < 5e-4; g *= 1.7) {
        const double gt_fd = (tau(g + 1e-9) - tau(g - 1e-9)) / 2e-9;
        const double gt = p.G0(pr) / std::pow(1.0 + 0.385 * g / p.gamma07, 2.0);  // uncapped Eq 7-8
        max_t = std::fmax(max_t, std::fabs(gt_fd - gt) / gt);
    }
    std::printf("  tangent G_t = d(tau)/d(gamma): max rel err = %.2e\n", max_t);
    check(max_t < 1e-4, "G_t = d(tau)/d(gamma) = G0/(1+0.385 g/g07)^2 (Eq 7-8)");

    // (4) Eq 7-10 cut-off: at gamma_cutoff the uncapped tangent equals G_ur; g_tangent clamps to G_ur beyond.
    const double gc = p.gamma_cutoff(pr);
    const double gt_uncapped_at_gc = p.G0(pr) / std::pow(1.0 + 0.385 * gc / p.gamma07, 2.0);
    std::printf("  gamma_cutoff=%.3e: uncapped G_t=%.1f vs G_ur=%.1f\n", gc, gt_uncapped_at_gc, p.Gur(pr));
    check(close(gt_uncapped_at_gc, p.Gur(pr), 1e-6), "uncapped G_t(gamma_cutoff) = G_ur (Eq 7-10)");
    check(close(p.g_tangent(2.0 * gc, pr), p.Gur(pr), 1e-12), "G_t clamps to G_ur beyond cut-off");
    check(p.G0(pr) > p.Gur(pr), "G0 > G_ur (small-strain stiffer than unload/reload)");

    // (5) Stress-level power law: G0 scales like the HS moduli (factor (sigma/p_ref)^m for c=0).
    const double s3 = 400.0;
    const double factor = std::pow(s3 / pr, p.m);   // c=0 -> (sigma3 sinphi)/(p_ref sinphi) = sigma3/p_ref
    std::printf("  stress dependence: G0(400)/G0(100) = %.4f (expected %.4f)\n", p.G0(s3) / p.G0(pr), factor);
    check(close(p.G0(s3) / p.G0(pr), factor, 1e-12), "G0 follows the HS stress power law");

    // (6) Et_small reduces to E_ur at large strain (HSsmall -> HS).
    check(close(p.Et_small(1.0, pr), p.Eur(pr), 1e-9), "Et_small -> E_ur at large strain (HSsmall -> HS)");
    const double E0 = 2.0 * (1.0 + p.nu_ur) * p.G0(pr);
    check(close(p.Et_small(1e-10, pr), E0, 1e-5), "Et_small -> E0 = 2(1+nu)G0 at very small strain");
    check(p.Et_small(1e-10, pr) > p.Eur(pr), "small-strain modulus stiffer than unload/reload E_ur");
}

// FE integration (material-point level): HSsmall makes the elastic (unload/reload) stiffness
// strain-dependent. At small accumulated strain it is governed by E0 (>> E_ur) -> stiffer; once
// gamma_hist exceeds the cut-off the response degrades exactly to HS. Verified via hs_forward.
void test_hssmall_fe() {
    HardeningSoilParams base;
    base.E50_ref = 3.0e4; base.Eur_ref = 9.0e4; base.Eoed_ref = 3.0e4; base.m = 0.5; base.p_ref = 100.0;
    base.cohesion = 0.0; base.friction = 35.0 * std::acos(-1.0) / 180.0; base.nu_ur = 0.2;
    base.cap_beta = 0.0;   // cap off (shear-only HS) for a clean elastic-stiffness comparison

    MaterialModel hs{MaterialType::HardeningSoil, 0, 0}; hs.hs = base;            // plain HS (G0_ref=0)
    MaterialModel hss = hs; hss.hs.G0_ref = 2.7e5; hss.hs.gamma07 = 1.5e-4;       // HSsmall

    // Committed: isotropic confining 100 kPa (compression-negative).
    GaussState comm; comm.stress << -100.0, -100.0, 0.0; comm.stress_zz = -100.0;
    const Eigen::Vector3d de(0.0, -2.0e-6, 0.0);   // tiny confined vertical compression

    GaussState t_hs, t_hss;
    katai::core::hs_forward(hs, comm, de, t_hs);
    katai::core::hs_forward(hss, comm, de, t_hss);   // gamma_hist=0 -> uses E0
    const double k_hs = std::fabs(t_hs.stress(1) - comm.stress(1)) / std::fabs(de(1));
    const double k_hss = std::fabs(t_hss.stress(1) - comm.stress(1)) / std::fabs(de(1));
    std::printf("  small strain: HSsmall stiffness=%.3e  HS=%.3e  ratio=%.2f\n", k_hss, k_hs, k_hss / k_hs);
    check(k_hss > 2.0 * k_hs, "HSsmall is much stiffer than HS at very small strain (E0 vs E_ur)");

    // gamma_hist accumulates -> degradation. Beyond the cut-off the OVERLAY is spent and the
    // unload/reload stiffness is E_ur again, bit for bit.
    check(katai::core::hs_small_strain_params(hss.hs, 0.05).Eur_ref == hss.hs.Eur_ref,
          "beyond the cut-off the small-strain overlay returns E_ur exactly");

    // ...and at the material point the two are then the same material -- ON AN ELASTIC STEP.
    // The probe has to be elastic to make that claim, because what degrades to HS beyond the
    // cut-off is the STIFFNESS, not the whole model: below phi_cv the two keep a different FLOW
    // RULE for good (sec 7.9.1, test_li_dafalias_fe). This check used to run a plastic step and
    // assert equality there, which was true only while HSsmall had no sec 7.9.1 branch -- a
    // claim wider than the thing it measured. Committed gamma_p is put beyond reach so the
    // increment stays inside the yield surface and only the stiffness is being compared.
    GaussState comm_deg = comm; comm_deg.gamma_hist = 0.05; comm_deg.gamma_p = 1.0;
    GaussState t_deg, t_hs_el;
    katai::core::hs_forward(hss, comm_deg, de, t_deg);
    katai::core::hs_forward(hs, comm_deg, de, t_hs_el);
    const double k_deg = std::fabs(t_deg.stress(1) - comm.stress(1)) / std::fabs(de(1));
    const double k_hs_el = std::fabs(t_hs_el.stress(1) - comm.stress(1)) / std::fabs(de(1));
    std::printf("  large gamma_hist (elastic step): HSsmall stiffness=%.6e  HS=%.6e\n", k_deg, k_hs_el);
    check(close(k_deg, k_hs_el, 1e-12), "HSsmall degrades exactly to HS beyond the cut-off strain");

    // gamma_hist accumulates monotonically through hs_forward.
    check(t_hss.gamma_hist > 0.0, "gamma_hist accumulates under deviatoric straining");
    check(t_hs.gamma_hist == 0.0, "plain HS carries gamma_hist = 0 (HSsmall off)");
}

// --- sec. 7.9.1: the mobilised dilatancy BELOW the phase-transformation line ----------------
// Rowe's formula returns a NEGATIVE psi_m wherever phi_m < phi_cv. Plain HS takes zero there.
// HSsmall does not: "bounding the lower value of psi_m may sometimes yield too little plastic
// volumetric strains", so below phi_cv it puts a small CONTRACTION there instead, after
// Li & Dafalias (2000) -- MMM Eq 7-19..7-23.
//
// The oracle below is written straight from those five equations and shares no line of code
// with the kernel's HsDilatancy, so agreement between them is a real comparison and not a
// tautology. The case is the one the manual itself plots in Figure 7-10: phi=35, psi=5.
void test_li_dafalias_dilatancy() {
    const double deg = std::acos(-1.0) / 180.0;
    const double phi = 35.0 * deg, psi = 5.0 * deg;

    const double sphi = std::sin(phi), sps = std::sin(psi);
    const double scv = (sphi - sps) / (1.0 - sphi * sps);  // Rowe's critical state
    const double Mc = 6.0 * scv / (3.0 - scv);             // Eq 7-20
    const double floor_s = sphi / (2.0 - sphi);            // Eq 7-23
    auto oracle = [&](double sphi_m) {
        const double s = std::fmax(sphi_m, floor_s);                              // Eq 7-23
        const double Md = 6.0 * s / (3.0 - s);                                    // Eq 7-21
        const double qqa = std::fmax((1.0 - scv) / scv * (s / (1.0 - s)), 1e-4);  // Eq 7-22
        return 0.1 * (-Mc * std::exp(std::log((Mc / Md) * qqa) / 15.0) + Md);     // Eq 7-19
    };
    auto rowe = [&](double s) { return (s - scv) / (1.0 - s * scv); };

    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.cohesion = 0.0; p.friction = phi; p.dilatancy = psi; p.nu_ur = 0.2;
    p.G0_ref = 1.2e5; p.gamma07 = 1.5e-4;         // HSsmall
    HardeningSoilParams ph = p; ph.G0_ref = 0.0;  // the SAME soil as plain HS

    const auto dil = katai::core::detail::hs_dilatancy(p);
    const auto hs = katai::core::detail::hs_dilatancy(ph);

    std::printf("  phi_cv = %.6f deg   Eq 7-23 floor = %.6f deg\n",
                std::asin(scv) / deg, std::asin(floor_s) / deg);

    // (1) THE TWO BRANCHES MEET EXACTLY. At phi_m = phi_cv we have M_d = M_c and q/q_a = 1, so
    //     the logarithm vanishes and Eq 7-19 returns zero -- which is where Rowe changes sign.
    //     The composite rule is continuous by construction; there is no step to tolerate. In
    //     double precision Eq 7-19 lands on 0 exactly, so this is pinned as an equality.
    check(oracle(scv) == 0.0, "Eq 7-19 vanishes EXACTLY at phi_m = phi_cv (M_d = M_c, q/q_a = 1)");
    check(dil(scv) == 0.0, "psi_m = 0 at phi_m = phi_cv");
    const double eps = 1e-9;
    std::printf("  continuity across phi_cv: psi_m(-eps) = %+.3e   psi_m(+eps) = %+.3e\n",
                dil(scv - eps), dil(scv + eps));
    check(std::fabs(dil(scv - eps)) < 1e-8, "no jump entering the Li & Dafalias branch");
    check(dil(scv + eps) > 0.0, "Rowe (dilatant) on the other side of phi_cv");

    // (2) The kernel reproduces Eq 7-19..7-23 across the whole branch -- and plain HS, the same
    //     soil with G0_ref = 0, keeps the zero cut-off at every one of those points.
    double worst = 0.0;
    bool hs_zero = true, mono = true, bounded = true;
    const double plateau = dil(floor_s);
    double prev = -1.0;
    for (double pm = 0.5; pm < std::asin(scv) / deg; pm += 0.05) {
        const double s = std::sin(pm * deg);
        worst = std::fmax(worst, std::fabs(dil(s) - oracle(s)));
        if (hs(s) != 0.0) hs_zero = false;
        if (dil(s) < prev - 1e-16) mono = false;   // psi_m rises monotonically towards phi_cv
        if (dil(s) < plateau - 1e-16) bounded = false;
        prev = dil(s);
    }
    std::printf("  vs Eq 7-19..7-23 oracle: worst |diff| = %.2e over the branch\n", worst);
    check(worst < 1e-15, "kernel reproduces Eq 7-19..7-23");
    check(hs_zero, "plain HS keeps psi_m = 0 below phi_cv (the sec 7.9.1 branch is HSsmall-only)");
    check(mono, "psi_m rises monotonically towards phi_cv");
    check(bounded, "Eq 7-23 bounds the contraction: no psi_m below the plateau");

    // (3) The plateau: Eq 7-23 floors phi_m, so psi_m is CONSTANT below 23.71 deg -- and that
    //     constant is the most contractant state the model can reach. Rowe, unbounded, would
    //     have asked for -7.96 deg there and -19.90 deg at phi_m = 12 deg; HS discards both.
    std::printf("  plateau psi_m = %.6f deg (Rowe there: %.4f deg; HS: 0)\n",
                std::asin(plateau) / deg, std::asin(rowe(floor_s)) / deg);
    check(close(std::asin(plateau) / deg, -1.678935, 1e-6), "plateau psi_m = -1.678935 deg");
    check(dil(std::sin(5.0 * deg)) == plateau && dil(std::sin(20.0 * deg)) == plateau,
          "psi_m is bit-for-bit constant below the Eq 7-23 floor");

    // (4) THE DILATANCY CUT-OFF STILL MEANS WHAT IT SAYS. Once the soil has dilated to e_max the
    //     cut-off sets the MOBILISED angle back to zero. The engine signals that by zeroing psi,
    //     which for plain HS is the same statement -- but sec. 7.9.1 reads psi only through
    //     phi_cv, so a zeroed psi alone would move phi_cv up to phi and switch the Li & Dafalias
    //     contraction ON across the whole pre-failure range: a "stop dilating" option that starts
    //     producing volume loss. HardeningSoilParams::dilatancy_cut is what prevents that.
    HardeningSoilParams pcut = p; pcut.dilatancy = 0.0; pcut.dilatancy_cut = true;
    HardeningSoilParams pzero = p; pzero.dilatancy = 0.0;  // psi zeroed, cut-off NOT declared
    const auto cut = katai::core::detail::hs_dilatancy(pcut);
    const auto zeroed = katai::core::detail::hs_dilatancy(pzero);
    bool cut_is_zero = true, zeroed_would_contract = false;
    for (double pm = 1.0; pm < 34.9; pm += 0.05) {
        const double s = std::sin(pm * deg);
        if (cut(s) != 0.0) cut_is_zero = false;
        if (zeroed(s) < 0.0) zeroed_would_contract = true;
    }
    std::printf("  dilatancy cut-off engaged: psi_m == 0 everywhere (%s); psi=0 alone would have "
                "contracted (%s)\n", cut_is_zero ? "yes" : "NO", zeroed_would_contract ? "yes" : "no");
    check(cut_is_zero, "with the dilatancy cut-off engaged psi_m is zero, not Li & Dafalias");
    check(zeroed_would_contract, "...and that is not vacuous: psi = 0 alone does reach the branch");

    // (5) SENTINEL. Above phi_cv the two models must be the same rule (Rowe), bit for bit: this
    //     change may only touch the region where Rowe was being discarded.
    bool same_above = true;
    for (double pm = std::asin(scv) / deg; pm <= 35.0; pm += 0.05) {
        const double s = std::sin(pm * deg);
        if (dil(s) != hs(s)) same_above = false;
    }
    check(same_above, "above phi_cv HSsmall and HS are bit-for-bit the same rule");
}

// The same thing at the material point, where it changes an answer. ISOLATION: at
// gamma_hist >> gamma_cutoff the small-strain overlay has degraded to E_ur exactly (pinned in
// test_hssmall_fe), so HSsmall and plain HS then share every stiffness in the model. Under a
// purely deviatoric strain increment the ONLY thing left that can separate them is the flow
// rule -- so whatever mean stress difference appears here IS sec 7.9.1 and nothing else.
void test_li_dafalias_fe() {
    const double deg = std::acos(-1.0) / 180.0;
    HardeningSoilParams base;
    base.E50_ref = 3.0e4; base.Eur_ref = 9.0e4; base.Eoed_ref = 3.0e4; base.m = 0.5;
    base.p_ref = 100.0; base.cohesion = 0.0; base.friction = 35.0 * deg;
    base.dilatancy = 5.0 * deg; base.nu_ur = 0.2;
    base.cap_beta = 0.0;  // cap off: the cap contracts too, and this probe is about the shear flow

    MaterialModel hs{MaterialType::HardeningSoil, 0, 0}; hs.hs = base;
    MaterialModel hss = hs; hss.hs.G0_ref = 2.7e5; hss.hs.gamma07 = 1.5e-4;

    GaussState comm; comm.stress << -100.0, -100.0, 0.0; comm.stress_zz = -100.0;
    comm.gamma_hist = 0.05;                        // >> cut-off: the overlay is spent
    const Eigen::Vector3d de(-3.0e-4, 3.0e-4, 0.0);  // pure deviatoric: eps_v = 0 exactly

    GaussState t_hs, t_hss;
    katai::core::hs_forward(hs, comm, de, t_hs);
    katai::core::hs_forward(hss, comm, de, t_hss);
    auto mean = [](const GaussState& g) { return (g.stress(0) + g.stress(1) + g.stress_zz) / 3.0; };

    // Where on the curve are we? The probe is only meaningful below phi_cv (= 30.798 deg).
    const double q = std::fabs(t_hss.stress(1) - t_hss.stress(0));
    const double sphi_m = q / std::fabs(t_hss.stress(1) + t_hss.stress(0));
    std::printf("  mobilised phi_m = %.3f deg (phi_cv = 30.798 deg -> Li & Dafalias branch)\n",
                std::asin(sphi_m) / deg);
    check(std::asin(sphi_m) / deg < 30.798, "the probe sits below phi_cv");

    // Plain HS: psi_m = 0 -> m_g = (1,-1/2,-1/2), trace 0 -> the plastic flow preserves volume.
    // The total increment is deviatoric, so the elastic volumetric strain stays zero too and
    // the mean stress cannot move off -100. That closed form is the baseline this is measured
    // against -- not the other model's number.
    std::printf("  mean stress: HS %.9f (closed form -100)   HSsmall %.9f   diff %+.3e kPa\n",
                mean(t_hs), mean(t_hss), mean(t_hss) - mean(t_hs));
    check(close(mean(t_hs), -100.0, 1e-9), "HS: volume-preserving flow leaves the mean stress at -100");
    check(mean(t_hss) > mean(t_hs) + 1e-6,
          "HSsmall CONTRACTS below phi_cv: plastic contraction unloads the mean effective stress");

    // Direction check. Contraction (comp-positive eps_v^p > 0) under a constant total volume
    // means the ELASTIC volumetric strain goes the other way, so the soil sheds mean effective
    // stress -- in undrained terms, it builds pore pressure. That is the realistic direction and
    // it is why the manual bothers: HS's zero simply produced no volumetric strain at all here.
    check(mean(t_hss) < 0.0, "still in compression (a small correction, not a sign flip)");
}

} // namespace

int main() {
    test_hssmall_law();
    test_hssmall_fe();
    test_li_dafalias_dilatancy();
    test_li_dafalias_fe();
    if (g_failures == 0) {
        std::printf("OK: HSsmall verified (Hardin-Drnevich, cut-off, stress law, Li & Dafalias sec 7.9.1)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
