// HSsmall -- small-strain stiffness law verification (Faz A.5). The Hardening Soil model with
// small-strain stiffness adds two parameters (G0, gamma0.7) and a Hardin-Drnevich degradation:
//   G_s/G0 = 1/(1+0.385 |g|/g07)          (Eq 7-3; G_s = 0.722 G0 at g=g07)
//   tau    = G_s g = G0 g/(1+0.385 g/g07) (Eq 7-7)
//   G_t    = dtau/dg = G0/(1+0.385 g/g07)^2, with cut-off G_t >= G_ur (Eq 7-8/7-9)
//   g_cut  = (1/0.385)(sqrt(G0/G_ur)-1) g07  (Eq 7-10)
// plus the stress-level power law (same as HS). All checked against the closed form (round-off).
// (See docs/references/hssmall-formulation.md; PLAXIS 2D Material Models Manual ch.7.)
#include <katai/materials/hardening_soil.hpp>
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

    // gamma_hist accumulates -> degradation. With large committed gamma_hist (>> cut-off) HSsmall == HS.
    GaussState comm_deg = comm; comm_deg.gamma_hist = 0.05;   // >> gamma_cutoff
    GaussState t_deg;
    katai::core::hs_forward(hss, comm_deg, de, t_deg);
    const double k_deg = std::fabs(t_deg.stress(1) - comm.stress(1)) / std::fabs(de(1));
    std::printf("  large gamma_hist: HSsmall stiffness=%.3e  HS=%.3e  (degraded to HS)\n", k_deg, k_hs);
    check(close(k_deg, k_hs, 1e-9), "HSsmall degrades exactly to HS beyond the cut-off strain");

    // gamma_hist accumulates monotonically through hs_forward.
    check(t_hss.gamma_hist > 0.0, "gamma_hist accumulates under deviatoric straining");
    check(t_hs.gamma_hist == 0.0, "plain HS carries gamma_hist = 0 (HSsmall off)");
}

} // namespace

int main() {
    test_hssmall_law();
    test_hssmall_fe();
    if (g_failures == 0) {
        std::printf("OK: HSsmall small-strain stiffness law verified (Hardin-Drnevich, cut-off, stress law)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
