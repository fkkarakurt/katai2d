// Hardening Soil -- robust multi-surface integrator (hs_integrate: explicit substepping
// + Koiter active-set, Sloan/Potts&Gens). Replaces the fragile implicit both-active
// solver. Verified to (1) reproduce the shear hyperbola (cap off) in a drained triaxial,
// (2) give the cap K_iso in isotropic compression, and (3) integrate the Berlin Sand III
// oedometer (E=105 MPa) ROBUSTLY -- the high-stiffness case where the previous implicit
// solver diverged. (See docs/references/hardening-soil-formulation.md.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::HsIntegrated;
using katai::core::hs_integrate;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

// One drained-triaxial step (sigma3 const): solve lateral strain via the tangent.
void triaxial_step(const HardeningSoilParams& p, double sigma3, double de1,
                   Eigen::Vector3d& sig, double& gp, double& pp, double& eps1) {
    double dlat = 0.0; HsIntegrated r;
    for (int it = 0; it < 25; ++it) {
        r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, dlat, dlat));
        const double g = r.stress(2) - sigma3;
        if (std::fabs(g) <= 1e-9 * (1.0 + std::fabs(sigma3))) break;
        const double dgdl = r.tangent(2, 1) + r.tangent(2, 2);
        dlat -= g / dgdl;
    }
    sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += de1;
}

void test_triaxial_hyperbola() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.dilatancy = 0.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.cap_beta = 0.0;  // cap off -> pure shear
    const double sigma3 = 100.0;
    const double Ei = p.Ei(sigma3), qa = p.q_asymptote(sigma3), qf = p.q_failure(sigma3);
    auto hyper = [&](double q) { return (1.0 / Ei) * q / (1.0 - q / qa); };

    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0, pp = 0, eps1 = 0, maxrel = 0; int n = 0; double frac = 0.25;
    for (int s = 0; s < 6000; ++s) {
        triaxial_step(p, sigma3, 5e-6, sig, gp, pp, eps1);
        const double q = sig(0) - sig(2);
        if (q >= frac * qf && frac < 0.76) {
            maxrel = std::fmax(maxrel, std::fabs(eps1 - hyper(q)) / hyper(q));
            ++n; frac += 0.25;
        }
    }
    std::printf("  triaxial (substep): checks=%d  max rel err vs hyperbola=%.3e\n", n, maxrel);
    check(n == 3 && maxrel < 1.5e-3, "robust integrator reproduces the shear hyperbola");
}

void test_isotropic_cap() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 30.0 * kPi / 180.0; p.nu_ur = 0.2;
    p.cap_alpha = 0.6; p.cap_beta = 2.0e-3;
    const double p0 = 100.0, e = 1e-6;
    Eigen::Vector3d sig(p0, p0, p0);
    const HsIntegrated r = hs_integrate(p, sig, 0.0, p0, Eigen::Vector3d(e, e, e));
    const double Ke = p.bulk(p0), Hcap = p.cap_hardening_modulus(r.pp);
    const double K_iso = Ke * Hcap / (Ke + Hcap);
    const double dp = (r.stress(0) + r.stress(1) + r.stress(2)) / 3.0 - p0;
    std::printf("  isotropic (substep): K_iso=%.1f  return secant=%.1f\n", K_iso, dp / (3 * e));
    check(close(dp / (3 * e), K_iso, 5e-3), "robust integrator gives cap K_iso");
}

void test_berlin_sand_oedometer_robust() {
    // Berlin Sand III (Rocscience/PLAXIS Table 15.1): the HIGH-stiffness case that broke
    // the previous implicit both-active solver. Here it must integrate robustly.
    HardeningSoilParams p;
    p.p_ref = 100.0; p.E50_ref = 105.0e3; p.Eur_ref = 315.0e3; p.Eoed_ref = 105.0e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38.0 * kPi / 180.0;
    p.dilatancy = 6.0 * kPi / 180.0; p.cohesion = 1.0; p.Rf = 0.9;
    p.cap_alpha = 1.0; p.cap_beta = 1.0e-3;  // any reasonable cap

    Eigen::Vector3d sig(5.0, 5.0, 5.0);
    double gp = 0, pp = 5.0, eps1 = 0;
    bool finite = true, monotone = true; double prev_gp = 0, K0 = 0, Eoed = 0, s1b = 5.0;
    for (int i = 0; i < 2000; ++i) {
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(5e-5, 0, 0));
        if (!r.stress.allFinite()) { finite = false; break; }
        if (r.gamma_p < prev_gp - 1e-10) monotone = false;
        prev_gp = r.gamma_p;
        if (s1b < p.p_ref && r.stress(0) >= p.p_ref) {
            Eoed = (r.stress(0) - s1b) / 5e-5; K0 = r.stress(2) / r.stress(0);
        }
        s1b = r.stress(0);
        sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += 5e-5;
        if (sig(0) > 3.0 * p.p_ref) break;
    }
    std::printf("  Berlin Sand oedometer: finite=%d monotone=%d  K0=%.3f Eoed=%.0f\n",
                finite, monotone, K0, Eoed);
    check(finite && monotone, "Berlin Sand III oedometer integrates robustly (no divergence)");
    check(K0 > 0.2 && K0 < 0.8 && Eoed > 0.0, "Berlin Sand oedometer gives physical K0, Eoed");
}

} // namespace

int main() {
    test_triaxial_hyperbola();
    test_isotropic_cap();
    test_berlin_sand_oedometer_robust();
    if (g_failures == 0) {
        std::printf("OK: robust multi-surface HS integrator (substepping + Koiter) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
