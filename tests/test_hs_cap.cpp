// Hardening Soil -- P2.3d (piece 2): cap (volumetric) yield surface return mapping.
// f_c = qtil^2/alpha^2 + p^2 - pp^2, associated flow, hardening pp_dot = Kc * evp_dot.
// Verified by the closed-form ISOTROPIC compression response: with qtil = 0 only the cap
// is active and the volumetric tangent is the series of the elastic and cap moduli,
//   K_iso = K_e * Kc / (K_e + Kc)
// (a clean, calibration-independent check of the cap return mapping + hardening). Also
// checks return-to-surface (f_c = 0) and the consistent (symmetric, associated) tangent
// against finite differences. (See docs/references/hardening-soil-formulation.md sec 4.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::HsCapStep;
using katai::core::hs_cap_step;

namespace {

constexpr double kPi = 3.14159265358979323846;

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

HardeningSoilParams make_params() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.p_ref = 100.0; p.friction = 30.0 * kPi / 180.0;
    p.nu_ur = 0.2; p.cap_alpha = 0.6; p.cap_beta = 2.0e-3;
    return p;
}

// f_c for assertions.
double fcap(const HardeningSoilParams& p, const Eigen::Vector3d& s, double pp) {
    const double sphi = std::sin(p.friction);
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double pm = (s(0) + s(1) + s(2)) / 3.0;
    const double qt = s(0) + (delta - 1.0) * s(1) - delta * s(2);
    return qt * qt / (p.cap_alpha * p.cap_alpha) + pm * pm - pp * pp;
}

void test_isotropic_compression() {
    const HardeningSoilParams p = make_params();
    const double p0 = 100.0;                 // start on the cap (pp = p0)
    const double Ke = p.bulk(p0);            // = 9e4/(3*0.6) = 5e4

    Eigen::Vector3d sigma(p0, p0, p0);
    const double e = 1.0e-6;                 // small isotropic compression increment
    const Eigen::Vector3d deps(e, e, e);
    const HsCapStep r = hs_cap_step(p, sigma, p0, deps);
    check(r.plastic, "isotropic loading beyond pp activates the cap");
    check(std::fabs(fcap(p, r.stress, r.pp)) < 1e-6 * p0 * p0,
          "stress returns to the cap surface (f_c = 0)");

    // Power-law cap hardening modulus at the (returned) state: H_cap = (pref/beta)(pc/pref)^m.
    const double Hcap = p.cap_hardening_modulus(r.pp);
    const double K_iso = Ke * Hcap / (Ke + Hcap);  // series of elastic and cap moduli

    // (1) Consistent tangent K = (1/9) sum(D^ep) matches K_iso (evaluated at the state).
    const double K_tan = r.tangent.sum() / 9.0;
    check(close(K_tan, K_iso, 1e-6), "consistent tangent = K_e H_cap/(K_e+H_cap)");

    // (2) Return secant ~ K_iso (looser: H_cap varies over the step).
    const double dp = (r.stress(0) + r.stress(1) + r.stress(2)) / 3.0 - p0;
    check(close(dp / (3.0 * e), K_iso, 1e-3), "isotropic return secant ~ K_iso");

    std::printf("  isotropic cap: K_e=%.0f H_cap=%.0f -> K_iso=%.1f (return %.1f, tangent %.1f)\n",
                Ke, Hcap, K_iso, dp / (3.0 * e), K_tan);
}

void test_cap_consistent_tangent() {
    const HardeningSoilParams p = make_params();
    Eigen::Vector3d sigma(100.0, 100.0, 100.0);
    const double pp = 100.0;
    const Eigen::Vector3d deps0(2.0e-4, 5.0e-5, 5.0e-5);  // volumetric + deviatoric
    const HsCapStep r = hs_cap_step(p, sigma, pp, deps0);
    check(r.plastic, "cap tangent check at a plastic state");

    Eigen::Matrix3d fd;
    const double h = 1.0e-9;
    for (int j = 0; j < 3; ++j) {
        Eigen::Vector3d dp = deps0;
        dp(j) += h;
        const HsCapStep rp = hs_cap_step(p, sigma, pp, dp);
        fd.col(j) = (rp.stress - r.stress) / h;
    }
    const double err = (fd - r.tangent).cwiseAbs().maxCoeff() /
                       r.tangent.cwiseAbs().maxCoeff();
    const double asym = (r.tangent - r.tangent.transpose()).cwiseAbs().maxCoeff() /
                        r.tangent.cwiseAbs().maxCoeff();
    std::printf("  cap tangent: max rel err vs FD=%.3e  asymmetry=%.3e\n", err, asym);
    check(err < 1e-5, "cap algorithmic tangent matches finite difference");
    // Associated flow -> near-symmetric; a small asymmetry is the algorithmic
    // hardening-stress coupling (pp depends on p within the step), vanishing as lambda->0.
    check(asym < 1e-2, "cap tangent is (near-)symmetric (associated flow)");
}

} // namespace

int main() {
    test_isotropic_compression();
    test_cap_consistent_tangent();
    if (g_failures == 0) {
        std::printf("OK: HS cap yield surface (isotropic compression + tangent) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
