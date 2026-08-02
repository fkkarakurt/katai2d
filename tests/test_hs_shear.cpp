// Hardening Soil -- P2.3d (piece 1): multiaxial SHEAR hardening return mapping verified
// by a single-Gauss-point drained triaxial test. The strain-driven predictor-corrector
// (the FE constitutive form) GENERATES the deviatoric stress-strain curve; with zero
// dilatancy (psi=0) the elastoplastic machine must reproduce the closed-form hyperbola
//   eps1 = (1/Ei) q/(1-q/qa)      (compression positive)
// to round-off (because sigma3 is held exactly constant -> Ei,Eur constant, and each
// step solves f=0 so eps1 = q/Eur + eps1^p = hyperbola(q) exactly, independent of step
// size). Also checks the consistent (asymmetric) tangent against finite differences and
// the MC failure plateau. (See docs/references/hardening-soil-formulation.md.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::HsShearStep;
using katai::core::hs_shear_step;

namespace {

constexpr double kPi = 3.14159265358979323846;

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

// Advance one drained-triaxial step: prescribe axial strain increment deps1 (compression
// positive), solve for the lateral strain that keeps sigma3 constant (mixed BC, inner
// Newton on the consistent tangent). Updates sigma, gamma_p, eps1 in place.
void triaxial_step(const HardeningSoilParams& p, double sigma3, double deps1,
                   Eigen::Vector3d& sigma, double& gamma_p, double& eps1) {
    double dlat = 0.0;
    HsShearStep r;
    for (int it = 0; it < 20; ++it) {
        const Eigen::Vector3d deps(deps1, dlat, dlat);
        r = hs_shear_step(p, sigma, gamma_p, deps);
        const double g = r.stress(2) - sigma3;
        if (std::fabs(g) <= 1e-10 * (1.0 + sigma3)) break;
        const double dgdl = r.tangent(2, 1) + r.tangent(2, 2);
        dlat -= g / dgdl;
    }
    sigma = r.stress;
    gamma_p = r.gamma_p;
    eps1 += deps1;
}

void test_triaxial_reproduces_hyperbola() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.dilatancy = 0.0;  // psi=0 -> exact hyperbola
    p.Rf = 0.9; p.nu_ur = 0.2;
    const double sigma3 = 100.0;
    const double qf = p.q_failure(sigma3), Ei = p.Ei(sigma3), qa = p.q_asymptote(sigma3);

    Eigen::Vector3d sigma(sigma3, sigma3, sigma3);  // isotropic start (q=0)
    double gamma_p = 0.0, eps1 = 0.0;

    auto hyperbola = [&](double q) { return (1.0 / Ei) * q / (1.0 - q / qa); };

    // Drive axial strain; at several deviator levels check eps1 == hyperbola(q).
    double max_rel = 0.0;
    int checks = 0;
    double next_frac = 0.2;  // check at 20%, 40%, 60%, 80% of qf
    for (int s = 0; s < 4000; ++s) {
        triaxial_step(p, sigma3, 5.0e-6, sigma, gamma_p, eps1);
        const double q = sigma(0) - sigma(2);
        if (q >= next_frac * qf && next_frac <= 0.81) {
            const double rel = std::fabs(eps1 - hyperbola(q)) / hyperbola(q);
            max_rel = std::fmax(max_rel, rel);
            ++checks;
            std::printf("  q=%.2f (%.0f%% qf): eps1=%.6e  hyperbola=%.6e  rel=%.2e\n",
                        q, 100 * q / qf, eps1, hyperbola(q), rel);
            next_frac += 0.2;
        }
    }
    check(checks == 4, "reached 20/40/60/80% of qf");
    check(max_rel < 1e-6, "triaxial return mapping reproduces the hyperbola (psi=0)");
}

void test_consistent_tangent() {
    // The algorithmic tangent D^ep must match a finite-difference of the stress update
    // at a clearly-plastic state (asymmetric, since flow is non-associated).
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.dilatancy = 8.0 * kPi / 180.0;  // psi>0
    p.Rf = 0.9; p.nu_ur = 0.2;
    const double sigma3 = 100.0;

    // Load to a plastic state first.
    Eigen::Vector3d sigma(sigma3, sigma3, sigma3);
    double gamma_p = 0.0, eps1 = 0.0;
    for (int s = 0; s < 300; ++s) triaxial_step(p, sigma3, 5.0e-6, sigma, gamma_p, eps1);

    // Reference plastic step.
    const Eigen::Vector3d deps0(1.0e-5, -3.0e-6, -3.0e-6);
    const HsShearStep r = hs_shear_step(p, sigma, gamma_p, deps0);
    check(r.plastic, "tangent check is at a plastic state");

    Eigen::Matrix3d fd;
    const double h = 1.0e-9;
    for (int j = 0; j < 3; ++j) {
        Eigen::Vector3d dp = deps0;
        dp(j) += h;
        const HsShearStep rp = hs_shear_step(p, sigma, gamma_p, dp);
        fd.col(j) = (rp.stress - r.stress) / h;
    }
    const double err = (fd - r.tangent).cwiseAbs().maxCoeff() /
                       r.tangent.cwiseAbs().maxCoeff();
    std::printf("  consistent tangent: max rel err vs FD = %.3e\n", err);
    check(err < 1e-5, "consistent tangent matches finite difference");
}

void test_failure_plateau() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.dilatancy = 0.0;
    p.Rf = 0.9; p.nu_ur = 0.2;
    const double sigma3 = 100.0, qf = p.q_failure(sigma3);

    Eigen::Vector3d sigma(sigma3, sigma3, sigma3);
    double gamma_p = 0.0, eps1 = 0.0;
    for (int s = 0; s < 20000; ++s) triaxial_step(p, sigma3, 5.0e-6, sigma, gamma_p, eps1);
    const double q = sigma(0) - sigma(2);
    std::printf("  failure plateau: q=%.4f  qf=%.4f  (eps1=%.3f)\n", q, qf, eps1);
    check(close(q, qf, 1e-4), "large strain -> q reaches qf (MC failure plateau)");
}

} // namespace

int main() {
    test_triaxial_reproduces_hyperbola();
    test_consistent_tangent();
    test_failure_plateau();
    if (g_failures == 0) {
        std::printf("OK: HS shear return mapping (triaxial hyperbola + tangent) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
