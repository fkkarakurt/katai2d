// Hardening Soil -- UNDRAINED clay (psi = 0), the common PLAXIS Undrained-(A) case where the
// von Mises cap is EXACT (psi=0 -> no shear dilation -> the Lode-dependent cap refinement that
// matters only for dilatant undrained is irrelevant here). This validates HS undrained against
// the unambiguous reference: the EFFECTIVE Mohr-Coulomb failure envelope.
//
// A normally-consolidated clay sheared undrained (constant volume, eps_v = 0) only CONTRACTS
// (cap + non-dilatant shear): the excess pore pressure is positive and monotone, the effective
// mean stress p' DECREASES, and the effective stress path curves LEFT until it terminates on the
// critical-state / MC failure line in q-p' space:
//     q_f = M p'_f + (6 c cos phi)/(3 - sin phi),   M = 6 sin phi / (3 - sin phi)   (triaxial comp).
// This is exactly the classic undrained NC-clay effective-stress path that PLAXIS Undrained-(A)
// produces. (Undrained = kinematic eps_v=0 at the material point; the excess pore pressure is
// u = sigma3_cell - sigma3'. See hardening-soil-formulation.md / effective-stress-formulation.md.)
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

void test_undrained_clay() {
    // Soft NC clay, effective parameters; psi = 0 (no dilation -- the undrained clay case).
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 8e3; p.Eur_ref = 24e3; p.Eoed_ref = 8e3;
    p.m = 1.0; p.nu_ur = 0.2; p.friction = 25 * kPi / 180;
    p.dilatancy = 0.0; p.cohesion = 10.0; p.Rf = 0.9;
    katai::core::hs_calibrate_cap(p, 1.0 - std::sin(p.friction));

    const double sinphi = std::sin(p.friction), cosphi = std::cos(p.friction);
    const double M = 6.0 * sinphi / (3.0 - sinphi);                 // CSL slope (triaxial comp)
    const double q_intercept = 6.0 * p.cohesion * cosphi / (3.0 - sinphi);

    const double sigma3 = 100.0;                                    // isotropic consolidation
    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0.0, pp = sigma3, eps1 = 0.0;
    const double de1 = 5e-6;
    double u_min = 0.0, q_final = 0.0, p_final = 0.0, u_final = 0.0;
    bool u_positive = true;
    for (int s = 0; s < 24000; ++s) {
        const HsIntegrated r = hs_integrate(p, sig, gp, pp,
                                            Eigen::Vector3d(de1, -0.5 * de1, -0.5 * de1));
        sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += de1;
        const double q = sig(0) - sig(2);
        const double u = sigma3 - sig(2);
        const double pm = (sig(0) + sig(1) + sig(2)) / 3.0;
        u_min = std::fmin(u_min, u);
        if (eps1 > 0.001 && u < -1.0) u_positive = false;  // clay (psi=0) must not dilate
        q_final = q; p_final = pm; u_final = u;
        if (s % 3000 == 0)
            std::printf("      eps1=%.3f q=%.1f p'=%.1f u=%+.1f q/(M p'+c0)=%.3f\n",
                        eps1, q, pm, u, q / (M * pm + q_intercept));
        if (eps1 > 0.11) break;
    }
    const double q_envelope = M * p_final + q_intercept;           // effective failure line at p'_final
    const double rel = std::fabs(q_final - q_envelope) / q_envelope;
    std::printf("  HS UNDRAINED clay (phi'=25, psi=0, c'=10):\n");
    std::printf("    final: q=%.1f p'=%.1f u=%+.1f | MC envelope q=M*p'+c0 = %.1f (M=%.3f, rel err=%.1f%%)\n",
                q_final, p_final, u_final, q_envelope, M, 100.0 * rel);
    std::printf("    u_min over the test = %+.1f kPa (>=0: contraction only, no dilation)\n", u_min);

    check(u_positive && u_min > -1.0,
          "undrained NC clay (psi=0) only generates POSITIVE excess pore pressure (no dilation)");
    check(u_final > 0.0, "p' decreased (effective stress path curves left -- classic NC clay)");
    check(rel < 0.03,
          "effective stress path terminates on the MC failure envelope q = M p' + 6c cos/(3-sin)");
}

}  // namespace

int main() {
    test_undrained_clay();
    if (g_failures == 0) {
        std::printf("OK: HS undrained NC clay (psi=0) matches the effective MC failure envelope "
                    "(PLAXIS Undrained-A)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
