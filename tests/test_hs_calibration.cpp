// Hardening Soil cap calibration (P2.3d-2b) -- the calibration step for the
// compression side. A user enters the standard HS inputs (E50_ref, Eoed_ref, Eur_ref, m,
// c, phi, psi, p_ref) and K0^NC; (alpha, beta) are derived so an oedometer reproduces
// BOTH K0^NC and Eoed_ref. The cap hardening modulus is anchored by the elastic-plastic
// closed form K_p = K1 K2/(K1-K2) (K1=Eur_ref/(3(1-2nu)), K2=Eoed_ref(1+2K0nc)/3), and the
// cap deviatoric measure is the SYMMETRIC von Mises q (q^2=3J2) -- the asymmetric
// delta-qtilde broke the oedometer K0. For a general soil (Eoed_ref < E50_ref) both
// targets are matched essentially exactly. (See hardening-soil-formulation.md sec 4.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

void test_calibrate_general_soil() {
    // A typical clay (Eoed_ref < E50_ref < Eur_ref), m=1.
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 1.2e4; p.Eur_ref = 6.0e4; p.Eoed_ref = 8.0e3;
    p.m = 1.0; p.nu_ur = 0.2; p.friction = 24 * kPi / 180;
    p.dilatancy = 0.0; p.cohesion = 5.0; p.Rf = 0.9;
    const double K0_NC = 1.0 - std::sin(p.friction);  // Jaky 0.5933

    katai::core::hs_calibrate_cap(p, K0_NC);
    double Eoed, K0;
    katai::core::hs_oedometer_probe(p, Eoed, K0);
    std::printf("  clay: alpha=%.4f beta=%.4e -> K0=%.4f (target %.4f) Eoed=%.0f (target %.0f)\n",
                p.cap_alpha, p.cap_beta, K0, K0_NC, Eoed, p.Eoed_ref);
    check(p.cap_alpha > 0 && p.cap_beta > 0, "calibrated alpha, beta physical");
    check(close(K0, K0_NC, 5e-3), "oedometer reproduces K0^NC (von Mises cap)");
    check(close(Eoed, p.Eoed_ref, 1e-2), "oedometer reproduces Eoed_ref (closed-form Kp)");
}

void test_closed_form_Kp() {
    // The cap hardening anchor K_p = K1 K2/(K1-K2): the plastic bulk modulus that, in series
    // with the elastic K1, gives the oedometric bulk modulus K2.
    HardeningSoilParams p;
    p.Eur_ref = 3.15e5; p.Eoed_ref = 1.05e5; p.nu_ur = 0.2;
    const double K0 = 0.38;
    const double K1 = p.Eur_ref / (3 * (1 - 2 * p.nu_ur));
    const double K2 = p.Eoed_ref * (1 + 2 * K0) / 3;
    const double Kp = K1 * K2 / (K1 - K2);
    check(close(katai::core::hs_cap_Kp(p, K0), Kp, 1e-9), "hs_cap_Kp = K1 K2/(K1-K2)");
}

// Berlin Sand III (Eoed_ref = E50_ref, very stiff cap, phi = 38 deg), K0^NC = 0.38.
//
// CORRECTED 2026-09-27. This used to assert that 0.38 is outside what the model can reproduce --
// "a property of the model, not a solver defect" -- and that the nearest the calibration could
// get was ~0.42. It was a property of the CALIBRATION: the oedometer probe it bisects on started
// from an isotropic state at 2% of p_ref and read K0 at p_ref as though the path had forgotten
// that start. It had not (at alpha = 2.4 the reading was 0.405 from a 2 kPa start, 0.453 from a
// 20 kPa one), and no alpha brought the reading down to 0.38. Started on the K0^NC line itself,
// the probe reads 0.38 exactly for alpha = 2.16, and a normally consolidated walk holds that K0
// at every stress level.
void test_berlin_k0nc_reached() {
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 105e3; p.Eur_ref = 315e3; p.Eoed_ref = 105e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38 * kPi / 180;
    p.dilatancy = 6 * kPi / 180; p.cohesion = 1.0; p.Rf = 0.9;
    katai::core::hs_calibrate_cap(p, 0.38);
    double Eoed, K0;
    katai::core::hs_oedometer_probe(p, Eoed, K0);
    std::printf("  Berlin Sand III: alpha=%.3f beta=%.3e -> K0=%.4f (input 0.38)"
                " Eoed=%.0f (target 105000)\n", p.cap_alpha, p.cap_beta, K0, Eoed);
    check(close(Eoed, p.Eoed_ref, 1e-2), "Berlin: Eoed_ref reproduced exactly");
    check(std::fabs(K0 - 0.38) < 1e-3, "Berlin: the input K0^NC = 0.38 is reproduced");
    // Away from the start and from p_ref: a normally consolidated oedometer walk keeps K0^NC,
    // which is what the calibration is for (c = 1 kPa, so only nearly scale-free).
    Eigen::Vector3d s(60.0, 60.0 * 0.38, 60.0 * 0.38);
    double gp = katai::core::hs_initial_gamma_p(p, s), pp = katai::core::hs_initial_pp(p, s);
    while (s(0) < 800.0) {
        const auto r = katai::core::hs_integrate(p, s, gp, pp, Eigen::Vector3d(2e-4, 0, 0));
        s = r.stress; gp = r.gamma_p; pp = r.pp;
    }
    std::printf("  Berlin Sand III: K0 at sigma1 = %.0f kPa on the NC line %.4f\n", s(0), s(2) / s(0));
    check(std::fabs(s(2) / s(0) - 0.38) < 5e-3, "Berlin: K0^NC held along the normally consolidated line");
}

} // namespace

int main() {
    test_closed_form_Kp();
    test_calibrate_general_soil();
    test_berlin_k0nc_reached();
    if (g_failures == 0) {
        std::printf("OK: HS cap calibration (K0^NC + Eoed_ref) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
