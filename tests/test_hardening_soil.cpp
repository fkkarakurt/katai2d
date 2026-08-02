// Hardening Soil (HS) -- P2.3a: stress-dependent stiffness laws + drained triaxial
// hyperbolic primary loading (Duncan-Chang core). Verified against the closed-form
// identities of Schanz, Vermeer & Bonnier (1999):
//   - stiffness power law E ~ ((c cos + sigma sin)/(c cos + p_ref sin))^m
//   - failure deviator qf = MC (= sigma3(Kp-1) + 2c sqrt(Kp))
//   - secant at q = qf/2 equals E50 (the definition of E50)
//   - the hyperbola -eps1 = (1/Ei) q/(1-q/qa) round-trips, asymptote qa = qf/Rf
//   - q is capped at qf (perfectly plastic failure plateau), Ei = 2 E50/(2-Rf)
// (See docs/references/hardening-soil-formulation.md.)
#include <katai/materials/hardening_soil.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;

namespace {

constexpr double kPi = 3.14159265358979323846;

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

void test_stiffness_power_law() {
    // Cohesionless sand: clean power law E50(sigma3) = E50_ref (sigma3/p_ref)^m.
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.p_ref = 100.0; p.friction = 35.0 * kPi / 180.0;

    check(close(p.E50(p.p_ref), p.E50_ref), "E50(p_ref) = E50_ref");
    check(close(p.E50(200.0) / p.E50(100.0), std::pow(2.0, p.m)),
          "E50 power law: doubling sigma3 -> 2^m");
    check(close(p.Eur(150.0) / p.E50(150.0), p.Eur_ref / p.E50_ref),
          "Eur/E50 = Eur_ref/E50_ref (3x stiffer unload-reload)");
    std::printf("  power law: E50(100)=%.0f  E50(200)=%.0f  ratio=%.4f (2^%.2f=%.4f)\n",
                p.E50(100.0), p.E50(200.0), p.E50(200.0) / p.E50(100.0), p.m,
                std::pow(2.0, p.m));
}

void test_failure_deviator_is_mc() {
    // qf must equal the Mohr-Coulomb deviator at failure, sigma1f = sigma3 Kp + 2c sqrt(Kp).
    HardeningSoilParams p;
    p.E50_ref = 2.0e4; p.m = 1.0; p.p_ref = 100.0;
    p.cohesion = 10.0; p.friction = 25.0 * kPi / 180.0;
    const double s = std::sin(p.friction);
    const double Kp = (1.0 + s) / (1.0 - s);
    for (double sigma3 : {50.0, 100.0, 300.0}) {
        const double qf_mc = sigma3 * (Kp - 1.0) + 2.0 * p.cohesion * std::sqrt(Kp);
        check(close(p.q_failure(sigma3), qf_mc),
              "qf = MC failure deviator (sigma3(Kp-1)+2c sqrt(Kp))");
    }
    std::printf("  failure: qf(100)=%.2f  MC=%.2f\n", p.q_failure(100.0),
                100.0 * (Kp - 1.0) + 2.0 * p.cohesion * std::sqrt(Kp));
}

void test_triaxial_hyperbola() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.Rf = 0.9;
    const double sigma3 = 100.0;
    const double qf = p.q_failure(sigma3), qa = p.q_asymptote(sigma3);
    const double Ei = p.Ei(sigma3);

    // (1) Initial tangent = Ei = 2 E50/(2-Rf).
    check(close(p.triaxial_tangent(sigma3, 1e-12), Ei, 1e-6),
          "initial triaxial tangent = Ei = 2 E50/(2-Rf)");

    // (2) Secant at q = qf/2 equals E50 (definition of E50). Find -eps1 at 50%.
    const double eps_50 = (qf / 2.0) / (Ei * (1.0 - p.Rf / 2.0));
    check(close(p.triaxial_q(sigma3, eps_50), qf / 2.0, 1e-9),
          "triaxial_q reaches qf/2 at the 50% strain");
    check(close((qf / 2.0) / eps_50, p.E50(sigma3), 1e-9),
          "secant at 50% strength = E50");

    // (3) Hyperbola round-trip for several pre-failure strains.
    double max_err = 0.0;
    for (double e : {1e-4, 5e-4, 1e-3, 3e-3}) {
        const double q = p.triaxial_q(sigma3, e);
        if (q >= qf) continue;  // skip the capped (failure) plateau
        const double eps_back = (1.0 / Ei) * q / (1.0 - q / qa);
        max_err = std::fmax(max_err, std::fabs(eps_back - e));
    }
    check(max_err < 1e-12, "hyperbola -eps1 = (1/Ei) q/(1-q/qa) round-trips");

    // (4) Failure plateau: very large strain caps q at qf (perfectly plastic).
    check(close(p.triaxial_q(sigma3, 1.0), qf, 1e-9), "q capped at qf (failure plateau)");
    check(p.triaxial_tangent(sigma3, 1.0) == 0.0, "post-failure tangent = 0 (perfect plasticity)");

    std::printf("  triaxial (sigma3=%.0f): qf=%.2f  qa=%.2f  Ei=%.0f  E50=%.0f\n",
                sigma3, qf, qa, Ei, p.E50(sigma3));
    std::printf("  hyperbola round-trip max err=%.3e\n", max_err);
}

void test_oedometer() {
    // (a) General c, m: the closed-form 1D compression strain must equal the numerical
    //     integral of the tangent law deps1 = dsigma1/Eoed(sigma1) -- a genuine check
    //     (closed form vs fine-step integration of the modulus law).
    HardeningSoilParams p;
    p.Eoed_ref = 3.0e4; p.m = 0.7; p.p_ref = 100.0;
    p.cohesion = 12.0; p.friction = 28.0 * kPi / 180.0;
    const double s0 = 50.0, s1 = 400.0;
    const int N = 200000;
    double eps_num = 0.0;
    const double ds = (s1 - s0) / N;
    for (int i = 0; i < N; ++i) {
        const double sm = s0 + (i + 0.5) * ds;  // midpoint rule
        eps_num += ds / p.Eoed(sm);
    }
    const double eps_cf = p.oedometer_strain(s0, s1);
    std::printf("  oedometer (c,m=0.7): closed=%.6e  numerical=%.6e\n", eps_cf, eps_num);
    check(close(eps_cf, eps_num, 1e-6), "oedometer closed form = integral of dsigma/Eoed");

    // (b) Tangent: d(eps1)/d(sigma1) = 1/Eoed(sigma1) (central difference).
    const double s = 250.0, h = 1e-3;
    const double dtan = (p.oedometer_strain(s0, s + h) - p.oedometer_strain(s0, s - h)) / (2 * h);
    check(close(dtan, 1.0 / p.Eoed(s), 1e-6), "oedometer tangent modulus = Eoed(sigma1)");

    // (c) m = 1, c = 0: the classic logarithmic (normally consolidated) compression
    //     eps1 = (p_ref/Eoed_ref) ln(sigma1/sigma0).
    HardeningSoilParams clay;
    clay.Eoed_ref = 2.0e3; clay.m = 1.0; clay.p_ref = 100.0;
    clay.friction = 22.0 * kPi / 180.0;  // c = 0
    const double eps_log = clay.oedometer_strain(100.0, 800.0);
    const double eps_log_exact = (clay.p_ref / clay.Eoed_ref) * std::log(800.0 / 100.0);
    check(close(eps_log, eps_log_exact, 1e-9), "m=1,c=0 oedometer is logarithmic");
    std::printf("  oedometer (clay m=1): eps1=%.4e  log-law=%.4e\n", eps_log, eps_log_exact);

    // (d) Eoed depends on the MAJOR stress sigma1 (vs E50/Eur on the minor sigma3).
    check(p.Eoed(400.0) > p.Eoed(50.0), "Eoed increases with vertical stress sigma1");
}

void test_unload_reload() {
    // The defining HS feature vs Duncan-Chang: primary loading is plastic (hyperbola),
    // unload/reload is elastic with the stiff Eur and leaves permanent plastic strain.
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.Rf = 0.9;
    const double sigma3 = 100.0;
    const double qf = p.q_failure(sigma3), Eur = p.Eur(sigma3);

    katai::core::HsTriaxialPath path{p, sigma3, 0.0};

    // Primary loading to q1: axial strain follows the hyperbola.
    const double q1 = 0.6 * qf;
    const double e1 = path.axial_strain(q1);
    check(close(e1, path.hyperbola(q1), 1e-9), "primary loading follows the hyperbola");

    // Unload to q2 < q1: elastic, slope 1/Eur.
    const double q2 = 0.2 * qf;
    const double e2 = path.axial_strain(q2);
    check(close((e1 - e2) / (q1 - q2), 1.0 / Eur, 1e-9), "unload slope = 1/Eur (stiff)");

    // Permanent plastic strain at q = 0 (the model is elastoplastic, not nonlinear elastic).
    const double e0 = path.axial_strain(0.0);
    const double e0_exact = path.hyperbola(q1) - q1 / Eur;
    check(close(e0, e0_exact, 1e-9) && e0 > 0.0, "permanent plastic strain remains at q=0");

    // Reload to q1: elastic back onto the primary curve (rejoin), then beyond -> primary.
    const double e1r = path.axial_strain(q1);
    check(close(e1r, path.hyperbola(q1), 1e-9), "reload rejoins the primary curve at q_max");
    const double q3 = 0.85 * qf;
    const double e3 = path.axial_strain(q3);
    check(close(e3, path.hyperbola(q3), 1e-9), "loading beyond q_max continues on primary");

    // Unload-reload stiffness greatly exceeds the primary secant at q1.
    const double secant1 = q1 / path.hyperbola(q1);
    check(Eur > 2.0 * secant1, "Eur >> primary secant (unload much stiffer)");
    std::printf("  unload-reload: perm plastic eps1=%.4e  Eur=%.0f  secant(q1)=%.0f\n",
                e0, Eur, secant1);
}

} // namespace

int main() {
    test_stiffness_power_law();
    test_failure_deviator_is_mc();
    test_triaxial_hyperbola();
    test_oedometer();
    test_unload_reload();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil stiffness laws + triaxial hyperbola verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
