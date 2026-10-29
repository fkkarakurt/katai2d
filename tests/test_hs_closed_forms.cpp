// Hardening Soil against the four results of the model that can be written down.
//
// The model has closed forms only where one of its two mechanisms works alone, or neither:
//
//   (1) DRAINED TRIAXIAL COMPRESSION with the cap far away and no dilatancy is the hyperbola the
//       model is built on: q = E_i eps_1 / (1 + E_i eps_1 / q_a), E_i = 2 E50 / (2 - R_f),
//       q_a = q_f / R_f. The confining stress is constant, so the stiffness is, the plastic
//       volume change is zero, and gamma_p = 2 eps_1^p exactly.
//       E_i is NOT a free choice of this file, and the sources do not write it one way. E50 is
//       defined in words as the secant at 50% of the failure deviator q_f (Schanz 1998, Sec. 4.2.4,
//       p. 81, Eqn 4.31); putting q = q_f / 2 into the hyperbola and asking for that secant gives
//       E_i = 2 E50 / (2 - R_f). The same text's yield function (Schanz 1998 Eqn 4.34; Benz 2007
//       Eqns 7.10-7.14) is written with q_a / E50, i.e. E_i = 2 E50 -- the secant at HALF THE
//       ASYMPTOTE q_a / 2, not at q_f / 2 -- and so contradicts its own definition by the factor
//       2 - R_f (1.1 at R_f = 0.9). This build keeps the definition in words, which is the one a
//       triaxial test is read by; (1) pins it, and fails by 10% if the other is ever taken.
//   (2) THE CAP'S STEEPNESS for m = 1 and c = 0. The shear function then depends on the stress
//       ratio alone, a normally consolidated oedometer path keeps that ratio, the shear mechanism
//       does not flow on it, and the two conditions of the oedometer (d eps_q / d eps_v = 2/3
//       and the lateral ratio K0nc) give alpha^2 = (1 - r) eta / (2/3 - r eta K / (3 G)) with
//       eta = 3 (1 - K0nc) / (1 + 2 K0nc), r = Kc / Ks = (E_oed / E_ur) (1 + 2 K0nc) (1 - 2 nu) / K0nc
//       and K / (3 G) = 2 (1 + nu) / (9 (1 - 2 nu)).
//   (3) UNDRAINED ELASTIC LOADING with m = 1. The volume is constant, so p' is; the confining
//       stress is p' - q / 3 and the shear modulus is proportional to it, which integrates to
//       sigma_3 = sigma_3,0 exp(-G_ref eps_1 / p_ref), q = 3 (p' - sigma_3).
//   (4) WHERE THE UNDRAINED PATH MEETS THE CAP. With no dilatancy the path is vertical in p'-q
//       until it does, at q^2 = q_c^2 + alpha^2 (p_c^2 - p_0^2), (p_c, q_c) the pre-consolidation
//       state on the normally consolidated line.
//
// Where both mechanisms work together -- undrained shearing past (4) -- there is no closed form,
// and that is where the published curves of other programs are the only comparison (KV-CST-019).
//
// verify: KV-CST-020
//   oracle:   closed_form
//   source:   Schanz, T. (1998). Zur Modellierung des mechanischen Verhaltens von Reibungsmaterialien. Mitteilung 45, Institut fuer Geotechnik, Universitaet Stuttgart, Sec. 4.2; Benz, T. (2007). Small-strain stiffness of soils and its numerical consequences. Mitteilung 55, Institut fuer Geotechnik, Universitaet Stuttgart, Sec. 7.2
//   locator:  Schanz (1998) Eqn 4.31 with the definition of E50 as the secant at 50% of q_f on p. 81 (from which E_i = 2 E50 / (2 - R_f), as derived at the head of this file; the yield function of Schanz Eqn 4.34 and Benz Eqns 7.10-7.14 reads E_i = 2 E50, the secant at q_a / 2, and is the one the test would fail by 10%), Benz (2007) Eqns 7.10-7.11 (the hyperbola, q_f, q_a) and Benz (2007) Eqns 7.13-7.26 with Box 7.1 (shear function, cap q~^2/alpha^2 + p^2 = p_p^2, associated cap flow, isotropic elasticity with E_ur proportional to sigma_3^m); the four closed forms follow from them as derived at the head of this file
//   quantity: (1) deviator q at five axial strains of a drained triaxial compression [kPa]; (2) cap steepness alpha [-]; (3) q and p' at four axial strains of an undrained elastic compression [kPa]; (4) the deviator at which the pre-consolidation pressure first grows in undrained compression [kPa]
//   expected: (1) 54.79 / 100.00 / 137.93 / 170.21 / 192.77 at 0.2 / 0.5 / 1 / 2 / 4 % (E50 = 20 MPa, phi = 30, sigma_3 = 100); (2) 0.6599 (K0nc = 0.53, nu = 0.15, E_ur / E_oed = 21000 / 5285); (3) 185.75 / 191.92 / 210.11 / 239.33 at 0.01 / 0.02 / 0.05 / 0.1 % with p' = 288.32 throughout; (4) 236.6
//   band:     (1) 0.1% of q; (2) 0.5% of alpha; (3) 0.01% of q and of p'; (4) the closed form lies within 0.1 kPa of the interval between the last deviator at which p_p has not moved and the first at which it has (the two are 0.05 kPa apart; measured 236.58 to 236.63 against 236.66)
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>

#include <cmath>
#include <cstdio>

using namespace katai::core;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what);
    if (!ok) ++g_failures;
}

MaterialModel soil(double E50, double Eoed, double Eur, double m, double nu, double phi_deg,
                   double k0nc) {
    MaterialParams p;
    p.E50_ref = E50; p.Eoed_ref = Eoed; p.Eur_ref = Eur; p.m = m; p.p_ref = 100;
    p.Rf = 0.9; p.nu_ur = nu; p.c = 0; p.phi_rad = phi_deg * kPi / 180; p.psi_rad = 0;
    p.k0nc_auto = false; p.k0nc = k0nc;
    p.tension_cutoff = true; p.tensile_strength = 0;
    p.E = Eur; p.nu = nu;
    return find_model("HardeningSoil")->build(p);
}

// Axisymmetric stress point, tension positive: stress = (radial, axial, shear), stress_zz = hoop.
GaussState state(double s_axial, double s_radial, double gamma_p, double pp) {
    GaussState gs;
    gs.stress = Eigen::Vector3d(-s_radial, -s_axial, 0.0);
    gs.stress_zz = -s_radial;
    gs.gamma_p = gamma_p;
    gs.pp = pp;
    return gs;
}
double q_of(const GaussState& gs) { return -gs.stress(1) + gs.stress(0); }
double p_of(const GaussState& gs) { return -(2.0 * gs.stress(0) + gs.stress(1)) / 3.0; }

GaussState step(const MaterialModel& m, const GaussState& gs, double de_axial, double de_radial) {
    GaussState tr = gs;
    Eigen::Matrix4d T;
    const Eigen::Vector4d de(-de_radial, -de_axial, 0.0, -de_radial);   // compression = negative strain
    integrate_point_axisym(m, gs, de, tr, T, TangentMode::kNone);
    return tr;
}

// (1) drained triaxial compression: the radial strain of each step is the one that brings the
// radial stress back to the cell pressure (a secant on the two-point response).
void test_drained_hyperbola() {
    std::printf("(1) drained triaxial compression, cap far away\n");
    const MaterialModel m = soil(20000, 20000, 60000, 0.5, 0.2, 30, 0.5);
    GaussState gs = state(100, 100, 0.0, 1.0e6);
    const double Ei = 2.0 * 20000 / (2.0 - 0.9), qa = 200.0 / 0.9;
    const double targets[] = {0.002, 0.005, 0.01, 0.02, 0.04};
    const double de = 2.0e-6;
    double eps = 0.0;
    for (double t : targets) {
        while (eps < t - 1e-12) {
            // ... measured against the cell pressure itself, so that nothing accumulates
            const double a = step(m, gs, de, 0.0).stress(0) + 100.0;
            const double b = step(m, gs, de, -1.0e-7).stress(0) + 100.0;
            const double d3 = -a * (-1.0e-7) / (b - a);
            gs = step(m, gs, de, d3);
            eps += de;
        }
        const double exact = Ei * t / (1.0 + Ei * t / qa);
        char msg[200];
        std::snprintf(msg, sizeof msg, "eps_1 = %.1f%%: q = %.2f, the hyperbola gives %.2f (%+.3f%%), sigma_3 = %.3f",
                      100 * t, q_of(gs), exact, 100 * (q_of(gs) / exact - 1), -gs.stress(0));
        check(std::fabs(q_of(gs) / exact - 1) < 1e-3 && std::fabs(-gs.stress(0) - 100.0) < 1e-3, msg);
    }
}

// (2) the cap's steepness for m = 1, c = 0
double alpha_closed(double k0nc, double nu, double Eur_over_Eoed) {
    const double eta = 3.0 * (1.0 - k0nc) / (1.0 + 2.0 * k0nc);
    const double r = (1.0 + 2.0 * k0nc) * (1.0 - 2.0 * nu) / (Eur_over_Eoed * k0nc);   // Kc / Ks
    const double k3g = 2.0 * (1.0 + nu) / (9.0 * (1.0 - 2.0 * nu));
    return std::sqrt((1.0 - r) * eta / (2.0 / 3.0 - r * eta * k3g));
}
void test_cap_steepness() {
    std::printf("(2) the cap's steepness, m = 1 and c = 0\n");
    const MaterialModel m = soil(10500, 5285, 21000, 1.0, 0.15, 30.5, 0.53);
    const double exact = alpha_closed(0.53, 0.15, 21000.0 / 5285.0);
    char msg[160];
    std::snprintf(msg, sizeof msg, "alpha calibrated on the oedometer %.4f, closed form %.4f (%+.2f%%)",
                  m.hs.cap_alpha, exact, 100 * (m.hs.cap_alpha / exact - 1));
    check(std::fabs(m.hs.cap_alpha / exact - 1) < 5e-3, msg);
}

// (3) undrained elastic compression, m = 1: both surfaces far away
void test_undrained_elastic() {
    std::printf("(3) undrained elastic compression, m = 1\n");
    const MaterialModel m = soil(10500, 5285, 21000, 1.0, 0.15, 30.5, 0.53);
    const double s1 = 408.0, s3 = 228.48, p0 = (s1 + 2 * s3) / 3.0;
    GaussState gs = state(s1, s3, 1.0e3, 1.0e6);
    const double Gref = 21000.0 / (2.0 * 1.15);
    const double targets[] = {1e-4, 2e-4, 5e-4, 1e-3};
    const double de = 1.0e-7;
    double eps = 0.0;
    for (double t : targets) {
        while (eps < t - 1e-13) { gs = step(m, gs, de, -de / 2.0); eps += de; }
        const double s3x = s3 * std::exp(-Gref * t / 100.0), qx = 3.0 * (p0 - s3x);
        char msg[200];
        std::snprintf(msg, sizeof msg, "eps_1 = %.2f%%: q = %.3f against %.3f, p' = %.3f against %.3f",
                      100 * t, q_of(gs), qx, p_of(gs), p0);
        check(std::fabs(q_of(gs) / qx - 1) < 1e-4 && std::fabs(p_of(gs) / p0 - 1) < 1e-4, msg);
    }
}

// (4) the deviator at which the undrained path meets the cap
void test_cap_knee() {
    std::printf("(4) where the undrained path meets the cap\n");
    const MaterialModel m = soil(10500, 5285, 21000, 1.0, 0.15, 30.5, 0.53);
    const double s1 = 408.0, s3 = 228.48, p0 = (s1 + 2 * s3) / 3.0;
    const double s1c = 465.1, s3c = 0.53 * s1c, pc = (s1c + 2 * s3c) / 3.0, qc = s1c - s3c;
    const double a = m.hs.cap_alpha;
    const double pp0 = std::sqrt(pc * pc + qc * qc / (a * a));
    const double exact = std::sqrt(qc * qc + a * a * (pc * pc - p0 * p0));
    // the shear surface far away: what is asked is the cap alone
    GaussState gs = state(s1, s3, 1.0e3, pp0);
    const double de = 1.0e-6;
    double q_before = q_of(gs), q_after = -1.0, p_at = 0.0;
    for (int i = 0; i < 4000 && q_after < 0.0; ++i) {
        const GaussState n = step(m, gs, de, -de / 2.0);
        if (n.pp > pp0 * (1.0 + 1e-9)) { q_after = q_of(n); p_at = p_of(gs); }
        else q_before = q_of(n);
        gs = n;
    }
    char msg[220];
    std::snprintf(msg, sizeof msg, "p_p first grows between q = %.2f and %.2f; the closed form is %.2f "
                  "(p_p = %.2f, alpha = %.4f); p' there %.3f against %.3f",
                  q_before, q_after, exact, pp0, a, p_at, p0);
    check(q_after > 0.0 && exact >= q_before - 0.1 && exact <= q_after + 0.1 && q_after - q_before < 0.5 &&
              std::fabs(p_at / p0 - 1) < 1e-4, msg);
}
}  // namespace

int main() {
    test_drained_hyperbola();
    test_cap_steepness();
    test_undrained_elastic();
    test_cap_knee();
    if (g_failures) std::printf("%d check(s) failed\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
