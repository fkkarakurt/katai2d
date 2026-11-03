// Soft Soil Creep where its creep direction is singular: at the top of the ellipse, q~ = M (p' + c cot phi).
//
// The rate is normalised on the volumetric part of the normal, so the deviatoric creep grows without bound as a stress
// nears the top, and it drives the stress away from it from either side. A step that kept the trial's direction moved by
// that enormous m: from a clay near the top of its ellipse one step gave q = 1.15e7 kPa and went on, reported as solved,
// and in a run of 200 steps q fell to 121 and 23 kPa where the answer is 297. The direction is now read where the step
// ends (soft_soil_creep.cpp), which is a question of one scalar.
//
// verify: KV-CST-028
//   oracle:   independent_path
//   source:   Vermeer, P.A. & Neher, H.P. (1999). A soft soil model that accounts for creep. Beyond 2000 in Computational Geotechnics, Balkema, 249-261, the rate equations integrated independently of the program by an implicit Runge-Kutta scheme (Radau IIA, relative tolerance 1e-10)
//   locator:  eps_v^c rate = (mu*/tau) (p_eq / p_p)^beta, p_eq = p' + q~^2 / (M^2 p'), the flow along d p_eq / d sigma normalised on its trace, p_p = p_p0 exp(eps_v^c / (lam* - kap*)), Mohr-Coulomb after creep; undrained triaxial compression, strain rate 2.16 per day. (a) lam* 0.105, kap* 0.016, mu* 0.004, nu 0.25, phi 32, K0nc 0.61638 (M 1.2916), p' 525, p_p 373; (b) K0nc 0.75 (M 0.934, below the friction line 1.287), nu 0.15, p' 100, p_p = OCR p'
//   quantity: deviator and mean stress at the end of the test [kPa]; the coefficient K0 on a steady one-dimensional path
//   expected: (a) q 298.64, p' 232.75 at 12% strain, peak 326.04; (b) at 6% strain, OCR 1: q 57.25, p' 63.12; OCR 2: q 92.73, p' 99.44; OCR 8: on the friction line, q 128.72, p' 100; K0 = K0nc
//   band:     0.1% on the end states at 12 to 1600 steps (measured +0.04% at 12, +0.01% at 1600); 0.3% on the peak; 0.1% on K0 (measured -0.03%; +2% when the direction was the trial's)
#include <katai/materials/soft_soil_creep.hpp>

#include <cmath>
#include <cstdio>
#include <string>

using namespace katai::core::softsoilcreep;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

struct End { double q = 0.0, p = 0.0, peak = 0.0; };

// Undrained triaxial compression at a rate, in n equal steps.
End undrained(const Params& P, double p0, double pp0, double rate, double strain, int n) {
    Eigen::Vector3d s = Eigen::Vector3d::Constant(p0);
    double pp = pp0;
    const double de = strain / n;
    End out;
    for (int i = 0; i < n; ++i) {
        const StepResult r = ssc_step(P, s, pp, Eigen::Vector3d(de, -de / 2, -de / 2), de / rate);
        s = r.sig;
        pp = r.pp;
        out.peak = std::max(out.peak, s(0) - s(2));
    }
    out.q = s(0) - s(2);
    out.p = s.mean();
    return out;
}

bool within(double got, double want, double band) { return std::fabs(got / want - 1.0) <= band; }
}  // namespace

int main() {
    std::printf("Soft Soil Creep at the top of its ellipse\n\n");
    Params P;
    P.lam_star = 0.105; P.kap_star = 0.016; P.mu_star = 0.004; P.nu_ur = 0.25;
    P.c = 0.0; P.phi = 32.0 * kPi / 180.0; P.psi = 0.0; P.K0nc = 0.6163803116868497; P.tau_day = 1.0;
    std::printf("(a) a clay whose top of the ellipse is just above its friction line\n");
    for (const int n : {12, 48, 100, 200, 400, 1600}) {
        const End e = undrained(P, 525.0, 373.0, 2.16, 0.12, n);
        std::printf("   %5d steps: q %.3f (%+.3f%%), p %.3f (%+.3f%%), peak %.3f (%+.3f%%)\n", n, e.q,
                    100.0 * (e.q / 298.637 - 1.0), e.p, 100.0 * (e.p / 232.745 - 1.0), e.peak, 100.0 * (e.peak / 326.039 - 1.0));
        check(within(e.q, 298.637, 1e-3) && within(e.p, 232.745, 1e-3), std::to_string(n) + " steps: the end state");
        check(within(e.peak, 326.039, 3e-3), std::to_string(n) + " steps: the peak");
    }
    // One long step from a state near the top: the same answer, and nothing like the 1.15e7 it was.
    {
        Eigen::Vector3d s = Eigen::Vector3d::Constant(525.0);
        double pp = 373.0;
        const double de = 0.06 / 2400;
        for (int i = 0; i < 2400; ++i) {
            const StepResult r = ssc_step(P, s, pp, Eigen::Vector3d(de, -de / 2, -de / 2), de / 2.16);
            s = r.sig;
            pp = r.pp;
        }
        const double D = 0.02;
        const StepResult r = ssc_step(P, s, pp, Eigen::Vector3d(D, -D / 2, -D / 2), D / 2.16);
        std::printf("   one step of 2%% from 6%% strain: q %.2f (rate equations 305.65)\n", r.sig(0) - r.sig(2));
        check(within(r.sig(0) - r.sig(2), 305.650, 1e-3), "one long step from near the top");
    }

    std::printf("\n(b) a clay whose top lies below its friction line: only creep that matters holds the stress there\n");
    Params Q = P;
    Q.nu_ur = 0.15; Q.K0nc = 0.75;
    struct Case { double ocr, q, p; };
    for (const Case c : {Case{1.0, 57.25, 63.12}, Case{2.0, 92.73, 99.44}, Case{8.0, 128.72, 100.0}}) {
        const End e = undrained(Q, 100.0, c.ocr * 100.0, 2.16, 0.06, 600);
        std::printf("   OCR %.0f: q %.2f (%.2f), p %.2f (%.2f)\n", c.ocr, e.q, c.q, e.p, c.p);
        check(within(e.q, c.q, 1e-3) && within(e.p, c.p, 1e-3), "OCR " + std::to_string((int)c.ocr));
    }

    std::printf("\n(c) K0 of a steady one-dimensional path, at the substeps of the rule\n");
    Params K;
    K.lam_star = 0.1; K.kap_star = 0.02; K.mu_star = 0.005; K.nu_ur = 0.15;
    K.c = 0.1; K.phi = 30.0 * kPi / 180.0; K.psi = 0.0; K.K0nc = 0.6; K.tau_day = 1.0;
    const katai::core::softsoil::Params KS = K.ss();
    for (const double rate : {0.01, 1.0, 100.0}) {
        Eigen::Vector3d s(100.0, 60.0, 60.0);
        double pp = p_eq(KS, s);
        const double de = 0.4 / 800;
        for (int i = 0; i < 800; ++i) {
            const StepResult r = ssc_step(K, s, pp, Eigen::Vector3d(de, 0.0, 0.0), de / rate);
            s = r.sig;
            pp = r.pp;
        }
        std::printf("   rate %.2f per day: K0 %.5f (K0nc 0.6)\n", rate, s(1) / s(0));
        check(within(s(1) / s(0), 0.6, 1e-3), "K0 at a rate of " + std::to_string(rate));
    }
    if (g_failures == 0) {
        std::printf("\nOK: at the top of the ellipse the integration is the rate equations' at every size of step\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
