// LIQUEFACTION TRIGGERING BY THE SIMPLIFIED METHOD OF TBDY 2018 ANNEX 16B, against a worked example
// published by the authority that issues the code.
//
// The oracle is the AFAD seminar of worked examples for TBDY 2018 Chapter 16 (IMO Ankara,
// 26.01.2023, slides 41-50): a boring with its water table at 1.80 m, unit weights 17 above and 20
// below it (water 9.8), S_DS = 0.75, Mw = 7.5, the rod-length factor C_R given per depth, and the
// whole chain printed at 6.00 m and summarised for every depth. Its values are printed to one or
// two decimals and its tau_R at 6.00 m was formed from CRR rounded to 0.304, so each quantity is
// compared at the precision it was printed with. One printed value is a slip of the slide's own
// arithmetic -- alpha = 2.739 where Eq. 16B.3b gives 2.767 -- and is checked against the equation,
// with the slip stated.
//
// The screening the slides apply is the code's: 4.50 m is not evaluated because N1,60 >= 30
// (16.6.5); 12.00 m (PI 16 %) is not a liquefiable soil (16.6.2). The slides evaluate 1.50 m too,
// although it lies above the water table, which 16.6.2 excludes: it is reported here as not
// liquefiable, and the departure is stated rather than copied.
//
// The magnitude factor away from 7.5 is checked against the values Kara & Erdem (2023) print for
// Mw 6.0 / 6.5 / 7.0 (1.77 / 1.44 / 1.19, two decimals).
//
// verify: KV-LIQ-001
//   oracle:   published_benchmark
//   source:   AFAD (2023). TBDY-2018 Chapter 16 (foundation soils and foundations under earthquake), worked examples; seminar slides in Turkish, IMO Ankara, 26.01.2023, slides 41-50 (https://ankara.imo.org.tr/Eklenti/8164,16-examples-imo-ankarapdf.pdf); TBDY 2018 (Resmi Gazete 30364 bis, 18.03.2018) Annex 16B Eq. 16B.1-16B.6 and Eq. 16.3; Kara, E. & Erdem, M. (2023). J. Sci. Tech. 2(2), 80-91, p. 89 for C_M
//   locator:  slide 44 (boring log, sigma_v, sigma'_v, C_N, C_R), slides 45-47 (the chain at 6.00 m), slide 48 (results at every depth), slide 49 (screening)
//   quantity: C_N, N1,60, N1,60f, CRR_M7.5, tau_R and tau_deprem [kPa], r_d, the factor tau_R / tau_deprem, and the verdict, at 3.00, 4.50, 6.00, 9.00, 10.50 and 12.00 m; C_M at Mw 6.0, 6.5, 7.0
//   expected: 6.00 m: C_N 1.141, N1,60 21.7, N1,60f 25.6, CRR 0.304, tau_R 22.3, r_d 0.954, tau_deprem 21.3, factor 1.05 (liquefies: below 1.10); 3.00 / 9.00 / 10.50 m: factor 0.90 / 0.83 / 0.79; 4.50 m not required; 12.00 m not a liquefiable soil; C_M 1.77 / 1.44 / 1.19
//   band:     half a unit of the last printed digit, plus 0.1 kPa on tau_R (the slides formed it from CRR rounded to 0.304); alpha and beta against Eq. 16B.3b itself (the slide's 2.739 is a slip: exp(1.0178) = 2.767)
#include <katai/analysis/liquefaction.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace lq = katai::core::liquefaction;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
bool near(double got, double printed, double half_unit) { return std::fabs(got - printed) <= half_unit + 1e-12; }

// The boring of slide 44: water table 1.80 m, gamma 17 above and 20 below, gamma_w 9.8.
lq::SptTest at(double z, double N, double fines, double CR) {
    lq::SptTest t;
    t.depth = z; t.N = N; t.fines = fines;
    const double wt = 1.80;
    t.sigma_v = 17.0 * std::fmin(z, wt) + 20.0 * std::fmax(0.0, z - wt);
    t.sigma_v_eff = t.sigma_v - 9.8 * std::fmax(0.0, z - wt);
    t.c.CR = CR;
    t.below_water = z > wt;
    return t;
}
}  // namespace

int main() {
    std::printf("TBDY 2018 Annex 16B against the AFAD worked example\n\n");
    lq::Design d; d.SDS = 0.75; d.Mw = 7.5;

    {   // the chain at 6.00 m
        const lq::SptTest t = at(6.00, 20, 16, 0.95);
        const lq::Result r = lq::evaluate(t, d);
        std::printf("   6.00 m: sigma_v %.1f, sigma'_v %.1f, C_N %.4f, N1,60 %.3f, alpha %.4f, beta %.4f, "
                    "N1,60f %.3f, CRR %.4f, C_M %.5f, tau_R %.3f, r_d %.4f, tau_eq %.3f, FS %.4f\n",
                    t.sigma_v, t.sigma_v_eff, r.CN, r.N160, r.alpha, r.beta, r.N160f, r.CRR75, r.CM,
                    r.tau_R, r.rd, r.tau_eq, r.FS);
        check(near(t.sigma_v, 114.6, 0.05) && near(t.sigma_v_eff, 73.4, 0.05), "6.00 m: sigma_v 114.6, sigma'_v 73.4");
        check(near(r.CN, 1.141, 0.0005), "6.00 m: C_N = 1.141 (Eq. 16B.2)");
        check(near(r.N160, 21.7, 0.05), "6.00 m: N1,60 = 21.7 (Eq. 16B.1)");
        // The slide prints alpha = 2.739 and beta = 1.053; Eq. 16B.3b gives exp(1.76 - 190/256)
        // = exp(1.0178) = 2.767 and 0.99 + 64/1000 = 1.054 -- a slip in the slide's arithmetic,
        // which N1,60f = 25.6 does not show at its printed precision.
        check(near(r.alpha, 2.7671, 0.00005) && near(r.beta, 1.054, 0.0000001),
              "6.00 m: alpha = exp(1.76 - 190/256) = 2.767, beta = 1.054 (Eq. 16B.3b; the slide's 2.739 / 1.053 are a slip)");
        check(near(r.N160f, 25.6, 0.05), "6.00 m: N1,60f = 25.6 (Eq. 16B.3a)");
        // The slide forms CRR from N1,60f rounded to 25.6; carried unrounded it is 0.3047.
        check(near(lq::crr75(25.6), 0.304, 0.0005) && near(r.CRR75, 0.304, 0.001),
              "6.00 m: CRR_M7.5 = 0.304 from the printed N1,60f (Eq. 16B.4b)");
        check(near(r.tau_R, 22.3, 0.15), "6.00 m: tau_R = 22.3 kPa (Eq. 16B.4a)");
        check(near(r.rd, 0.954, 0.0005), "6.00 m: r_d = 0.954 (Eq. 16B.6)");
        check(near(r.tau_eq, 21.3, 0.05), "6.00 m: tau_deprem = 21.3 kPa (Eq. 16B.5)");
        check(near(r.FS, 1.05, 0.005) && r.verdict == lq::Verdict::Triggered,
              "6.00 m: tau_R / tau_deprem = 1.05 < 1.10: liquefies (Eq. 16.3)");
    }
    struct Row { double z, N, fines, CR, N160, N160f, FS; };
    for (const Row w : {Row{3.00, 16, 5, 0.85, 20.3, 20.3, 0.90}, Row{9.00, 18, 22, 1.00, 17.3, 22.8, 0.83},
                        Row{10.50, 16, 32, 1.00, 14.3, 21.6, 0.79}}) {
        const lq::Result r = lq::evaluate(at(w.z, w.N, w.fines, w.CR), d);
        char label[96];
        std::snprintf(label, sizeof(label), "%.2f m: N1,60 %.1f, N1,60f %.1f, factor %.2f", w.z, w.N160, w.N160f, w.FS);
        std::printf("   %.2f m: N1,60 %.3f, N1,60f %.3f, tau_R %.3f, tau_eq %.3f, FS %.4f\n", w.z, r.N160,
                    r.N160f, r.tau_R, r.tau_eq, r.FS);
        check(near(r.N160, w.N160, 0.05) && near(r.N160f, w.N160f, 0.05) && near(r.FS, w.FS, 0.005) &&
                  r.verdict == lq::Verdict::Triggered,
              std::string(label) + ", liquefies");
    }
    {   // the screening
        const lq::Result dense = lq::evaluate(at(4.50, 39, 16, 0.85), d);
        std::printf("   4.50 m: N1,60 %.2f\n", dense.N160);
        check(near(dense.N160, 42.5, 0.05) && dense.verdict == lq::Verdict::NotRequired,
              "4.50 m: N1,60 = 42.5 >= 30, no evaluation required (16.6.5)");
        lq::SptTest clay = at(12.00, 20, 89, 1.00);
        clay.susceptible = false;   // ML, PI 16 %
        check(lq::evaluate(clay, d).verdict == lq::Verdict::NotSusceptible,
              "12.00 m: ML with PI 16 % is not a liquefiable soil (16.6.2)");
        check(lq::evaluate(at(1.50, 13, 17, 0.75), d).verdict == lq::Verdict::NotSusceptible,
              "1.50 m: above the water table, not liquefiable (16.6.2; the slides evaluate it anyway)");
        lq::SptTest deep = at(21.0, 10, 10, 1.00);
        check(lq::evaluate(deep, d).verdict == lq::Verdict::NotSusceptible, "21.0 m: deeper than 20 m (16.6.2)");
    }
    {   // the magnitude factor away from 7.5
        bool ok = true;
        for (const auto& [Mw, cm] : {std::pair{6.0, 1.77}, std::pair{6.5, 1.44}, std::pair{7.0, 1.19}}) {
            const double got = lq::magnitude_factor(Mw);
            std::printf("   C_M(%.1f) = %.4f (printed %.2f)\n", Mw, got, cm);
            ok = ok && near(got, cm, 0.005);
        }
        check(ok, "C_M = 10^2.24 / Mw^2.56 = 1.77 / 1.44 / 1.19 at Mw 6.0 / 6.5 / 7.0 (Eq. 16B.4c)");
    }
    {   // beyond the curve: N1,60 below 30 with enough fines to push N1,60f past it
        lq::SptTest t = at(8.0, 26, 30, 1.00);
        const lq::Result r = lq::evaluate(t, d);
        std::printf("   8.00 m, N 26, fines 30 %%: N1,60 %.2f, N1,60f %.2f\n", r.N160, r.N160f);
        check(r.N160 < 30.0 && r.N160f >= 30.0 && r.verdict == lq::Verdict::BeyondCurve && r.FS == 0.0,
              "N1,60 < 30 but N1,60f >= 30: beyond the curve, no factor invented");
    }
    if (g_failures == 0) {
        std::printf("\nOK: Annex 16B reproduces the authority's worked example\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", g_failures);
    return 1;
}
