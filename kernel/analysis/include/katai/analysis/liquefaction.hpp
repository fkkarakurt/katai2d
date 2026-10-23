#pragma once
// LIQUEFACTION TRIGGERING BY THE SIMPLIFIED METHOD (TBDY 2018, Annex 16B -- mandatory annex of the
// Turkish Building Earthquake Code, Resmi Gazete 30364 bis, 18.03.2018, printed pp. 362-364).
//
// One SPT test at a time: the field blow count corrected to N1,60 (Eq. 16B.1-16B.2, Table 16B.1)
// and for its fines content to N1,60f (Eq. 16B.3), the cyclic resistance ratio of a magnitude 7.5
// earthquake (Eq. 16B.4b) scaled to the design magnitude (Eq. 16B.4c), and the cyclic shear stress
// of the simplified method (Eq. 16B.5, 16B.6). The ground is safe against triggering where
// tau_R / tau_eq >= 1.10 (Eq. 16.3).
//
// What the code does not say is not filled in here; it is said. The CRR curve of Eq. 16B.4b is the
// clean-sand curve of Youd et al. (2001) and is undefined at N1,60f = 34 (its first term); the
// annex gives no upper limit. Youd et al. state the curve for N1,60cs < 30 and class clean
// granular soil denser than that as too dense to liquefy; 16.6.5 asks for an evaluation only
// where N1,60 < 30. A test whose N1,60f reaches 30 is therefore reported as beyond the curve,
// with both statements, and no factor of safety is invented for it.

#include <algorithm>
#include <cmath>
#include <string>

namespace katai::core::liquefaction {

// Table 16B.1 -- the corrections other than depth. CR is read from the rod length; the table
// starts at 3 m, and a shorter rod is flagged rather than extrapolated.
struct SptCorrections {
    double CR = 1.0;   // rod length
    double CS = 1.0;   // sampler: 1.00 standard (with liner), 1.10-1.30 without
    double CB = 1.0;   // borehole diameter: 1.00 for 65-115 mm, 1.05 for 150 mm, 1.15 for 200 mm
    double CE = 1.0;   // energy ratio: safety hammer 0.60-1.17, donut 0.45-1.00, automatic 0.90-1.60
};
inline double rod_length_factor(double rod_length_m, bool& in_table) {
    in_table = rod_length_m >= 3.0;
    if (rod_length_m < 4.0) return 0.75;
    if (rod_length_m < 6.0) return 0.85;
    if (rod_length_m < 10.0) return 0.95;
    return 1.00;
}

struct SptTest {
    double depth = 0.0;         // z below the ground surface [m]
    double N = 0.0;             // field blow count [blows / 30 cm]
    double fines = 0.0;         // fines content IDI [%]
    double sigma_v = 0.0;       // total vertical stress at the test depth, as the site will be [kPa]
    double sigma_v_eff = 0.0;   // effective vertical stress at the test depth [kPa]
    // Eq. 16B.2 reads the effective stress "as the site was when the test was made"; fills,
    // footings or excavations after it are not taken into account. 0 = the same as sigma_v_eff.
    double sigma_v_eff_at_test = 0.0;
    SptCorrections c;
    bool below_water = true;    // 16.6.2 / 16.6.4: only ground below the water table liquefies
    bool susceptible = true;    // the engineer's classification: cohesionless or PI < 12 % (16.6.2,
                                // 16.6.4) and not excluded by 16.6.6
};

struct Design {
    double SDS = 0.0;           // short-period design spectral acceleration coefficient
    double Mw = 7.5;            // design earthquake moment magnitude
    double required = 1.10;     // Eq. 16.3
};

enum class Verdict {
    Safe,            // tau_R / tau_eq >= required
    Triggered,       // below it
    NotSusceptible,  // above the water table, deeper than 20 m, or not a liquefiable soil
    NotRequired,     // N1,60 >= 30 (16.6.5)
    BeyondCurve,     // N1,60f >= 30: Eq. 16B.4b gives no defined resistance there
};

struct Result {
    double CN = 0.0, N160 = 0.0, alpha = 0.0, beta = 1.0, N160f = 0.0;
    double CRR75 = 0.0, CM = 0.0, tau_R = 0.0, rd = 0.0, tau_eq = 0.0, FS = 0.0;
    Verdict verdict = Verdict::NotSusceptible;
    std::string note;
};

// Eq. 16B.2 (sigma'_vo in kN/m2).
inline double depth_factor(double sigma_v_eff) {
    return std::min(9.78 * std::sqrt(1.0 / sigma_v_eff), 1.70);
}
// Eq. 16B.3b.
inline void fines_terms(double IDI, double& alpha, double& beta) {
    if (IDI <= 5.0) { alpha = 0.0; beta = 1.0; return; }
    if (IDI >= 35.0) { alpha = 5.0; beta = 1.2; return; }
    alpha = std::exp(1.76 - 190.0 / (IDI * IDI));
    beta = 0.99 + std::pow(IDI, 1.5) / 1000.0;
}
// Eq. 16B.4b.
inline double crr75(double N160f) {
    const double d = 10.0 * N160f + 45.0;
    return 1.0 / (34.0 - N160f) + N160f / 135.0 + 50.0 / (d * d) - 1.0 / 200.0;
}
// Eq. 16B.4c.
inline double magnitude_factor(double Mw) { return std::pow(10.0, 2.24) / std::pow(Mw, 2.56); }
// Eq. 16B.6.
inline double stress_reduction(double z) {
    if (z <= 9.15) return 1.0 - 0.00765 * z;
    if (z <= 23.0) return 1.174 - 0.0267 * z;
    if (z <= 30.0) return 0.744 - 0.008 * z;
    return 0.50;
}

inline Result evaluate(const SptTest& t, const Design& d) {
    Result r;
    if (!t.below_water || t.depth > 20.0 || !t.susceptible) {
        r.verdict = Verdict::NotSusceptible;
        r.note = !t.below_water ? "above the water table (16.6.2)"
                 : t.depth > 20.0 ? "deeper than 20 m (16.6.2)"
                                  : "not a liquefiable soil by the engineer's classification (16.6.4)";
        return r;
    }
    const double s_test = t.sigma_v_eff_at_test > 0.0 ? t.sigma_v_eff_at_test : t.sigma_v_eff;
    r.CN = depth_factor(s_test);
    r.N160 = t.N * r.CN * t.c.CR * t.c.CS * t.c.CB * t.c.CE;
    if (r.N160 >= 30.0) {
        r.verdict = Verdict::NotRequired;
        r.note = "N1,60 >= 30: no triggering evaluation is required (16.6.5)";
        return r;
    }
    fines_terms(t.fines, r.alpha, r.beta);
    r.N160f = r.alpha + r.beta * r.N160;
    r.rd = stress_reduction(t.depth);
    r.tau_eq = 0.65 * t.sigma_v * (0.4 * d.SDS) * r.rd;
    if (r.N160f >= 30.0) {
        r.verdict = Verdict::BeyondCurve;
        r.note = "N1,60f >= 30: Eq. 16B.4b is the clean-sand curve stated for N1,60cs < 30 "
                 "(Youd et al. 2001), undefined at 34, and the annex gives no value beyond it";
        return r;
    }
    r.CRR75 = crr75(r.N160f);
    r.CM = magnitude_factor(d.Mw);
    r.tau_R = r.CRR75 * r.CM * t.sigma_v_eff;
    r.FS = r.tau_eq > 0.0 ? r.tau_R / r.tau_eq : 1e300;
    r.verdict = r.FS >= d.required ? Verdict::Safe : Verdict::Triggered;
    return r;
}

}  // namespace katai::core::liquefaction
