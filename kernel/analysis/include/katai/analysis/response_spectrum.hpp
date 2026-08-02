#pragma once
// TBDY 2018 seismic design spectrum + response spectrum from an acceleration record.
// Turkish Building Earthquake Code (TBDY) 2018 §2.3 (horizontal elastic design spectrum +
// local site coefficients Tables 2.1/2.2). The response spectrum = reuse of the dynamics
// engine (analysis/dynamics.hpp solve_newmark) on an SDOF oscillator. Official values +
// formulation: docs/references/tbdy-2018-seismic.md (AFAD TBDY 2018).

#include <algorithm>
#include <cmath>
#include <vector>

namespace katai::core {

// TBDY 2018 local site class (Table 16.1). ZF = site-specific analysis (no table → coefficients undefined).
enum class SiteClass { ZA, ZB, ZC, ZD, ZE };

struct SiteCoefficients { double F_S = 1.0, F_1 = 1.0; };

namespace detail {
// Linear interpolation; x outside the anchor range clamps to the end value. n anchor points.
inline double interp_clamp(const double* xs, const double* ys, int n, double x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    int i = 0;
    while (i < n - 1 && x > xs[i + 1]) ++i;
    const double t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    return ys[i] + t * (ys[i + 1] - ys[i]);
}
}  // namespace detail

// Local site coefficients F_S, F_1 (TBDY 2018 Tables 2.1/2.2 — OFFICIAL AFAD values; 6
// columns). S_S = short-period, S_1 = 1.0 s map spectral acceleration coefficient. Linear
// interpolation in between.
inline SiteCoefficients tbdy_site_coefficients(SiteClass cls, double S_S, double S_1) {
    static const double ss[6] = {0.25, 0.50, 0.75, 1.00, 1.25, 1.50};
    static const double s1[6] = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60};
    static const double FS[5][6] = {
        {0.8, 0.8, 0.8, 0.8, 0.8, 0.8},  // ZA
        {0.9, 0.9, 0.9, 0.9, 0.9, 0.9},  // ZB
        {1.3, 1.3, 1.2, 1.2, 1.2, 1.2},  // ZC
        {1.6, 1.4, 1.2, 1.1, 1.0, 1.0},  // ZD
        {2.4, 1.7, 1.3, 1.1, 0.9, 0.8},  // ZE
    };
    static const double F1[5][6] = {
        {0.8, 0.8, 0.8, 0.8, 0.8, 0.8},  // ZA
        {0.8, 0.8, 0.8, 0.8, 0.8, 0.8},  // ZB
        {1.5, 1.5, 1.5, 1.5, 1.5, 1.4},  // ZC
        {2.4, 2.2, 2.0, 1.9, 1.8, 1.7},  // ZD
        {4.2, 3.3, 2.8, 2.4, 2.2, 2.0},  // ZE
    };
    const int c = static_cast<int>(cls);
    return {detail::interp_clamp(ss, FS[c], 6, S_S), detail::interp_clamp(s1, F1[c], 6, S_1)};
}

// Design spectral acceleration coefficients S_DS = S_S·F_S, S_D1 = S_1·F_1 (TBDY §2.3.2.2).
inline SiteCoefficients tbdy_design_coefficients(SiteClass cls, double S_S, double S_1) {
    const auto f = tbdy_site_coefficients(cls, S_S, S_1);
    return {S_S * f.F_S, S_1 * f.F_1};  // {S_DS, S_D1}
}

// TBDY 2018 horizontal elastic design spectrum S_ae(T) [in g] (§2.3.4.1, Eqs. 2.2-2.3).
// Corner periods T_A=0.2·S_D1/S_DS, T_B=S_D1/S_DS; constant-displacement transition T_L=6 s.
inline double tbdy_elastic_spectrum(double S_DS, double S_D1, double T, double T_L = 6.0) {
    if (S_DS <= 0.0) return 0.0;
    const double T_A = 0.2 * S_D1 / S_DS;
    const double T_B = S_D1 / S_DS;
    if (T <= T_A) return (0.4 + 0.6 * T / T_A) * S_DS;
    if (T <= T_B) return S_DS;
    if (T <= T_L) return S_D1 / T;
    return S_D1 * T_L / (T * T);
}

// TBDY 2018 horizontal elastic design DISPLACEMENT spectrum S_de(T) [m] (§2.3.4.2, Eq. 2.4):
// S_de = (T²/4π²)·g·S_ae. g = gravitational acceleration [m/s²].
inline double tbdy_displacement_spectrum(double S_DS, double S_D1, double T, double g = 9.81,
                                         double T_L = 6.0) {
    constexpr double kPi = 3.14159265358979323846;
    return (T * T / (4.0 * kPi * kPi)) * g * tbdy_elastic_spectrum(S_DS, S_D1, T, T_L);
}

// ---- EC8 (EN 1998-1:2004) horizontal elastic response spectrum --------------------------------
// §3.2.2.2, Eqs. (3.2)-(3.5) + damping correction η (Eq. 3.6). Ground types A-E (Table
// 3.1); recommended S, T_B, T_C, T_D values for Type 1 (M_s > 5.5, Table 3.2) and Type 2
// (M_s ≤ 5.5, Table 3.3). The table values were cross-checked against independent
// engineering documentation (see the EC8 appendix of docs/references/tbdy-2018-seismic.md).
// National annexes may change these values — the ones here are the EN RECOMMENDED values
// (declared in the report).
enum class Ec8GroundType { A, B, C, D, E };
enum class Ec8SpectrumType { Type1, Type2 };

struct Ec8Params { double S, T_B, T_C, T_D; };

inline Ec8Params ec8_spectrum_params(Ec8GroundType ground, Ec8SpectrumType type) {
    // Table 3.2 (Type 1) and Table 3.3 (Type 2), columns S, T_B, T_C, T_D; rows A..E.
    static const Ec8Params kType1[5] = {
        {1.00, 0.15, 0.40, 2.0},   // A
        {1.20, 0.15, 0.50, 2.0},   // B
        {1.15, 0.20, 0.60, 2.0},   // C
        {1.35, 0.20, 0.80, 2.0},   // D
        {1.40, 0.15, 0.50, 2.0},   // E
    };
    static const Ec8Params kType2[5] = {
        {1.00, 0.05, 0.25, 1.2},   // A
        {1.35, 0.05, 0.25, 1.2},   // B
        {1.50, 0.10, 0.25, 1.2},   // C
        {1.80, 0.10, 0.30, 1.2},   // D
        {1.60, 0.05, 0.25, 1.2},   // E
    };
    const int g = static_cast<int>(ground);
    return type == Ec8SpectrumType::Type1 ? kType1[g] : kType2[g];
}

// Damping correction factor η = sqrt(10/(5+ξ)) ≥ 0.55 (Eq. 3.6; ξ in PERCENT, η = 1 at 5%).
inline double ec8_eta(double xi_percent) {
    return std::max(std::sqrt(10.0 / (5.0 + xi_percent)), 0.55);
}

// S_e(T), in the SAME unit as a_g (a_g = γ_I·a_gR — the type-A-ground design acceleration).
// Eqs. (3.2)-(3.5):
//   0 ≤ T ≤ T_B : a_g·S·[1 + (T/T_B)·(2.5η − 1)]
//   T_B ≤ T ≤ T_C: a_g·S·2.5η
//   T_C ≤ T ≤ T_D: a_g·S·2.5η·(T_C/T)
//   T_D ≤ T      : a_g·S·2.5η·(T_C·T_D/T²)    (EN defines up to 4 s; the tail continues the same formula)
inline double ec8_elastic_spectrum(double a_g, Ec8GroundType ground, Ec8SpectrumType type,
                                   double T, double xi_percent = 5.0) {
    if (a_g <= 0.0) return 0.0;
    const Ec8Params p = ec8_spectrum_params(ground, type);
    const double eta = ec8_eta(xi_percent);
    if (T <= 0.0) return a_g * p.S;
    if (T <= p.T_B) return a_g * p.S * (1.0 + (T / p.T_B) * (2.5 * eta - 1.0));
    if (T <= p.T_C) return a_g * p.S * 2.5 * eta;
    if (T <= p.T_D) return a_g * p.S * 2.5 * eta * (p.T_C / T);
    return a_g * p.S * 2.5 * eta * (p.T_C * p.T_D / (T * T));
}

// Elastic pseudo-response spectrum at ξ% damping: for each period T a single-DOF oscillator
// (m=1, ω=2π/T, c=2ξω) is driven by the base excitation −a_g(t) with Newmark-β (γ=½,β=¼);
// pseudo-spectral acceleration S_a(T)=ω²·max|u_rel|. `accel` = the acceleration record
// (uniform dt); `periods` = the requested T (T≤0 → rigid = follows the ground = PGA). The
// SCHEME is the SAME as `dynamics.hpp` `solve_newmark` (cross-checked exactly in the test),
// specialized to SDOF: no history stored (only max|u|). **SUB-STEPPING:** interpolates the
// record so each oscillator gets ≥25 steps/period → removes the silent Newmark
// period-elongation error at SHORT periods (correct even on a coarse-dt record). Reference:
// docs/validation/seismic-verification.md.
inline std::vector<double> response_spectrum(const std::vector<double>& accel, double dt,
                                             const std::vector<double>& periods, double xi = 0.05) {
    constexpr double kPi = 3.14159265358979323846;
    const int n = static_cast<int>(accel.size());
    std::vector<double> Sa(periods.size(), 0.0);
    if (n < 2 || dt <= 0.0) return Sa;
    double pga = 0.0;
    for (double a : accel) pga = std::max(pga, std::fabs(a));
    for (size_t idx = 0; idx < periods.size(); ++idx) {
        const double T = periods[idx];
        if (T <= 0.0) { Sa[idx] = pga; continue; }           // rigid oscillator follows the ground
        const int sub = std::max(1, static_cast<int>(std::ceil(25.0 * dt / T)));  // >=25 steps/period
        const double h = dt / sub;
        const double w = 2.0 * kPi / T, k = w * w, c = 2.0 * xi * w;  // m = 1
        // Newmark-β constants (γ=½, β=¼): a0=1/(βh²)=4/h², a1=γ/(βh)=2/h, a2=1/(βh)=4/h; a3=a4=1, a5=0.
        const double a0 = 4.0 / (h * h), a1 = 2.0 / h, a2 = 4.0 / h;
        const double keff = k + a0 + a1 * c;
        double u = 0.0, v = 0.0, a = -accel[0], umax = 0.0;   // a(0) = (p0 - c v0 - k u0)/m = -a_g(0)
        const long total = static_cast<long>(n - 1) * sub;
        for (long s = 0; s < total; ++s) {
            const double t = (s + 1) * h, x = t / dt;         // interpolate a_g at the next step
            int i = static_cast<int>(x); if (i > n - 2) i = n - 2;
            const double f = x - i, ag = accel[i] * (1.0 - f) + accel[i + 1] * f;
            const double peff = -ag + (a0 * u + a2 * v + a) + c * (a1 * u + v);
            const double un = peff / keff;
            const double an = a0 * (un - u) - a2 * v - a;
            v = v + h * 0.5 * (a + an);
            u = un; a = an;
            umax = std::max(umax, std::fabs(u));
        }
        Sa[idx] = w * w * umax;  // pseudo-spectral acceleration
    }
    return Sa;
}

}  // namespace katai::core
