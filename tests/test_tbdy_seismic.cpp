// TBDY 2018 seismic design spectrum + site coefficients + response spectrum (D4a).
// Verified against the OFFICIAL AFAD TBDY 2018 text (Sec. 2.3.4.1 Eq. 2.2-2.3, Tables 2.1/2.2) and,
// for the response spectrum, the closed-form SDOF steady-state amplification. Reference + locked
// values: docs/references/tbdy-2018-seismic.md.
#include <katai/analysis/response_spectrum.hpp>
#include <katai/analysis/dynamics.hpp>
#include <katai/math/sparse_matrix.hpp>

#include <Eigen/Core>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::SiteClass;
using katai::core::tbdy_site_coefficients;
using katai::core::tbdy_design_coefficients;
using katai::core::tbdy_elastic_spectrum;
using katai::core::response_spectrum;
using katai::core::Ec8GroundType;
using katai::core::Ec8SpectrumType;
using katai::core::ec8_spectrum_params;
using katai::core::ec8_eta;
using katai::core::ec8_elastic_spectrum;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

// (a) Horizontal elastic design spectrum: the four branches + corner continuity (TBDY Eq. 2.2-2.3).
void test_spectrum() {
    std::printf("-- (a) TBDY 2018 horizontal elastic design spectrum Sae(T) --\n");
    const double SDS = 1.2, SD1 = 0.4, TL = 6.0;
    const double TA = 0.2 * SD1 / SDS, TB = SD1 / SDS;  // 0.06667, 0.33333
    std::printf("   SDS=%.2f SD1=%.2f -> TA=%.4f TB=%.4f TL=%.1f\n", SDS, SD1, TA, TB, TL);
    // T=0 -> 0.4*SDS (spectrum ordinate at zero period).
    check(close(tbdy_elastic_spectrum(SDS, SD1, 0.0), 0.4 * SDS, 1e-12), "Sae(0) = 0.4 SDS");
    // Rising branch midpoint T=TA/2 -> (0.4+0.3)*SDS = 0.7 SDS.
    check(close(tbdy_elastic_spectrum(SDS, SD1, TA / 2), 0.7 * SDS, 1e-12), "rising branch (0.4+0.6 T/TA) SDS");
    // Plateau: TA <= T <= TB -> SDS.
    check(close(tbdy_elastic_spectrum(SDS, SD1, TA), SDS, 1e-12), "Sae(TA) = SDS (plateau start)");
    check(close(tbdy_elastic_spectrum(SDS, SD1, 0.2), SDS, 1e-12), "Sae plateau = SDS");
    check(close(tbdy_elastic_spectrum(SDS, SD1, TB), SDS, 1e-12), "Sae(TB) = SDS");
    // Corner continuity at TB: SD1/TB == SDS.
    check(close(SD1 / TB, SDS, 1e-12), "SD1/TB = SDS (continuity at TB)");
    // Constant-velocity branch TB < T <= TL -> SD1/T.
    check(close(tbdy_elastic_spectrum(SDS, SD1, 1.0), SD1 / 1.0, 1e-12), "Sae = SD1/T (constant velocity)");
    check(close(tbdy_elastic_spectrum(SDS, SD1, TL), SD1 / TL, 1e-12), "Sae(TL) = SD1/TL");
    // Constant-displacement branch T > TL -> SD1 TL / T^2; continuity at TL.
    check(close(tbdy_elastic_spectrum(SDS, SD1, 8.0), SD1 * TL / (8.0 * 8.0), 1e-12), "Sae = SD1 TL/T^2");
    check(close(SD1 / TL, SD1 * TL / (TL * TL), 1e-12), "continuity at TL");
}

// (b) Local soil effect coefficients FS, F1 (TBDY Tables 2.1/2.2 -- official values + interpolation).
void test_site_coefficients() {
    std::printf("-- (b) TBDY 2018 site coefficients FS, F1 (Tables 2.1/2.2) --\n");
    // Table anchors (exact).
    check(close(tbdy_site_coefficients(SiteClass::ZA, 0.5, 0.2).F_S, 0.8, 1e-12), "ZA FS = 0.8");
    check(close(tbdy_site_coefficients(SiteClass::ZC, 0.75, 0.3).F_S, 1.2, 1e-12), "ZC FS(SS=0.75) = 1.2");
    check(close(tbdy_site_coefficients(SiteClass::ZD, 1.00, 0.4).F_S, 1.1, 1e-12), "ZD FS(SS=1.00) = 1.1");
    check(close(tbdy_site_coefficients(SiteClass::ZE, 0.25, 0.1).F_S, 2.4, 1e-12), "ZE FS(SS=0.25) = 2.4");
    check(close(tbdy_site_coefficients(SiteClass::ZC, 0.5, 0.30).F_1, 1.5, 1e-12), "ZC F1(S1=0.30) = 1.5");
    check(close(tbdy_site_coefficients(SiteClass::ZE, 0.5, 0.10).F_1, 4.2, 1e-12), "ZE F1(S1=0.10) = 4.2");
    check(close(tbdy_site_coefficients(SiteClass::ZD, 0.5, 0.20).F_1, 2.2, 1e-12), "ZD F1(S1=0.20) = 2.2");
    // Linear interpolation: ZD FS at SS=0.875 (between 0.75->1.2 and 1.00->1.1) = 1.15.
    check(close(tbdy_site_coefficients(SiteClass::ZD, 0.875, 0.3).F_S, 1.15, 1e-9), "ZD FS interpolated (SS=0.875) = 1.15");
    // ZE F1 at S1=0.35 (between 0.30->2.8 and 0.40->2.4) = 2.6.
    check(close(tbdy_site_coefficients(SiteClass::ZE, 0.5, 0.35).F_1, 2.6, 1e-9), "ZE F1 interpolated (S1=0.35) = 2.6");
    // Clamp beyond the table: SS>=1.50 and SS<=0.25.
    check(close(tbdy_site_coefficients(SiteClass::ZE, 2.0, 0.1).F_S, 0.8, 1e-12), "ZE FS clamped (SS=2.0) = 0.8");
    check(close(tbdy_site_coefficients(SiteClass::ZD, 0.1, 0.05).F_S, 1.6, 1e-12), "ZD FS clamped (SS<0.25) = 1.6");
    // Design coefficients SDS = SS FS, SD1 = S1 F1 (worked example ZC, SS=1.0, S1=0.3).
    const auto d = tbdy_design_coefficients(SiteClass::ZC, 1.0, 0.3);
    check(close(d.F_S, 1.0 * 1.2, 1e-9) && close(d.F_1, 0.3 * 1.5, 1e-9), "SDS=SS FS=1.20, SD1=S1 F1=0.45 (ZC)");
    std::printf("   worked ZC (SS=1.0,S1=0.3): SDS=%.3f  SD1=%.3f\n", d.F_S, d.F_1);
}

// (c) Response spectrum from a cosine-tapered sine a_g = A w(t) sin(w0 t): the peak SDOF response
// converges to the closed-form 5%-damped steady state S_a(T) = A Rd(r,xi), r = w0/wn = T/T0. The
// ramp-up window suppresses the sudden-onset transient so max|u| == the steady amplitude. T->0 ->
// PGA = A; T=T0 (resonance) -> A/(2 xi). Verifies the SDOF Newmark reuse + the omega^2 max|u| formula.
void test_response_spectrum() {
    std::printf("-- (c) response spectrum from a tapered sine (S_a = A Rd(r,xi)) --\n");
    const double A = 1.0, f0 = 1.0, w0 = 2 * kPi * f0, xi = 0.05, T0 = 1.0 / f0;
    const double dt = T0 / 120.0, t_ramp = 30 * T0;
    const int ns = static_cast<int>(70 * T0 / dt);  // 30-cycle ramp + 40 cycles steady
    std::vector<double> accel(ns + 1);
    for (int i = 0; i <= ns; ++i) {
        const double t = i * dt, w = t < t_ramp ? 0.5 * (1 - std::cos(kPi * t / t_ramp)) : 1.0;
        accel[i] = A * w * std::sin(w0 * t);
    }
    const std::vector<double> periods = {0.0, 0.25, 0.5, 1.0, 2.0};
    const auto Sa = response_spectrum(accel, dt, periods, xi);
    auto Rd = [&](double r) { return 1.0 / std::sqrt((1 - r * r) * (1 - r * r) + (2 * xi * r) * (2 * xi * r)); };
    for (size_t k = 0; k < periods.size(); ++k) {
        const double T = periods[k], r = T / T0;
        const double exact = (T <= 0.0) ? A : A * Rd(r);
        std::printf("   T=%.3f (r=%.2f): S_a=%.4f  exact %s=%.4f  (err %.1f%%)\n",
                    T, r, Sa[k], T <= 0 ? "PGA" : "A*Rd", exact, exact > 0 ? 100 * (Sa[k] - exact) / exact : 0);
        const double tol = (T <= 0.0) ? 1e-9 : 0.05 * exact;
        check(close(Sa[k], exact, tol), "S_a(T) = A Rd(r,xi) (SDOF steady-state)");
    }
    // The response spectrum peaks at resonance T = T0.
    double peak = 0.0; size_t pk = 0;
    for (size_t k = 0; k < Sa.size(); ++k) if (Sa[k] > peak) { peak = Sa[k]; pk = k; }
    check(std::fabs(periods[pk] - T0) < 1e-9, "response spectrum peaks at the resonant period T0");
}

// (d) BROADBAND verification: a multi-sine record a_g = w(t) sum_j A_j sin(j w0 t) (realistic, not a
// single tone). The reference is ANALYTICAL: the steady SDOF response is the superposition of each
// tone's exact steady response, S_a(T) = wn^2 max_t |sum_j (A_j/wn^2) Rd(r_j) sin(j w0 t - phi_j)|,
// evaluated densely over the common period. This tests response_spectrum on broadband input against a
// closed-form reference (no numerical reference method -> no shared error).
void test_broadband() {
    std::printf("-- (d) response spectrum, broadband multi-sine vs analytical superposition --\n");
    const double xi = 0.05, f0 = 0.7, w0 = 2 * kPi * f0;
    const int J = 4; const double Aj[4] = {1.0, 0.7, 0.5, 0.35};
    const double dt = (1.0 / (J * f0)) / 40.0;       // resolve the highest tone (J*f0) well
    const double t_ramp = 25.0 / f0, tend = 60.0 / f0;
    const int ns = static_cast<int>(tend / dt);
    std::vector<double> accel(ns + 1);
    for (int i = 0; i <= ns; ++i) {
        const double t = i * dt, wn = t < t_ramp ? 0.5 * (1 - std::cos(kPi * t / t_ramp)) : 1.0;
        double s = 0.0; for (int j = 1; j <= J; ++j) s += Aj[j - 1] * std::sin(j * w0 * t);
        accel[i] = wn * s;
    }
    // Test periods (avoid landing exactly on a sharp tone resonance).
    const std::vector<double> periods = {0.30, 0.55, 0.9, 1.4, 2.1};
    const auto Sa = response_spectrum(accel, dt, periods, xi);
    const double P = 2 * kPi / w0;                    // common period of the superposition
    for (size_t k = 0; k < periods.size(); ++k) {
        const double T = periods[k], wn = 2 * kPi / T;
        double mx = 0.0;
        for (int m = 0; m <= 4000; ++m) {             // dense max over one common period (steady state)
            const double t = P * m / 4000.0; double u = 0.0;
            for (int j = 1; j <= J; ++j) {
                const double r = j * w0 / wn, Rd = 1.0 / std::sqrt((1 - r * r) * (1 - r * r) + (2 * xi * r) * (2 * xi * r));
                const double phi = std::atan2(2 * xi * r, 1 - r * r);
                u += (Aj[j - 1] / (wn * wn)) * Rd * std::sin(j * w0 * t - phi);
            }
            mx = std::max(mx, std::fabs(u));
        }
        const double ref = wn * wn * mx;
        std::printf("   T=%.2f: S_a=%.4f  analytic superpos=%.4f  (err %.2f%%)\n",
                    T, Sa[k], ref, 100 * (Sa[k] - ref) / ref);
        check(close(Sa[k], ref, 0.02 * ref), "broadband S_a = analytical steady-state superposition");
    }
}

// Naive single-step (no sub-stepping) SDOF Newmark response-spectrum ordinate, to demonstrate the
// short-period error that response_spectrum's sub-stepping removes.
double naive_Sa(const std::vector<double>& accel, double dt, double T, double xi) {
    const double w = 2 * kPi / T, k = w * w, c = 2 * xi * w;
    const double a0 = 4.0 / (dt * dt), a1 = 2.0 / dt, a2 = 4.0 / dt, keff = k + a0 + a1 * c;
    double u = 0, v = 0, a = -accel[0], umax = 0;
    for (size_t s = 0; s + 1 < accel.size(); ++s) {
        const double peff = -accel[s + 1] + (a0 * u + a2 * v + a) + c * (a1 * u + v);
        const double un = peff / keff, an = a0 * (un - u) - a2 * v - a;
        v += dt * 0.5 * (a + an); u = un; a = an; umax = std::max(umax, std::fabs(u));
    }
    return w * w * umax;
}

// (e) Short-period SILENT-ERROR guard: a coarse-dt record (few points per oscillator period). The
// analytical resonant ordinate is A/(2 xi). response_spectrum (sub-stepped) hits it; a naive single-
// step Newmark at the coarse dt is badly wrong (period elongation detunes the resonance). This is the
// exact "silently wrong short-period result" a naive tool would produce -- proven prevented here.
void test_short_period_guard() {
    std::printf("-- (e) short-period guard: coarse-dt record, sub-stepping vs naive Newmark --\n");
    const double A = 1.0, xi = 0.05, Ts = 0.2, fs = 1.0 / Ts, ws = 2 * kPi * fs;
    const double dt = 0.02;                            // 10 points per period -> coarse for T=0.2
    const double t_ramp = 25 * Ts, tend = 60 * Ts;
    const int ns = static_cast<int>(tend / dt);
    std::vector<double> accel(ns + 1);
    for (int i = 0; i <= ns; ++i) {
        const double t = i * dt, wn = t < t_ramp ? 0.5 * (1 - std::cos(kPi * t / t_ramp)) : 1.0;
        accel[i] = A * wn * std::sin(ws * t);
    }
    const double exact = A / (2 * xi);                 // resonant ordinate = 10 A
    const double sub = response_spectrum(accel, dt, {Ts}, xi)[0];
    const double naive = naive_Sa(accel, dt, Ts, xi);
    std::printf("   T=%.2f (dt=%.3f, %.0f pts/period): sub-stepped S_a=%.3f  naive S_a=%.3f  exact A/(2xi)=%.3f\n",
                Ts, dt, Ts / dt, sub, naive, exact);
    std::printf("     sub-stepped err=%.1f%%   naive err=%.1f%% (silent error a naive tool would give)\n",
                100 * (sub - exact) / exact, 100 * (naive - exact) / exact);
    // Sub-stepping fixes the INTEGRATION error (naive ~-19%); the residual (~-4%) is the INPUT-sampling
    // limit -- a 10-pts/period record cannot represent the 5 Hz tone exactly (physical, not a bug; real
    // records are sampled finer). So: sub-stepping cuts the short-period error several-fold, and the
    // naive single-step is badly wrong -- the exact silent error a naive tool would report.
    check(std::fabs(naive - exact) > 0.10 * exact, "naive single-step IS wrong at short period (guard justified)");
    check(std::fabs(sub - exact) < 0.4 * std::fabs(naive - exact), "sub-stepping cuts the short-period error >2x");
    check(std::fabs(sub - exact) < 0.06 * exact, "sub-stepped ordinate accurate to the input-sampling limit");
    // With an adequately-sampled record (dt <= T/20) the ordinate is accurate to <1% (typical real records).
    std::vector<double> fine(4 * ns + 1);
    for (size_t i = 0; i < fine.size(); ++i) {
        const double t = i * (dt / 4), wf = t < t_ramp ? 0.5 * (1 - std::cos(kPi * t / t_ramp)) : 1.0;
        fine[i] = A * wf * std::sin(ws * t);
    }
    const double sub_fine = response_spectrum(fine, dt / 4, {Ts}, xi)[0];
    std::printf("     fine record (dt=%.4f, %.0f pts/period): S_a=%.3f  err=%.1f%%\n",
                dt / 4, Ts / (dt / 4), sub_fine, 100 * (sub_fine - exact) / exact);
    check(std::fabs(sub_fine - exact) < 0.01 * exact, "adequately-sampled record: short-period ordinate <1%");
}

// (f) Engine link: the inline SDOF Newmark inside response_spectrum matches the general solve_newmark
// (dynamics.hpp) to round-off on the same 1-DOF system -> the response spectrum uses the SAME verified
// scheme as the 2D dynamic engine.
void test_engine_link() {
    std::printf("-- (f) response_spectrum inline Newmark == solve_newmark (engine) --\n");
    const double A = 1.0, xi = 0.05, T = 0.7, w = 2 * kPi / T, c = 2 * xi * w;
    const double f_drive = 1.0, wd = 2 * kPi * f_drive, dt = T / 200.0;  // fine dt, sub=1 -> same steps
    const int ns = 3000;
    std::vector<double> accel(ns + 1);
    for (int i = 0; i <= ns; ++i) accel[i] = A * std::sin(wd * i * dt);
    const double sa_rs = response_spectrum(accel, dt, {T}, xi)[0];
    // solve_newmark on the identical SDOF (m=1, k=w^2, c) with force -a_g.
    katai::math::SparseMatrixBuilder bM(1), bC(1), bK(1);
    bM.add_entry(0, 0, 1.0); bC.add_entry(0, 0, c); bK.add_entry(0, 0, w * w);
    auto force = [&](int s) { Eigen::VectorXd f(1); f << -accel[std::min(s, ns)]; return f; };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(1);
    const auto R = katai::core::solve_newmark(bM.build(), bC.build(), bK.build(), force, dt, ns, z, z);
    double umax = 0.0; for (const auto& u : R.u) umax = std::max(umax, std::fabs(u[0]));
    const double sa_nm = w * w * umax;
    std::printf("   S_a(response_spectrum)=%.6f  S_a(solve_newmark)=%.6f  (rel diff %.2e)\n",
                sa_rs, sa_nm, std::fabs(sa_rs - sa_nm) / sa_nm);
    check(close(sa_rs, sa_nm, 1e-6 * sa_nm), "response_spectrum inline Newmark == solve_newmark engine");
}

}  // namespace

// (f) EC8 (EN 1998-1:2004) horizontal elastic response spectrum -- sec 3.2.2.2 Eq. (3.2)-(3.5),
// Tables 3.2/3.3 (EN recommended values), damping correction Eq. (3.6). The parameter tables are
// RETYPED here from the verified source (independent transcription: a typo in either copy fails).
void test_ec8_spectrum() {
    std::printf("-- (f) EC8 (EN 1998-1) horizontal elastic response spectrum Se(T) --\n");
    // Independent retype of Tables 3.2 (Type 1) and 3.3 (Type 2): {S, TB, TC, TD} for A..E.
    const double t1[5][4] = {{1.00, 0.15, 0.40, 2.0}, {1.20, 0.15, 0.50, 2.0},
                             {1.15, 0.20, 0.60, 2.0}, {1.35, 0.20, 0.80, 2.0},
                             {1.40, 0.15, 0.50, 2.0}};
    const double t2[5][4] = {{1.00, 0.05, 0.25, 1.2}, {1.35, 0.05, 0.25, 1.2},
                             {1.50, 0.10, 0.25, 1.2}, {1.80, 0.10, 0.30, 1.2},
                             {1.60, 0.05, 0.25, 1.2}};
    bool tables_ok = true;
    for (int g = 0; g < 5; ++g) {
        const auto p1 = ec8_spectrum_params((Ec8GroundType)g, Ec8SpectrumType::Type1);
        const auto p2 = ec8_spectrum_params((Ec8GroundType)g, Ec8SpectrumType::Type2);
        tables_ok = tables_ok && p1.S == t1[g][0] && p1.T_B == t1[g][1] && p1.T_C == t1[g][2] &&
                    p1.T_D == t1[g][3] && p2.S == t2[g][0] && p2.T_B == t2[g][1] &&
                    p2.T_C == t2[g][2] && p2.T_D == t2[g][3];
    }
    check(tables_ok, "Tables 3.2/3.3 match the independent retype, all 10 rows exactly");

    // Damping correction (Eq. 3.6): eta(5%) = 1 exactly; eta(10%) = sqrt(10/15); floor 0.55.
    check(close(ec8_eta(5.0), 1.0, 1e-15), "eta(5%) = 1 (the reference damping)");
    check(close(ec8_eta(10.0), std::sqrt(10.0 / 15.0), 1e-15), "eta(10%) = sqrt(10/15)");
    check(close(ec8_eta(50.0), 0.55, 1e-15), "eta floors at 0.55 (Eq. 3.6 lower bound)");

    // Branches + corners on ground C, Type 1 (S=1.15, TB=0.2, TC=0.6, TD=2.0), a_g = 3.0 m/s^2.
    const double ag = 3.0, S = 1.15, TB = 0.20, TC = 0.60, TD = 2.0;
    const auto gC = Ec8GroundType::C; const auto T1 = Ec8SpectrumType::Type1;
    check(close(ec8_elastic_spectrum(ag, gC, T1, 0.0), ag * S, 1e-12), "Se(0) = a_g S (Eq. 3.2 at T=0)");
    check(close(ec8_elastic_spectrum(ag, gC, T1, 0.5 * (TB + TC)), ag * S * 2.5, 1e-12),
          "plateau Se = 2.5 a_g S over TB..TC (Eq. 3.3, eta=1)");
    check(close(ec8_elastic_spectrum(ag, gC, T1, 1.0), ag * S * 2.5 * TC / 1.0, 1e-12),
          "1/T branch: Se(1.0 s) = 2.5 a_g S TC/T (Eq. 3.4)");
    check(close(ec8_elastic_spectrum(ag, gC, T1, 3.0), ag * S * 2.5 * TC * TD / 9.0, 1e-12),
          "1/T^2 branch: Se(3.0 s) = 2.5 a_g S TC TD/T^2 (Eq. 3.5)");
    // Corner continuity: the branch pairs must agree AT TB, TC, TD (left/right evaluation).
    for (double Tc : {TB, TC, TD}) {
        const double lo = ec8_elastic_spectrum(ag, gC, T1, Tc - 1e-9);
        const double hi = ec8_elastic_spectrum(ag, gC, T1, Tc + 1e-9);
        check(std::fabs(lo - hi) < 1e-6 * hi, "branch continuity at a corner period");
    }
    // Linear scaling in a_g (= gamma_I a_gR): doubling a_g doubles every ordinate.
    check(close(ec8_elastic_spectrum(2 * ag, gC, T1, 0.3), 2 * ec8_elastic_spectrum(ag, gC, T1, 0.3), 1e-12),
          "Se scales linearly with a_g = gamma_I a_gR");
    // Type 2 differs materially from Type 1 on the same ground (D: S 1.35 vs 1.80, TC 0.8 vs 0.3).
    const double d1 = ec8_elastic_spectrum(ag, Ec8GroundType::D, Ec8SpectrumType::Type1, 0.5);
    const double d2 = ec8_elastic_spectrum(ag, Ec8GroundType::D, Ec8SpectrumType::Type2, 0.5);
    std::printf("   ground D at T=0.5 s: Type1 Se = %.4f, Type2 Se = %.4f m/s^2\n", d1, d2);
    check(std::fabs(d1 - d2) > 0.1, "Type 1 and Type 2 are genuinely different spectra");
}

int main() {
    std::printf("TBDY 2018 seismic design spectrum + response spectrum (D4a)\n\n");
    test_spectrum();
    test_site_coefficients();
    test_response_spectrum();
    test_broadband();
    test_short_period_guard();
    test_engine_link();
    std::printf("\n");
    test_ec8_spectrum();
    if (g_failures == 0) {
        std::printf("\nOK: TBDY 2018 spectrum + site coefficients + response spectrum verified\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
