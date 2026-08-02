// BENCHMARK Wave-1 (ROADMAP Track 7): MULTI-LAYER 1D site response against the SHAKE-class
// layered viscoelastic transfer-matrix solution -- the first INDEPENDENT dynamic benchmark.
//
// THE ORACLE (shares no code with the solver): the exact steady-state harmonic solution of a
// layered shear column on a rigid base driven by within motion a_g = A e^{iwt}. Per layer,
// relative-frame equilibrium with the EXACT continuum analog of the FE's Rayleigh damping
// (M u'' + (aM + bK) u' + K u = -M r a_g):
//
//     -w^2 rho* u - G* u'' = -rho A,   G* = G (1 + i w beta),   rho* = rho (1 - i alpha / w)
//
// so the comparison is like-for-like: no "hysteretic vs viscous" approximation gap -- the oracle
// damps EXACTLY as the FE does (SHAKE's own frequency-independent hysteretic damping differs from
// Rayleigh by construction; that difference is a MODELING choice, documented in the validation
// note, not an error band we hide behind). General solution per layer j (local z from the layer
// bottom): u_j = a_j cos(k_j z) + b_j sin(k_j z) + P_j with k_j = w sqrt(rho*_j / G*_j) and the
// particular constant P_j = A rho_j / (w^2 rho*_j). Boundary/interface conditions:
//   base:      u = 0                          (rigid base, relative frame)
//   interface: u and tau = G* u' continuous
//   surface:   tau = 0
// Transfer function (total surface motion / ground motion): T(w) = 1 - w^2 u(surface) / A.
// Closed-form check built into the test: single undamped layer collapses to the classic
// T = 1 / cos(wH/Vs) (resonances at f = (2n-1) Vs / 4H) -- pinning the oracle's own signs/BCs.
//
// THE MEASUREMENT: the GUI compute path (mesh_from_project + solve_phases, PhaseType::Dynamic)
// solves the same column in 2D (SH boundary conditions), and |T|_FE = steady-state amplitude of
// the reported total surface acceleration / A. Discretization (mesh + Newmark dt + leftover
// transient) is the only legitimate gap; tolerance bands are stated per point and kept honest.
//
// Independent Rayleigh cross-check: alpha/beta are recomputed here from the two-frequency
// closed form (alpha = 2 xi w1 w2/(w1+w2), beta = 2 xi/(w1+w2)); if the core mapped the phase's
// (f1, f2, xi) input to different coefficients, every damped comparison below would fail.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// One soil layer of the 1D column, listed TOP-DOWN (surface layer first).
struct Layer {
    double h;    // thickness [m]
    double rho;  // mass density [Mg/m^3]
    double G;    // shear modulus [kPa]
};

// |T(w)|: layered viscoelastic transfer function (see file header). Independent implementation:
// small dense complex system (2 unknowns per layer), Eigen partial-pivot LU.
double layered_transfer(const std::vector<Layer>& top_down, double alpha, double beta, double w) {
    using C = std::complex<double>;
    const int n = (int)top_down.size();
    std::vector<Layer> up(top_down.rbegin(), top_down.rend());  // bottom-up (layer 0 sits on the base)
    std::vector<C> k(n), Gs(n), P(n);
    const C i01(0.0, 1.0);
    const double A = 1.0;  // T is amplitude-independent; use unit ground acceleration
    for (int j = 0; j < n; ++j) {
        const C Gstar = up[j].G * (1.0 + i01 * (w * beta));
        const C rstar = up[j].rho * (1.0 - i01 * (alpha / w));
        Gs[j] = Gstar;
        k[j] = w * std::sqrt(rstar / Gstar);
        P[j] = A * up[j].rho / (w * w * rstar);
    }
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
    Eigen::VectorXcd r = Eigen::VectorXcd::Zero(2 * n);
    int eq = 0;
    S(eq, 0) = 1.0; r(eq) = -P[0]; ++eq;                       // base: a_0 + P_0 = 0
    for (int j = 0; j + 1 < n; ++j) {                          // interfaces
        const C ch = std::cos(k[j] * up[j].h), sh = std::sin(k[j] * up[j].h);
        S(eq, 2 * j) = ch; S(eq, 2 * j + 1) = sh;              // u continuity
        S(eq, 2 * (j + 1)) = -1.0;
        r(eq) = P[j + 1] - P[j]; ++eq;
        S(eq, 2 * j) = -Gs[j] * k[j] * sh;                     // tau continuity
        S(eq, 2 * j + 1) = Gs[j] * k[j] * ch;
        S(eq, 2 * (j + 1) + 1) = -Gs[j + 1] * k[j + 1];
        r(eq) = 0.0; ++eq;
    }
    {                                                          // surface: tau = 0
        const int t = n - 1;
        const C ch = std::cos(k[t] * up[t].h), sh = std::sin(k[t] * up[t].h);
        S(eq, 2 * t) = -sh; S(eq, 2 * t + 1) = ch;
        r(eq) = 0.0; ++eq;
    }
    const Eigen::VectorXcd x = S.partialPivLu().solve(r);
    const int t = n - 1;
    const C u_surf = x(2 * t) * std::cos(k[t] * up[t].h) + x(2 * t + 1) * std::sin(k[t] * up[t].h) + P[t];
    return std::abs(1.0 - w * w * u_surf / A);
}

// Independent two-frequency Rayleigh closed form (equal xi at w1, w2).
void rayleigh_ab(double f1, double f2, double xi, double& alpha, double& beta) {
    const double w1 = 2 * kPi * f1, w2 = 2 * kPi * f2;
    alpha = 2.0 * xi * w1 * w2 / (w1 + w2);
    beta = 2.0 * xi / (w1 + w2);
}

// Steady-state |T| from the FE run: max |total surface acceleration| over the LAST quarter of the
// record (the transient has decayed by then; see the per-run cycle counts), divided by A.
double fe_transfer(const katai::app::SolveResult& R, double A) {
    if (R.dyn_time.empty()) return -1.0;
    const double t_end = R.dyn_time.back(), t_from = 0.75 * t_end;
    double amp = 0.0;
    for (size_t i = 0; i < R.dyn_time.size(); ++i)
        if (R.dyn_time[i] >= t_from) amp = std::fmax(amp, std::fabs(R.dyn_surface_ax[i]));
    return amp / A;
}

// Build and solve the column (single- or two-layer) at one driving frequency through the REAL
// user path (mesh_from_project + solve_phases). Returns |T|_FE (or -1 on failure).
struct ColumnSpec {
    std::vector<Layer> top_down;   // 1 or 2 layers
    double nu = 0.3;
    double ray_f1 = 1.0, ray_f2 = 3.0, xi = 0.05;
};
double run_column(const ColumnSpec& cs, double freq, std::string* msg = nullptr,
                  katai::app::SolveResult* out = nullptr) {
    constexpr double kW = 2.0, kGrav = 9.81, kA = 1.0;
    m::Project pr;
    double ytop = 0.0;
    for (const auto& L : cs.top_down) ytop += L.h;
    // Materials, top-down order (material i = layer i).
    for (const auto& L : cs.top_down) {
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 2.0 * (1.0 + cs.nu) * L.G;    // E from the layer's G
        s.nu = cs.nu;
        s.gamma_unsat = L.rho * kGrav; s.gamma_sat = s.gamma_unsat; s.e_init = 0.5;
        pr.materials.push_back(s);
    }
    // Polygons, top-down; SH column BCs: base FullyFixed, sides VerticallyFixed (pure shear),
    // internal interfaces Free.
    double y1 = ytop;
    for (size_t li = 0; li < cs.top_down.size(); ++li) {
        const double y0 = y1 - cs.top_down[li].h;
        m::SoilPolygon P; P.material = (int)li;
        P.x = {0, kW, kW, 0}; P.y = {y0, y0, y1, y1};
        const bool base = li + 1 == cs.top_down.size();
        P.edge_bc = {base ? (int)m::BCType::FullyFixed : (int)m::BCType::Free,
                     (int)m::BCType::VerticallyFixed, (int)m::BCType::Free,
                     (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P);
        y1 = y0;
    }
    pr.has_water = false;
    // 30 cycles at 60 steps/cycle: with xi = 5% the transient is ~e^{-0.05 * 2pi * 22} ~ 1e-3 of
    // the response by the final quarter of the record where the amplitude is read.
    m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Sweep";
    p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = kA; p.seismic_freq = freq;
    p.damping_ratio = cs.xi; p.rayleigh_f1 = cs.ray_f1; p.rayleigh_f2 = cs.ray_f2;
    p.duration = 30.0 / freq; p.time_steps = 1800;
    pr.phases.push_back(p);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    if (!M.ok) { if (msg) *msg = M.message; return -1.0; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2 || !res.back().ok) {
        if (msg && !res.empty()) *msg = res.back().message;
        return -1.0;
    }
    if (out) *out = res.back();
    return fe_transfer(res.back(), kA);
}

// ============================================================================================
// (a) ORACLE SELF-CHECK: the single undamped layer must collapse to T = 1/cos(wH/Vs).
// ============================================================================================
void test_oracle_closed_form() {
    std::printf("-- (a) oracle self-check: undamped single layer == 1/cos(wH/Vs) --\n");
    const Layer L{20.0, 2.0, 32000.0};                    // Vs = sqrt(G/rho) = 126.49 m/s
    const double Vs = std::sqrt(L.G / L.rho);
    double worst = 0.0;
    for (double f : {0.3, 0.9, 1.27, 2.2, 3.6}) {         // spans below/near/above f1 = 1.58 Hz
        const double w = 2 * kPi * f;
        const double Tor = layered_transfer({L}, 0.0, 0.0, w);
        const double Tcf = std::fabs(1.0 / std::cos(w * L.h / Vs));
        worst = std::fmax(worst, std::fabs(Tor - Tcf) / Tcf);
    }
    std::printf("   worst |oracle - 1/cos| / |1/cos| = %.2e\n", worst);
    check(worst < 1e-12, "oracle reduces to the classic closed form (signs/BCs pinned)");
}

// ============================================================================================
// (b) SINGLE LAYER, DAMPED: FE (GUI path) vs oracle across the transfer curve.
// ============================================================================================
void test_single_layer() {
    std::printf("-- (b) single damped layer: FE == layered viscoelastic solution --\n");
    ColumnSpec cs;
    cs.top_down = {{20.0, 2.0, 32000.0}};                 // Vs = 126.49, f1 = Vs/4H = 1.581 Hz
    const double f1 = std::sqrt(cs.top_down[0].G / cs.top_down[0].rho) / (4 * cs.top_down[0].h);
    cs.ray_f1 = f1; cs.ray_f2 = 3 * f1; cs.xi = 0.05;
    double alpha, beta; rayleigh_ab(cs.ray_f1, cs.ray_f2, cs.xi, alpha, beta);
    double worst = 0.0;
    katai::app::SolveResult last;
    for (double f : {0.5 * f1, f1, 2.0 * f1}) {
        std::string msg;
        const double Tfe = run_column(cs, f, &msg, &last);
        if (Tfe < 0) { check(false, ("FE run solved (" + msg + ")").c_str()); continue; }
        const double Tor = layered_transfer(cs.top_down, alpha, beta, 2 * kPi * f);
        const double err = (Tfe - Tor) / Tor;
        std::printf("   f = %.3f Hz:  |T|_FE = %.4f   |T|_oracle = %.4f   err = %+.2f%%\n",
                    f, Tfe, Tor, 100 * err);
        worst = std::fmax(worst, std::fabs(err));
    }
    check(worst < 0.05, "single-layer FE transfer matches the viscoelastic solution (< 5%)");
    // The auto-computed model fundamental (inverse-power on K,M) must be the column's Vs/4H.
    std::printf("   model f1 (auto) = %.4f Hz   (Vs/4H = %.4f, err %+.2f%%)\n",
                last.dyn_model_f1, f1, 100 * (last.dyn_model_f1 - f1) / f1);
    check(std::fabs(last.dyn_model_f1 - f1) < 0.02 * f1,
          "auto f1 == Vs/4H on the uniform column (< 2%)");
    check(last.message.find("Model fundamental frequency") != std::string::npos,
          "the phase message reports the model's own f1");
    check(last.message.find("NOTE: the Rayleigh band") == std::string::npos,
          "a properly bracketing Rayleigh band raises NO warning (no false alarm)");
}

// ============================================================================================
// (c) TWO-LAYER PROFILE (soft over stiff): the SHAKE-class benchmark proper. The impedance
//     contrast shifts the resonances and amplitudes away from any single-layer formula -- only
//     the layered solution predicts them (fundamental below the soft layer's own quarter-wave,
//     second mode set by the interference of the two layers).
// ============================================================================================
void test_two_layer() {
    std::printf("-- (c) two-layer profile vs the layered transfer-matrix solution --\n");
    ColumnSpec cs;
    // Soft over stiff -- a classic amplification profile with a real impedance contrast:
    //   top:    8 m,  rho 1.8 (gamma 17.66), Vs 120  -> G = 25920 kPa
    //   bottom: 12 m, rho 2.1 (gamma 20.60), Vs 300  -> G = 189000 kPa
    // Quarter-wave travel-time estimate f1 ~ 1 / (4 (8/120 + 12/300)) = 2.34 Hz (exact from oracle).
    cs.top_down = {{8.0, 1.8, 25920.0}, {12.0, 2.1, 189000.0}};
    const double f1_est = 1.0 / (4.0 * (8.0 / 120.0 + 12.0 / 300.0));
    cs.ray_f1 = f1_est; cs.ray_f2 = 3 * f1_est; cs.xi = 0.05;
    double alpha, beta; rayleigh_ab(cs.ray_f1, cs.ray_f2, cs.xi, alpha, beta);

    // Locate the first two resonances of the ORACLE (peak search on a fine grid).
    auto peak_in = [&](double flo, double fhi) {
        double fb = flo, Tb = 0.0;
        for (double f = flo; f <= fhi; f += 0.005) {
            const double T = layered_transfer(cs.top_down, alpha, beta, 2 * kPi * f);
            if (T > Tb) { Tb = T; fb = f; }
        }
        return std::pair<double, double>(fb, Tb);
    };
    const auto [f_pk1, T_pk1] = peak_in(0.5, 4.0);
    const auto [f_pk2, T_pk2] = peak_in(4.0, 10.0);
    std::printf("   oracle: f1 = %.3f Hz (|T| = %.3f)   f2 = %.3f Hz (|T| = %.3f)   "
                "[travel-time est. f1 ~ %.3f]\n", f_pk1, T_pk1, f_pk2, T_pk2, f1_est);
    // ORACLE PIN, hand-derivable closed form: for two elastic layers on a rigid base the
    // resonances satisfy   tan(k1 h1) tan(k2 h2) = 1/alpha_imp,  alpha_imp = rho1 Vs1/(rho2 Vs2)
    // (free surface -> cosine mode in the top layer, rigid base -> sine in the bottom, continuity
    // of u and tau at the interface). The travel-time rule f1 ~ 1/(4 sum h/Vs) is only an
    // APPROXIMATION and degrades with impedance contrast (here alpha = 0.343: est. 2.34 Hz vs the
    // true 3.0 Hz -- the soft layer leans toward its own rigid-base quarter-wave 120/(4*8) = 3.75).
    // So the pin is the exact equation, and the estimate only brackets from below.
    {
        const double Vs1 = std::sqrt(cs.top_down[0].G / cs.top_down[0].rho);
        const double Vs2 = std::sqrt(cs.top_down[1].G / cs.top_down[1].rho);
        const double a_imp = cs.top_down[0].rho * Vs1 / (cs.top_down[1].rho * Vs2);
        // Undamped oracle peak (pole) on a fine grid -> must satisfy the characteristic equation.
        double f0 = 0.5, T0 = 0.0;
        for (double f = 0.5; f <= 4.0; f += 0.001) {
            const double T = layered_transfer(cs.top_down, 0.0, 0.0, 2 * kPi * f);
            if (T > T0) { T0 = T; f0 = f; }
        }
        const double w0 = 2 * kPi * f0;
        const double lhs = std::tan(w0 * cs.top_down[0].h / Vs1) *
                           std::tan(w0 * cs.top_down[1].h / Vs2) * a_imp;
        std::printf("   oracle pin: undamped f1 = %.3f Hz, tan*tan*alpha = %.4f (must be 1)\n", f0, lhs);
        check(std::fabs(lhs - 1.0) < 0.02,
              "oracle fundamental satisfies the two-layer characteristic equation (exact pin)");
        check(f0 > f1_est && f0 < Vs1 / (4.0 * cs.top_down[0].h),
              "fundamental sits between the travel-time estimate and the soft-layer quarter-wave");
    }
    check(T_pk1 > 3.0, "the profile really amplifies at its fundamental (the benchmark has teeth)");

    // FE at four frequencies spanning the curve: half the fundamental, both resonances, and the
    // trough between them. Peaks are the SAFEST comparison points (d|T|/dw = 0 there), the flank
    // point is the most sensitive one -- band chosen per point, stated honestly.
    struct Pt { double f, tol; const char* what; };
    const Pt pts[] = {
        {0.5 * f_pk1, 0.04, "below the fundamental (flank)"},
        {f_pk1, 0.05, "first resonance (peak value)"},
        {0.5 * (f_pk1 + f_pk2), 0.06, "trough between the modes"},
        {f_pk2, 0.07, "second resonance (peak value)"},
    };
    katai::app::SolveResult last;
    for (const auto& pt : pts) {
        std::string msg;
        const double Tfe = run_column(cs, pt.f, &msg, &last);
        if (Tfe < 0) { check(false, ("FE run solved (" + msg + ")").c_str()); continue; }
        const double Tor = layered_transfer(cs.top_down, alpha, beta, 2 * kPi * pt.f);
        const double err = (Tfe - Tor) / Tor;
        std::printf("   f = %.3f Hz (%s):\n      |T|_FE = %.4f   |T|_oracle = %.4f   err = %+.2f%%"
                    "  (band %.0f%%)\n", pt.f, pt.what, Tfe, Tor, 100 * err, 100 * pt.tol);
        char label[160];
        std::snprintf(label, sizeof(label), "two-layer FE == layered solution at %s", pt.what);
        check(std::fabs(err) < pt.tol, label);
    }

    // The auto f1 must agree with the ORACLE's resonance -- two fully independent routes to the
    // same number (inverse-power on the assembled 2D K,M vs the layered transfer-matrix peak).
    std::printf("   model f1 (auto) = %.4f Hz   (oracle peak = %.4f, err %+.2f%%)\n",
                last.dyn_model_f1, f_pk1, 100 * (last.dyn_model_f1 - f_pk1) / f_pk1);
    check(std::fabs(last.dyn_model_f1 - f_pk1) < 0.02 * f_pk1,
          "auto f1 == layered-oracle resonance on the two-layer profile (< 2%)");
    // The 22% travel-time trap, productized: the informational line always appears; the WARNING
    // appears only when the Rayleigh band genuinely mis-damps the fundamental. Here the (deliberate)
    // travel-time-based target still BRACKETS f1 with f2 = 3 f_R1, so no warning may fire...
    check(last.message.find("NOTE: the Rayleigh band") == std::string::npos,
          "bracketing band (travel-time f_R1, 3x span): informational only, no warning");
    // ...but a band placed entirely below the true fundamental MUST warn.
    ColumnSpec bad = cs;
    bad.ray_f1 = 0.5; bad.ray_f2 = 1.0;   // true f1 = 3.0 Hz sits far outside
    std::string msg2;
    katai::app::SolveResult res_bad;
    const double Tbad = run_column(bad, f_pk1, &msg2, &res_bad);
    check(Tbad >= 0, "mistuned-band run solved");
    if (Tbad >= 0) {
        std::printf("   mistuned band [0.5, 1.0] Hz -> message warns: %s\n",
                    res_bad.message.find("NOTE: the Rayleigh band") != std::string::npos ? "yes" : "NO");
        check(res_bad.message.find("NOTE: the Rayleigh band") != std::string::npos,
              "a band that misses f1 raises the Rayleigh warning (the guardrail has teeth)");
    }
}

}  // namespace

int main() {
    std::printf("BENCHMARK Wave-1: multi-layer 1D site response vs the SHAKE-class\n"
                "layered viscoelastic transfer-matrix solution (independent oracle)\n\n");
    test_oracle_closed_form();
    std::printf("\n");
    test_single_layer();
    std::printf("\n");
    test_two_layer();
    if (g_failures == 0) {
        std::printf("\nOK: the GUI dynamic path reproduces the layered site-response solution "
                    "(single- and two-layer, peaks and flanks)\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
