// COMPLIANT (absorbing) BASE -- V&V against the layered solution with a RADIATION base condition
// (formulation locked in docs/references/dynamic-seismic-formulation.md sec 11; PLAXIS Sci 6.3.2
// Joyner & Chen factor 2 + Tut 17.8.5 half-of-within input convention).
//
// THE ORACLE (total motion, shares no code with the solver): layers u_j = a_j cos(k_j z) +
// b_j sin(k_j z) (no body force -- the input enters through the base), with the FE's Rayleigh
// damping mapped exactly into the layers (G* = G(1+iw beta), rho* = rho(1-i alpha/w)) and the
// halfspace as a REAL impedance rho_r Vs_r (a constant dashpot -- exactly what the FE assembles).
// Base condition (hand-derived from tau balance at the bottom face, z up, U_up = 1):
//     -i w rho_r Vs_r a_0 + G*_0 k_0 b_0 = -2 i w rho_r Vs_r
// Interfaces: continuity of u and G* u'. Surface: u' = 0. |T_outcrop| = |u(surface)| / 2.
// Closed-form pin: the undamped single layer collapses to the classic
//     |T| = 1 / sqrt(cos^2(kH) + alpha^2 sin^2(kH)),   alpha = rho1 Vs1 / (rho_r Vs_r),
// with the FINITE resonance 1/alpha -- radiation damping. That finiteness is the entire point of
// the feature: a rigid base reflects all downgoing energy and (with light material damping) rings
// far above 1/alpha; the compliant base must SIT on it.
//
// FE measurement: the GUI path with Phase.seismic_compliant_base = true. a_g is the within motion,
// so the outcrop-equivalent amplitude is A and |T|_FE = steady-state |a_surface,total| / A.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
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

struct Layer { double h, rho, G; };   // top-down

// |T_outcrop(w)| for layers over an elastic halfspace (impedance rho_r Vs_r), Rayleigh-damped soil.
double radiating_transfer(const std::vector<Layer>& top_down, double rhoVs_half,
                          double alpha, double beta, double w) {
    using C = std::complex<double>;
    const C i01(0.0, 1.0);
    const int n = (int)top_down.size();
    std::vector<Layer> up(top_down.rbegin(), top_down.rend());   // bottom-up
    std::vector<C> k(n), Gs(n);
    for (int j = 0; j < n; ++j) {
        const C Gstar = up[j].G * (1.0 + i01 * (w * beta));
        const C rstar = up[j].rho * (1.0 - i01 * (alpha / w));
        Gs[j] = Gstar;
        k[j] = w * std::sqrt(rstar / Gstar);
    }
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
    Eigen::VectorXcd r = Eigen::VectorXcd::Zero(2 * n);
    int eq = 0;
    // Base (radiation, U_up = 1): -i w Z a_0 + G*_0 k_0 b_0 = -2 i w Z,  Z = rho_r Vs_r.
    S(eq, 0) = -i01 * w * rhoVs_half;
    S(eq, 1) = Gs[0] * k[0];
    r(eq) = -2.0 * i01 * w * rhoVs_half; ++eq;
    for (int j = 0; j + 1 < n; ++j) {
        const C ch = std::cos(k[j] * up[j].h), sh = std::sin(k[j] * up[j].h);
        S(eq, 2 * j) = ch; S(eq, 2 * j + 1) = sh; S(eq, 2 * (j + 1)) = -1.0; r(eq) = 0.0; ++eq;
        S(eq, 2 * j) = -Gs[j] * k[j] * sh;
        S(eq, 2 * j + 1) = Gs[j] * k[j] * ch;
        S(eq, 2 * (j + 1) + 1) = -Gs[j + 1] * k[j + 1];
        r(eq) = 0.0; ++eq;
    }
    {
        const int t = n - 1;
        const C ch = std::cos(k[t] * up[t].h), sh = std::sin(k[t] * up[t].h);
        S(eq, 2 * t) = -sh; S(eq, 2 * t + 1) = ch; r(eq) = 0.0; ++eq;
    }
    const Eigen::VectorXcd x = S.partialPivLu().solve(r);
    const int t = n - 1;
    const C u_surf = x(2 * t) * std::cos(k[t] * up[t].h) + x(2 * t + 1) * std::sin(k[t] * up[t].h);
    return std::abs(u_surf) / 2.0;
}

void rayleigh_ab(double f1, double f2, double xi, double& alpha, double& beta) {
    const double w1 = 2 * kPi * f1, w2 = 2 * kPi * f2;
    alpha = 2.0 * xi * w1 * w2 / (w1 + w2);
    beta = 2.0 * xi / (w1 + w2);
}

double fe_transfer(const katai::app::SolveResult& R, double A) {
    if (R.dyn_time.empty()) return -1.0;
    const double t_end = R.dyn_time.back(), t_from = 0.75 * t_end;
    double amp = 0.0;
    for (size_t i = 0; i < R.dyn_time.size(); ++i)
        if (R.dyn_time[i] >= t_from) amp = std::fmax(amp, std::fabs(R.dyn_surface_ax[i]));
    return amp / A;
}

struct ColumnSpec {
    std::vector<Layer> top_down;
    double nu = 0.3;
    double ray_f1 = 1.0, ray_f2 = 3.0, xi = 0.05;
    bool compliant = true, free_field = false, nonlinear = false;
    bool free_sides = false;   // Free instead of VerticallyFixed (the free-field-side setup)
    double mc_c = -1.0;        // > 0: Mohr-Coulomb (phi = 0) with this cohesion instead of LE
    double width = 2.0;
};
double run_column(const ColumnSpec& cs, double freq, std::string* msg = nullptr,
                  katai::app::SolveResult* out = nullptr) {
    const double kW = cs.width;
    constexpr double kGrav = 9.81, kA = 1.0;
    m::Project pr;
    double ytop = 0.0;
    for (const auto& L : cs.top_down) ytop += L.h;
    for (const auto& L : cs.top_down) {
        m::Material s;
        s.model = cs.mc_c > 0.0 ? m::SoilModel::MohrCoulomb : m::SoilModel::LinearElastic;
        s.E = 2.0 * (1.0 + cs.nu) * L.G; s.nu = cs.nu;
        s.c = cs.mc_c > 0.0 ? cs.mc_c : 1.0; s.phi = 0.0; s.psi = 0.0;
        s.gamma_unsat = L.rho * kGrav; s.gamma_sat = s.gamma_unsat; s.e_init = 0.5;
        pr.materials.push_back(s);
    }
    const int side_bc = cs.free_sides ? (int)m::BCType::Free : (int)m::BCType::VerticallyFixed;
    double y1 = ytop;
    for (size_t li = 0; li < cs.top_down.size(); ++li) {
        const double y0 = y1 - cs.top_down[li].h;
        m::SoilPolygon P; P.material = (int)li;
        P.x = {0, kW, kW, 0}; P.y = {y0, y0, y1, y1};
        const bool base = li + 1 == cs.top_down.size();
        P.edge_bc = {base ? (int)m::BCType::FullyFixed : (int)m::BCType::Free,
                     side_bc, (int)m::BCType::Free, side_bc};
        pr.polygons.push_back(P);
        y1 = y0;
    }
    pr.has_water = false;
    m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "CB";
    p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = kA; p.seismic_freq = freq;
    p.damping_ratio = cs.xi; p.rayleigh_f1 = cs.ray_f1; p.rayleigh_f2 = cs.ray_f2;
    p.duration = 20.0 / freq; p.time_steps = 1200;
    p.seismic_compliant_base = cs.compliant;
    p.seismic_free_field = cs.free_field;
    p.dynamic_nonlinear = cs.nonlinear;
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

// (a) Oracle pin: undamped single layer over a halfspace == the classic closed form.
void test_oracle_pin() {
    std::printf("-- (a) oracle pin: undamped layer over halfspace == 1/sqrt(cos^2 + a^2 sin^2) --\n");
    const Layer L{20.0, 2.0, 32000.0};                 // Vs = 126.49
    const double Vs = std::sqrt(L.G / L.rho), a_imp = 0.35;
    const double Z = L.rho * Vs / a_imp;               // halfspace impedance for that alpha
    double worst = 0.0;
    for (double f : {0.4, 1.0, 1.581, 2.6, 3.9}) {
        const double w = 2 * kPi * f, kH = w * L.h / Vs;
        const double Tor = radiating_transfer({L}, Z, 0.0, 0.0, w);
        const double Tcf = 1.0 / std::sqrt(std::pow(std::cos(kH), 2) +
                                           a_imp * a_imp * std::pow(std::sin(kH), 2));
        worst = std::fmax(worst, std::fabs(Tor - Tcf) / Tcf);
    }
    std::printf("   worst |oracle - closed form| = %.2e   (resonance cap 1/alpha = %.3f)\n",
                worst, 1.0 / a_imp);
    check(worst < 1e-12, "radiation oracle == the classic layer-over-halfspace closed form");
}

// (b) UNIFORM column + compliant base: the halfspace CONTINUES the column (alpha = 1) -> no
//     impedance contrast, no trapped energy, |T| ~ 1 at EVERY frequency -- including the rigid-
//     base resonance, where a missing factor 2 (|T| -> 0.5) or a missing dashpot (|T| >> 1) would
//     be unmistakable. The sharpest single check of the whole boundary condition.
void test_uniform_column() {
    std::printf("-- (b) uniform column: perfect radiation, |T| ~ 1 at every frequency --\n");
    ColumnSpec cs;
    cs.top_down = {{20.0, 2.0, 32000.0}};
    const double f1 = std::sqrt(cs.top_down[0].G / cs.top_down[0].rho) / (4 * cs.top_down[0].h);
    cs.ray_f1 = f1; cs.ray_f2 = 3 * f1; cs.xi = 0.05;
    double alpha, beta; rayleigh_ab(cs.ray_f1, cs.ray_f2, cs.xi, alpha, beta);
    const double Z = cs.top_down[0].rho * std::sqrt(cs.top_down[0].G / cs.top_down[0].rho);
    for (double f : {0.7 * f1, f1, 1.8 * f1}) {
        std::string msg;
        katai::app::SolveResult R;
        const double Tfe = run_column(cs, f, &msg, &R);
        if (Tfe < 0) { check(false, ("FE run solved (" + msg + ")").c_str()); continue; }
        const double Tor = radiating_transfer(cs.top_down, Z, alpha, beta, 2 * kPi * f);
        const double err = (Tfe - Tor) / Tor;
        std::printf("   f = %.3f Hz:  |T|_FE = %.4f   |T|_oracle = %.4f   err = %+.2f%%\n",
                    f, Tfe, Tor, 100 * err);
        check(std::fabs(err) < 0.05, "uniform compliant column == radiation oracle (< 5%)");
        check(Tor > 0.75 && Tor < 1.05, "oracle itself stays ~1 (perfect radiation sanity)");
        if (f == f1)
            check(R.message.find("COMPLIANT BASE") != std::string::npos,
                  "the phase message declares the compliant base and its input convention");
    }
}

// (c) TWO-LAYER profile (soft over stiff; the halfspace continues the stiff layer) + the TEETH:
//     the same profile on a RIGID base rings at |T| ~ 17 (site-response benchmark); the compliant
//     base must sit on the radiation oracle's FINITE peak, several times lower.
void test_two_layer_and_teeth() {
    std::printf("-- (c) two-layer + teeth: rigid base rings, compliant base radiates --\n");
    ColumnSpec cs;
    cs.top_down = {{8.0, 1.8, 25920.0}, {12.0, 2.1, 189000.0}};   // Vs 120 / 300
    const double f1_est = 1.0 / (4.0 * (8.0 / 120.0 + 12.0 / 300.0));
    cs.ray_f1 = f1_est; cs.ray_f2 = 3 * f1_est; cs.xi = 0.05;
    double alpha, beta; rayleigh_ab(cs.ray_f1, cs.ray_f2, cs.xi, alpha, beta);
    const double Z = cs.top_down[1].rho * std::sqrt(cs.top_down[1].G / cs.top_down[1].rho);

    // Oracle peak of the RADIATING profile (finite; the damped scan).
    double f_pk = 0.5, T_pk = 0.0;
    for (double f = 0.5; f <= 6.0; f += 0.005) {
        const double T = radiating_transfer(cs.top_down, Z, alpha, beta, 2 * kPi * f);
        if (T > T_pk) { T_pk = T; f_pk = f; }
    }
    std::printf("   radiating oracle peak: f = %.3f Hz, |T| = %.3f (rigid-base peak was ~17)\n",
                f_pk, T_pk);
    check(T_pk > 1.5 && T_pk < 8.0, "radiation caps the resonance at a finite, moderate value");

    struct Pt { double f; const char* what; };
    for (const Pt& pt : {Pt{0.5 * f_pk, "flank"}, Pt{f_pk, "peak"}, Pt{1.6 * f_pk, "above"}}) {
        std::string msg;
        const double Tfe = run_column(cs, pt.f, &msg);
        if (Tfe < 0) { check(false, ("FE run solved (" + msg + ")").c_str()); continue; }
        const double Tor = radiating_transfer(cs.top_down, Z, alpha, beta, 2 * kPi * pt.f);
        const double err = (Tfe - Tor) / Tor;
        std::printf("   f = %.3f Hz (%s):  |T|_FE = %.4f   |T|_oracle = %.4f   err = %+.2f%%\n",
                    pt.f, pt.what, Tfe, Tor, 100 * err);
        char label[128];
        std::snprintf(label, sizeof(label), "two-layer compliant == radiation oracle at the %s", pt.what);
        check(std::fabs(err) < 0.06, label);
    }

    // TEETH: rigid vs compliant at the compliant peak's neighbourhood -- the rigid base must ring
    // far above the radiating cap (that overshoot IS the box-resonance error the feature removes).
    ColumnSpec rigid = cs; rigid.compliant = false;
    std::string msg;
    const double T_rigid = run_column(rigid, 3.005, &msg);   // rigid-base fundamental (benchmark)
    const double T_comp = radiating_transfer(cs.top_down, Z, alpha, beta, 2 * kPi * 3.005);
    std::printf("   at 3.005 Hz: rigid |T|_FE = %.3f  vs  radiating cap = %.3f  (ratio %.1fx)\n",
                T_rigid, T_comp, T_rigid / T_comp);
    check(T_rigid > 3.0 * T_comp,
          "the rigid base rings >= 3x above the radiating solution (the error this feature removes)");
}

// (d) The FORMER scope refusals are now capabilities ("an honest refusal is not a capability"):
//     (e) compliant + NONLINEAR: on an elastic column the per-step-Newton path must reproduce the
//         linear compliant solution (the linear-limit identity -- any total-motion baseline error
//         in the nonlinear wiring breaks it); a low-cohesion Tresca column must then DEPART
//         (the nonlinearity genuinely engages under absorbing-base shaking).
//     (g) compliant + FREE-FIELD sides: on a laterally-uniform column the free-field sides must
//         not disturb the 1D radiating physics -- FE == the radiation oracle, same as SH sides
//         (the 1D side columns now solve on a compliant base too; a rigid-base side column would
//         re-impose the resonance the base just absorbed).
void test_former_refusals_are_capabilities() {
    std::printf("-- (d) former refusals are now capabilities (nonlinear + free-field sides) --\n");
    ColumnSpec cs;
    cs.top_down = {{20.0, 2.0, 32000.0}};
    const double f1 = std::sqrt(cs.top_down[0].G / cs.top_down[0].rho) / (4 * cs.top_down[0].h);
    cs.ray_f1 = f1; cs.ray_f2 = 3 * f1; cs.xi = 0.05;
    double alpha, beta; rayleigh_ab(cs.ray_f1, cs.ray_f2, cs.xi, alpha, beta);
    const double Z = cs.top_down[0].rho * std::sqrt(cs.top_down[0].G / cs.top_down[0].rho);

    // (e) linear limit: nonlinear compliant == linear compliant on an elastic column.
    std::string msg;
    const double T_lin = run_column(cs, f1, &msg);
    ColumnSpec nl = cs; nl.nonlinear = true;
    const double T_nl = run_column(nl, f1, &msg);
    if (T_lin < 0 || T_nl < 0) { check(false, ("compliant runs solved (" + msg + ")").c_str()); return; }
    std::printf("   nonlinear linear-limit: |T|_lin = %.6f  |T|_nl = %.6f  (rel diff %.2e)\n",
                T_lin, T_nl, std::fabs(T_nl - T_lin) / T_lin);
    check(std::fabs(T_nl - T_lin) < 1e-6 * T_lin,
          "compliant + NONLINEAR reduces to the linear compliant solution on elastic soil");
    // ...and genuinely yields when the strength is low (Tresca c; the nonlinearity has teeth).
    ColumnSpec mc = nl; mc.mc_c = 8.0;
    const double T_mc = run_column(mc, f1, &msg);
    if (T_mc < 0) { check(false, ("MC compliant run solved (" + msg + ")").c_str()); return; }
    std::printf("   MC (c=8, phi=0):        |T|_mc = %.4f  (departs %.1f%% from elastic)\n",
                T_mc, 100 * std::fabs(T_mc - T_lin) / T_lin);
    check(std::fabs(T_mc - T_lin) > 0.01 * T_lin,
          "a weak Tresca column departs from elastic under absorbing-base shaking (yield engages)");

    // (g) free-field sides on the compliant column: the 1D radiating identity must survive.
    ColumnSpec ff = cs; ff.free_field = true; ff.free_sides = true; ff.width = 10.0;
    const double T_ff = run_column(ff, f1, &msg);
    if (T_ff < 0) { check(false, ("compliant+FF run solved (" + msg + ")").c_str()); return; }
    const double Tor = radiating_transfer(cs.top_down, Z, alpha, beta, 2 * kPi * f1);
    std::printf("   FF sides + compliant:   |T|_FE = %.4f  |T|_oracle = %.4f  (err %+.2f%%)\n",
                T_ff, Tor, 100 * (T_ff - Tor) / Tor);
    check(std::fabs(T_ff - Tor) < 0.05 * Tor,
          "compliant base + free-field sides reproduces the 1D radiating solution (< 5%)");
}

}  // namespace

int main() {
    std::printf("COMPLIANT (absorbing) BASE vs the layered radiation solution (independent oracle)\n\n");
    test_oracle_pin();
    std::printf("\n");
    test_uniform_column();
    std::printf("\n");
    test_two_layer_and_teeth();
    std::printf("\n");
    test_former_refusals_are_capabilities();
    if (g_failures == 0) {
        std::printf("\nOK: the compliant base radiates -- uniform column |T| ~ 1, layered profile on "
                    "the finite radiating solution, rigid base rings above it, and the nonlinear + "
                    "free-field combinations hold the same identities\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
