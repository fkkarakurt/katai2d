// Transient (time-dependent) saturated groundwater flow — W1 verification (head form).
// PDE (1D): S_s ∂h/∂t = k ∂²h/∂y²,  diffusion coefficient D = k/S_s.  Three closed-form oracles:
//   (A) Terzaghi-style isochrone dissipation: uniform initial excess head, top-drained column
//       u(Z,Tv)/u0 = Σ (2/M) sin(M·Z) exp(-M²Tv), M=(2m+1)π/2, Z=(Hc-y)/Hc (depth from drainage),
//       Tv = D t / Hc² (single drainage). (Carslaw-Jaeger/Terzaghi Fourier series.)
//   (B) Sudden head step (Carslaw-Jaeger/Bear): semi-infinite column, h=0 initially, h=hs at the
//       end for t>0; h(d,t)/hs = erfc( d / (2√(Dt)) ), d = distance to the end. (Time-constant
//       Dirichlet jump.)
//   (C) Steady-state regression: constant head gradient → long time → linear profile (= steady seepage).
// See docs/references/transient-unsaturated-flow-formulation.md (W1); PLAXIS 2D Sci.Man §3.1-3.3.
#include <katai/analysis/transient_flow.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::HeadBoundary;
using katai::core::Permeability;
using katai::core::solve_transient_flow;
using katai::core::TransientFlowResult;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// Terzaghi point value (isochrone): excess-head ratio u/u0 at normalized depth Z.
double terzaghi_uz(double Z, double Tv) {
    double s = 0.0;
    for (int m = 0; m < 200; ++m) {
        const double M = (2 * m + 1) * kPi / 2.0;
        s += (2.0 / M) * std::sin(M * Z) * std::exp(-M * M * Tv);
    }
    return s;
}

// For a column (nx=1), returns the left-edge nodes sorted by y (a vertical line).
std::vector<int> sorted_left_line(const Mesh& mesh) {
    std::vector<int> line = mesh.left_nodes;
    std::sort(line.begin(), line.end(), [&](int a, int b) { return mesh.y[a] < mesh.y[b]; });
    return line;
}

// (A) Terzaghi-style: top-drained column, uniform initial excess head dissipates.
void test_terzaghi_dissipation() {
    constexpr double W = 0.1, Hc = 1.0;
    constexpr double k = 1.0e-4, Ss = 1.0e-2;   // D = k/Ss = 1e-2
    const double D = k / Ss;
    constexpr double u0 = 10.0;                  // uniform initial excess head

    RectangularDomain domain{0.0, 0.0, W, Hc, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 1, 40);  // slender column
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> storage = {Ss};

    // Top nodes drained (h=0); the rest of the boundary impermeable (natural). Initial h=u0.
    std::vector<char> presc(mesh.node_count, 0);
    for (int n : mesh.top_nodes) presc[n] = 1;
    HeadBoundary bc = [](int, double) { return 0.0; };
    std::vector<double> h0(mesh.node_count, u0);
    for (int n : mesh.top_nodes) h0[n] = 0.0;  // t=0 boundary (IC↔BC consistency)

    const double dt = 2.0;        // Tv step = D·dt/Hc² = 0.02
    const int nsteps = 100;       // Tv 0..2.0 (at Tv≈2, U≈0.99 → excess head negligible)
    const TransientFlowResult r = solve_transient_flow(mesh, perm, storage, presc, bc, h0, dt, nsteps);

    const std::vector<int> line = sorted_left_line(mesh);
    std::printf("  (A) Terzaghi dissipation: D=%.4f  Hc=%.1f  u0=%.1f\n", D, Hc, u0);
    int checked = 0;
    // Interior depths (away from the drainage boundary layer) × a few Tv.
    for (double Tv : {0.2, 0.4, 0.7}) {
        const double t = Tv * Hc * Hc / D;
        // The nearest recorded time step.
        int si = 0; double best = 1e30;
        for (size_t i = 0; i < r.times.size(); ++i)
            if (std::fabs(r.times[i] - t) < best) { best = std::fabs(r.times[i] - t); si = (int)i; }
        const double Tv_act = D * r.times[si] / (Hc * Hc);
        for (int n : line) {
            const double Z = (Hc - mesh.y[n]) / Hc;        // depth from the (top) drainage
            if (Z < 0.3 || Z > 0.7) continue;              // skip the boundary layer/bottom
            const double ufe = r.head[si][n] / u0;
            const double uth = terzaghi_uz(Z, Tv_act);
            const double err = std::fabs(ufe - uth);
            if (Z > 0.45 && Z < 0.55)
                std::printf("    Tv=%.3f Z=%.2f  u/u0_FE=%.4f  u/u0_th=%.4f  |Δ|=%.4f\n",
                            Tv_act, Z, ufe, uth, err);
            check(err < 0.03, "transient u(Z,Tv) matches Terzaghi isochrone within 0.03");
            ++checked;
        }
    }
    check(checked >= 6, "(A) checked several depth/time points");
    // Steady state: t→large, excess head ~0.
    double maxh = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) maxh = std::fmax(maxh, std::fabs(r.head.back()[n]));
    std::printf("    final max|h| = %.4e (→0)\n", maxh);
    check(maxh < 0.02 * u0, "(A) excess head fully dissipated at large Tv");
}

// (B) Sudden head step: semi-infinite column, h=0 initially, h=hs at the top end for t>0 → erfc.
void test_erfc_step() {
    constexpr double W = 0.1, L = 2.0;
    constexpr double k = 1.0e-4, Ss = 1.0e-2;   // D = 1e-2
    const double D = k / Ss;
    constexpr double hs = 5.0;

    RectangularDomain domain{0.0, 0.0, W, L, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 1, 80);  // resolution (sharp front)
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> storage = {Ss};

    std::vector<char> presc(mesh.node_count, 0);
    for (int n : mesh.top_nodes) presc[n] = 1;          // top end h=hs (t>0)
    HeadBoundary bc = [&](int, double) { return hs; };
    std::vector<double> h0(mesh.node_count, 0.0);       // initial h=0; bottom impermeable (natural)

    const double dt = 0.25;
    const int nsteps = 36;            // t_end = 9.0 → 2√(Dt)=0.6 (front in the upper half of the column)
    const TransientFlowResult r = solve_transient_flow(mesh, perm, storage, presc, bc, h0, dt, nsteps);

    const double t = r.times.back();
    const double front = 2.0 * std::sqrt(D * t);
    std::printf("  (B) erfc step: D=%.4f  hs=%.1f  t=%.2f  2sqrt(Dt)=%.3f\n", D, hs, t, front);
    const std::vector<int> line = sorted_left_line(mesh);
    int checked = 0;
    for (int n : line) {
        const double d = L - mesh.y[n];                // distance to the top end
        if (d < 0.15 || d > 0.9) continue;             // skip too-close (mesh) / beyond-front (≈0)
        const double hfe = r.head[t == 0 ? 0 : (int)r.times.size() - 1][n];
        const double hth = hs * std::erfc(d / (2.0 * std::sqrt(D * t)));
        const double rel = std::fabs(hfe - hth) / hs;
        if (d > 0.35 && d < 0.55)
            std::printf("    d=%.3f  h_FE=%.4f  h_erfc=%.4f  relΔ=%.3f%%\n", d, hfe, hth, 100 * rel);
        check(rel < 0.03, "transient h(d,t) matches Carslaw-Jaeger erfc within 3% of hs");
        ++checked;
    }
    check(checked >= 5, "(B) checked several depths along the erfc front");
}

// (C) Steady-state regression: constant head gradient → linear profile (= the steady seepage solution).
void test_steady_regression() {
    constexpr double W = 0.2, Hc = 1.0;
    constexpr double k = 1.0e-3, Ss = 1.0e-2;
    constexpr double Htop = 3.0, Hbot = 1.0;     // constant Dirichlet ends

    RectangularDomain domain{0.0, 0.0, W, Hc, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 1, 20);
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> storage = {Ss};

    std::vector<char> presc(mesh.node_count, 0);
    std::vector<double> bcval(mesh.node_count, 0.0);
    for (int n : mesh.top_nodes)    { presc[n] = 1; bcval[n] = Htop; }
    for (int n : mesh.bottom_nodes) { presc[n] = 1; bcval[n] = Hbot; }
    HeadBoundary bc = [&](int n, double) { return bcval[n]; };
    std::vector<double> h0(mesh.node_count, 0.0);  // a bad initial guess

    const double dt = 50.0;          // large dt → fast to steady state
    const int nsteps = 60;
    const TransientFlowResult r = solve_transient_flow(mesh, perm, storage, presc, bc, h0, dt, nsteps);

    std::printf("  (C) steady regression: linear head profile Htop=%.1f Hbot=%.1f\n", Htop, Hbot);
    double maxerr = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double exact = Hbot + (Htop - Hbot) * (mesh.y[n] / Hc);
        maxerr = std::fmax(maxerr, std::fabs(r.head.back()[n] - exact));
    }
    std::printf("    max|h_FE - h_linear| at large t = %.3e\n", maxerr);
    check(maxerr < 1e-3, "(C) transient steady-state recovers the linear (seepage) head profile");
}

}  // namespace

int main() {
    test_terzaghi_dissipation();
    test_erfc_step();
    test_steady_regression();
    if (g_failures == 0) {
        std::printf("OK: transient flow verified (Terzaghi isochrone + Carslaw-Jaeger erfc + steady regression)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
