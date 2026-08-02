// Transient unsaturated (Richards) flow — W2b verification (head form, Celia 1990
// modified Picard).
//   (1) Saturated limit: a column that stays fully saturated → the excess head (h−y)
//       follows the Terzaghi isochrone (closed form; exercises the unsaturated solver's
//       saturated branch).
//   (2) Global mass conservation (the Celia criterion): infiltration into a dry column →
//       the total water-volume gain = the time-integrated boundary Darcy flux (the
//       signature property of modified Picard). |1 − ratio| small.
//   (3) Philip √t sorptivity: early-time cumulative infiltration I(t) ∝ √t (unsaturated
//       dynamics).
// Kaynaklar: Celia, Bouloutas & Zarba (1990) WRR 26:1483; Philip (1957) sorptivite; van Genuchten/Mualem.
#include <katai/analysis/transient_flow.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::HeadBoundary;
using katai::core::Permeability;
using katai::core::solve_transient_unsaturated_flow;
using katai::core::UnsaturatedFlowResult;
using katai::core::WaterRetention;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
double terzaghi_uz(double Z, double Tv) {
    double s = 0.0;
    for (int m = 0; m < 200; ++m) {
        const double M = (2 * m + 1) * kPi / 2.0;
        s += (2.0 / M) * std::sin(M * Z) * std::exp(-M * M * Tv);
    }
    return s;
}
std::vector<int> sorted_left_line(const Mesh& mesh) {
    std::vector<int> line = mesh.left_nodes;
    std::sort(line.begin(), line.end(), [&](int a, int b) { return mesh.y[a] < mesh.y[b]; });
    return line;
}

// (1) Saturated limit: the column stays saturated (ψ=h−y≥0) → excess head follows the Terzaghi isochrone (D=k/S_s).
void test_saturated_terzaghi() {
    constexpr double W = 0.1, Hc = 1.0;
    constexpr double k = 1.0e-4, ss = 1.0e-2;        // D = k/ss = 1e-2
    const double D = k / ss;
    constexpr double e0 = 10.0;                       // uniform initial excess head (h−y)
    RectangularDomain domain{0.0, 0.0, W, Hc, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 1, 40);
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> n = {0.4}, sss = {ss};
    const std::vector<WaterRetention> ret = {WaterRetention{}};

    // Total-head form: the dissipating field is e = h − y_top. Initial h ≡ y_top+e0
    // (uniform → ψ=h−y≥e0>0 saturated); top drainage h = y_top (e=0). Steady h ≡ y_top
    // (e→0) = Terzaghi. H·h=H·e (y_top constant → null space).
    const double y_top = Hc;
    std::vector<char> presc(mesh.node_count, 0);
    for (int nd : mesh.top_nodes) presc[nd] = 1;
    std::vector<double> h0(mesh.node_count, y_top + e0);
    for (int nd : mesh.top_nodes) h0[nd] = y_top;                            // top drainage e=0
    HeadBoundary bc = [&](int, double) { return y_top; };                    // h=y_top → e=0

    const double dt = 2.0; const int nsteps = 100;   // Tv 0..2
    const UnsaturatedFlowResult r = solve_transient_unsaturated_flow(mesh, perm, n, sss, ret, presc, bc, h0, dt, nsteps);
    check(r.converged, "(1) unsaturated solver converged (saturated branch)");
    std::printf("  (1) saturated-limit Terzaghi: D=%.4f  max Picard iters=%d\n", D, r.max_picard_iters);
    const std::vector<int> line = sorted_left_line(mesh);
    int checked = 0;
    for (double Tv : {0.2, 0.4, 0.7}) {
        const double t = Tv * Hc * Hc / D;
        int si = 0; double best = 1e30;
        for (size_t i = 0; i < r.times.size(); ++i) if (std::fabs(r.times[i] - t) < best) { best = std::fabs(r.times[i] - t); si = (int)i; }
        const double Tva = D * r.times[si] / (Hc * Hc);
        for (int nd : line) {
            const double Z = (Hc - mesh.y[nd]) / Hc;
            if (Z < 0.3 || Z > 0.7) continue;
            const double efe = (r.head[si][nd] - y_top) / e0;
            const double eth = terzaghi_uz(Z, Tva);
            if (Z > 0.45 && Z < 0.55)
                std::printf("    Tv=%.3f Z=%.2f  e/e0_FE=%.4f  e/e0_th=%.4f  |Δ|=%.4f\n", Tva, Z, efe, eth, std::fabs(efe - eth));
            check(std::fabs(efe - eth) < 0.03, "saturated-branch excess head matches Terzaghi");
            ++checked;
        }
    }
    check(checked >= 6, "(1) checked several depth/time points");
}

// (2)+(3) Ponded infiltration into a dry column: mass conservation + Philip √t.
void test_infiltration() {
    constexpr double W = 0.1, L = 1.0;
    constexpr double k = 1.0e-3, ss = 0.0;           // pure Richards (specific storage neglected — standard)
    constexpr double n = 0.40, psi_dry = 1.0;        // initial suction 1 m (dry)
    RectangularDomain domain{0.0, 0.0, W, L, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 1, 60);
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> nv = {n}, sss = {ss};
    const std::vector<WaterRetention> ret = {WaterRetention{2.0, 2.0, 0.5, 0.0, 1.0}};  // g_a=2,g_n=2

    std::vector<char> presc(mesh.node_count, 0);
    std::vector<double> bcval(mesh.node_count, 0.0);
    for (int nd : mesh.top_nodes)    { presc[nd] = 1; bcval[nd] = mesh.y[nd]; }            // ponded h_p=0 (S=1)
    for (int nd : mesh.bottom_nodes) { presc[nd] = 1; bcval[nd] = mesh.y[nd] - psi_dry; }  // kuru tut
    HeadBoundary bc = [&](int nd, double) { return bcval[nd]; };
    std::vector<double> h0(mesh.node_count);
    for (int nd = 0; nd < mesh.node_count; ++nd) h0[nd] = mesh.y[nd] - psi_dry;            // dry start
    for (int nd : mesh.top_nodes) h0[nd] = mesh.y[nd];

    const double dt = 0.5; const int nsteps = 40;
    const UnsaturatedFlowResult r = solve_transient_unsaturated_flow(mesh, perm, nv, sss, ret, presc, bc, h0, dt, nsteps);
    check(r.converged, "(2) infiltration solver converged");
    std::printf("  (2) infiltration: max Picard iters=%d  steps=%zu\n", r.max_picard_iters, r.times.size());

    // Mass conservation: Σ_step dt·darcy_influx  ==  the total water-volume gain.
    double cum_influx = 0.0;
    for (size_t i = 1; i < r.times.size(); ++i) cum_influx += dt * r.darcy_influx[i];
    const double dvol = r.water_volume.back() - r.water_volume.front();
    const double ratio = cum_influx / dvol;
    std::printf("    cum boundary Darcy influx=%.6e  water-vol increase=%.6e  ratio=%.6f\n", cum_influx, dvol, ratio);
    check(dvol > 0.0, "(2) water content increased (infiltration)");
    check(std::fabs(ratio - 1.0) < 2e-3, "(2) global mass conservation (influx == storage change)");

    // (3) Philip early time: I(t) ∝ √t. I(t)=water-volume gain(t). t1=2, t2=8 → I2/I1 ≈ √4 = 2.
    auto I_at = [&](double t) {
        int si = 0; double best = 1e30;
        for (size_t i = 0; i < r.times.size(); ++i) if (std::fabs(r.times[i] - t) < best) { best = std::fabs(r.times[i] - t); si = (int)i; }
        return r.water_volume[si] - r.water_volume.front();
    };
    const double I1 = I_at(2.0), I2 = I_at(8.0);
    const double rt = I2 / I1, expect = std::sqrt(8.0 / 2.0);
    std::printf("  (3) Philip sqrt(t): I(8)/I(2)=%.3f  (sqrt-law %.3f; linear would be 4.0)\n", rt, expect);
    check(rt > 1.8 && rt < 2.8, "(3) early-time infiltration scales ~ sqrt(t) (not linear)");
}

}  // namespace

int main() {
    test_saturated_terzaghi();
    test_infiltration();
    if (g_failures == 0) {
        std::printf("OK: transient unsaturated (Richards) flow verified "
                    "(saturated Terzaghi + mass conservation + Philip sqrt(t))\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
