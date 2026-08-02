// LBB / inf-sup DIAGNOSTIC (audit P2b): do the equal-order u-p Biot elements produce spurious
// checkerboard pore-pressure oscillations under 2D undrained (near-incompressible, t=0+) loading?
//
// consolidation.hpp uses the SAME shape functions N for both displacement and pore pressure (tri6:
// P2-P2, tri15: P4-P4). That violates the LBB/inf-sup condition; the storage block S = int (n/Kw) N^T N
// only stabilizes it by O(1/Kw), which vanishes as the pore fluid becomes incompressible. The audit
// flagged 2D undrained pore realism as UNPROVEN (all consolidation tests are 1D columns, which barely
// excite the mode, and skip the early-time boundary layer Tv<0.15).
//
// Setup: a laterally-confined, base-fixed, SEALED (undrained) saturated block with a partial-width
// strip pressure on the crest, applied at t=0+. A PARTIAL (strip) load gives a 2D non-uniform undrained
// pore field with an edge stress-concentration -- exactly where checkerboard appears.
//
// THREE probes:
//  (1) dt->0 (pure undrained saddle point [K L; L^T -S]) at realistic vs soft Kw. If the node-to-node
//      pore oscillation amplitude is much larger at the incompressible (realistic) Kw, the LBB mode is
//      real. Amplitude-normalized metric: osc[n] = |p[n] - mean(p over element-neighbours)| / range.
//  (2) dt sweep at realistic Kw: the backward-Euler diffusion Dt*H smooths high-frequency modes, so a
//      larger time step damps checkerboard. Shows whether NORMAL time-stepping hides it in practice.
//  (3) displacement (settlement) oscillation: is the SETTLEMENT field polluted too, or is the
//      checkerboard confined to the pore pressure (the pressure null-space)?  Determines practical impact.
#include <katai/analysis/consolidation.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using namespace katai::core;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {

ConsolidationSolveFactory pardiso_sym_indef() {
    return [](const katai::math::CsrMatrix& A) {
        std::shared_ptr<linsolve::DirectSolver> s =
            linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricIndefinite);
        s->factorize(A);
        return std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(
            [s](const Eigen::VectorXd& b) { return s->solve(b); });
    };
}

// Amplitude-normalized node-to-node oscillation of a nodal field: |f[n] - mean(neighbours)| / range.
// range-normalized, so smooth fields -> ~0.02-0.05 (discretization), a checkerboard -> O(0.1-1).
double node_osc_rms(const Mesh& mesh, const std::vector<std::vector<int>>& nbr,
                    const std::vector<double>& f, double* pmin = nullptr, double* pmax = nullptr,
                    int* checker = nullptr) {
    double lo = 1e300, hi = -1e300;
    for (int n = 0; n < mesh.node_count; ++n) { lo = std::min(lo, f[n]); hi = std::max(hi, f[n]); }
    const double range = std::max(1e-30, hi - lo);
    double ss = 0.0; int cnt = 0, chk = 0;
    for (int n = 0; n < mesh.node_count; ++n) {
        if (nbr[n].empty()) continue;
        double m = 0.0; for (int j : nbr[n]) m += f[j]; m /= (double)nbr[n].size();
        const double o = std::fabs(f[n] - m) / range; ss += o * o; ++cnt; if (o > 0.10) ++chk;
    }
    if (pmin) *pmin = lo; if (pmax) *pmax = hi; if (checker) *checker = chk;
    return std::sqrt(ss / std::max(1, cnt));
}

struct Result { double pore_osc, pore_min, pore_max, pore_range, disp_osc; int checker; };

// One undrained solve; returns pore + settlement oscillation at the last recorded step.
Result run(int order, double kw_over_n, double q, double dt) {
    constexpr double W = 12.0, H = 8.0, xf0 = 4.0, xf1 = 8.0;
    constexpr double E = 1.0e4, nu = 0.3, gamma_w = 10.0, k = 1.0e-3;
    Mesh mesh = (order == 15) ? katai::mesh::generate_structured_tri15({0.0, 0.0, W, H, 0}, 24, 16)
                              : katai::mesh::generate_structured_tri6({0.0, 0.0, W, H, 0}, 24, 16);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    std::vector<char> drained(mesh.node_count, 0);   // sealed -> undrained limit at t=0+
    std::vector<double> p0(mesh.node_count, 0.0);
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const std::vector<Permeability> perm = {{k, k}};

    std::vector<int> ft;
    for (int n : mesh.top_nodes) if (mesh.x[n] >= xf0 - 1e-9 && mesh.x[n] <= xf1 + 1e-9) ft.push_back(n);
    std::sort(ft.begin(), ft.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    for (size_t i = 0; i < ft.size(); ++i) {
        const double xL = (i == 0) ? mesh.x[ft[i]] : mesh.x[ft[i - 1]];
        const double xR = (i + 1 == ft.size()) ? mesh.x[ft[i]] : mesh.x[ft[i + 1]];
        const double trib = 0.5 * (mesh.x[ft[i]] - xL) + 0.5 * (xR - mesh.x[ft[i]]);
        const int eq = dofs.equation(dofs.global_dof(ft[i], 1));
        if (eq >= 0) f(eq) += -q * trib;
    }

    const auto r = solve_consolidation(mesh, dofs, mm, perm, gamma_w, kw_over_n, drained, p0, dt, 1,
                                       {}, &f, pardiso_sym_indef());

    std::vector<std::vector<int>> nbr(mesh.node_count);
    { std::vector<std::vector<char>> seen(mesh.node_count);
      for (int e = 0; e < mesh.element_count; ++e)
        for (int a = 0; a < mesh.nodes_per_element; ++a) for (int b = 0; b < mesh.nodes_per_element; ++b)
            if (a != b) nbr[mesh.node_of(e, a)].push_back(mesh.node_of(e, b));
      for (auto& v : nbr) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); } }

    Result R{};
    std::vector<double> pore(r.pore.back().data(), r.pore.back().data() + mesh.node_count);
    R.pore_osc = node_osc_rms(mesh, nbr, pore, &R.pore_min, &R.pore_max, &R.checker);
    R.pore_range = R.pore_max - R.pore_min;
    std::vector<double> uy(mesh.node_count);
    for (int n = 0; n < mesh.node_count; ++n) uy[n] = r.displacement.back()[dofs.global_dof(n, 1)];
    R.disp_osc = node_osc_rms(mesh, nbr, uy);
    return R;
}

void report(int order) {
    const char* nm = (order == 15) ? "tri15 (P4-P4)" : "tri6 (P2-P2)";
    constexpr double q = 100.0;
    const double kw_real = 2.0e6 / 0.3, kw_soft = 1.0e2;
    std::printf("\n=== %s, strip footing q=%.0f kPa, undrained ===\n", nm, q);
    const Result a = run(order, kw_real, q, 1e-8);
    const Result b = run(order, kw_soft, q, 1e-8);
    std::printf("  (1) t=0+ (dt->0):  realistic Kw/n=%.2e   pore rms-osc=%.3f  checker-nodes=%d  pore[%.1f..%.1f]  disp rms-osc=%.3f\n",
                kw_real, a.pore_osc, a.checker, a.pore_min, a.pore_max, a.disp_osc);
    std::printf("                     soft Kw/n=%.2e (control)   pore rms-osc=%.3f  checker-nodes=%d  disp rms-osc=%.3f\n",
                kw_soft, b.pore_osc, b.checker, b.disp_osc);
    const double ratio = a.pore_osc / std::max(1e-9, b.pore_osc);
    std::printf("      pore-osc amplitude ratio (incompressible / compressible) = %.1fx\n", ratio);
    // (2) dt sweep at realistic Kw: does backward-Euler diffusion damp the checkerboard?
    std::printf("  (2) dt sweep at realistic Kw (Dt*H diffusion damping):\n");
    for (double dt : {1e-8, 1e-2, 1.0, 100.0}) {
        const Result s = run(order, kw_real, q, dt);
        std::printf("        dt=%9.1e  pore rms-osc=%.3f  checker-nodes=%d\n", dt, s.pore_osc, s.checker);
    }
    // Verdict from the amplitude-normalized metric at the undrained instant.
    const bool lbb = a.pore_osc > 0.06 && ratio > 2.0;
    std::printf("  VERDICT: %s\n", lbb
        ? "CHECKERBOARD pore field at the incompressible undrained limit -> LBB instability is REAL."
        : "pore field smooth at the incompressible limit -> equal-order adequate here.");
    if (lbb) std::printf("           (settlement rms-osc=%.3f -> displacements %s.)\n",
        a.disp_osc, a.disp_osc > 0.06 ? "ALSO polluted" : "stay smooth; artefact confined to pore pressure");
}

}  // namespace

int main() {
    std::printf("LBB / inf-sup diagnostic for equal-order u-p Biot consolidation (audit P2b)\n");
    report(6);
    report(15);
    std::printf("\n(rms-osc: RMS of |p_node - mean(element-neighbours)| / range; checkerboard alternates node-to-node.)\n");
    return 0;
}
