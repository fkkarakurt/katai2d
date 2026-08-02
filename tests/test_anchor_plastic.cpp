// Anchor elastoplasticity (P2.4): the unidirectional axial spring with a force capacity
// F_max,tens / F_max,comp (PLAXIS MMM sec 18.1). Elastic-perfectly-plastic: N grows with EA/L
// stiffness until |N| reaches the cap, then plateaus while a permanent plastic elongation U_p
// accumulates (a support that yields -- the critical supported-excavation behaviour).
//
// CLOSED-FORM validation via a controlled TWO-SPRING system. One free DOF (ux of a single node;
// every other DOF fixed, soil made negligibly soft) loaded by F, held by two horizontal fixed-
// end anchors: A1 (elastoplastic, stiffness k1, cap Fc) toward +x, A2 (elastic, k2) toward -x.
// Solver sign convention: anchor elongation U = -ux for a +x fixed point. Equilibrium in x:
//   -N1 + k2 ux = f_ext_x.
// Pull (-x, f_ext_x=-F):  elastic N1 = k1 F/(k1+k2).  If that exceeds Fc, A1 yields:
//   N1 = Fc,  ux = (Fc - F)/k2,  U_p = -ux - Fc/k1   (return-mapping identity, exact).
// Push (+x, f_ext_x=+F):  A1 compresses; cap Fcc -> N1 = -Fcc, ux = (F - Fcc)/k2, U_p = -ux + Fcc/k1.
// (See docs/references/structural-plate-formulation.md sec 7a.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::AnchorElement;
using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Structures;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(k);
    return s->solve(r);
}
int top_mid_node(const Mesh& mesh) {
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    return top[top.size() / 2];
}

// One free DOF (ux of `node`), two horizontal fixed-end anchors. Returns (ux, U_p of A1).
// `Fc1` is A1's cap (Fmax_tens if tension>0 case, applied to both directions via the fields);
// `push` selects the load direction (+x if true). Soil made negligibly soft.
struct Result2 { double ux, U1p, N1; bool converged; };
Result2 run_two_spring(double EA1, double L1, double EA2, double L2, double F, double cap,
                       bool push) {
    constexpr double W = 6.0, H = 3.0;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 6, 3);
    const int node = top_mid_node(mesh);

    // Fix every DOF except ux of `node`: a clean single-DOF system dominated by the anchors.
    DofMap dofs(mesh.node_count, 2);
    for (int n = 0; n < mesh.node_count; ++n) {
        dofs.fix_node_component(n, 1);                 // all y fixed
        if (n != node) dofs.fix_node_component(n, 0);  // all x fixed but the chosen node
    }
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const int eq = dofs.equation(dofs.global_dof(node, 0));
    f(eq) = push ? F : -F;

    AnchorElement a1;  // elastoplastic, fixed point toward +x
    a1.node_a = node; a1.node_b = -1;
    a1.fixed_point = Eigen::Vector2d(mesh.x[node] + L1, mesh.y[node]);
    a1.EA = EA1; a1.L = L1;
    a1.Fmax_tens = cap; a1.Fmax_comp = cap;  // symmetric cap (only the active side matters)

    AnchorElement a2;  // elastic, fixed point toward -x
    a2.node_a = node; a2.node_b = -1;
    a2.fixed_point = Eigen::Vector2d(mesh.x[node] - L2, mesh.y[node]);
    a2.EA = EA2; a2.L = L2;  // Fmax default <=0 -> unlimited (elastic)

    Structures st; st.anchors = {a1, a2};
    const MaterialModel m{MaterialType::LinearElastic, 1.0e-3, 0.3};  // negligible soil
    const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                {4, 50, 1e-12}, {}, {}, st);
    const double ux = r.displacement[dofs.global_dof(node, 0)];
    const double U1 = -ux, U1p = r.anchor_plastic.empty() ? 0.0 : r.anchor_plastic[0];
    const double N1 = (EA1 / L1) * (U1 - U1p);
    return {ux, U1p, N1, r.converged};
}

// Tension yield of A1 (closed form).
void test_anchor_tension_cap() {
    const double EA1 = 1.0e6, L1 = 5.0, k1 = EA1 / L1;   // 2e5
    const double EA2 = 1.0e6, L2 = 10.0, k2 = EA2 / L2;  // 1e5
    const double F = 300.0;
    const double N1_elastic = k1 * F / (k1 + k2);        // 200
    const double Fc = 120.0;                             // < 200 -> yields

    // Elastic reference (cap above yield -> no plasticity).
    const auto el = run_two_spring(EA1, L1, EA2, L2, F, /*cap*/ 1.0e9, /*push*/ false);
    check(el.converged, "tension elastic-reference solve converged");
    check(close(el.ux, -F / (k1 + k2), 1e-3), "elastic ux = -F/(k1+k2)");
    check(close(el.N1, N1_elastic, 1e-3), "elastic N1 = k1 F/(k1+k2)");
    check(std::fabs(el.U1p) < 1e-9, "elastic run leaves no plastic elongation");

    // Capped: A1 yields at Fc.
    const auto pl = run_two_spring(EA1, L1, EA2, L2, F, Fc, /*push*/ false);
    const double ux_exact = (Fc - F) / k2;          // -1.8e-3
    const double U1p_exact = -ux_exact - Fc / k1;   // 1.2e-3
    std::printf("  tension: ux=%.6e (exact %.6e)  N1=%.3f (Fc=%.0f)  U_p=%.6e (exact %.6e)\n",
                pl.ux, ux_exact, pl.N1, Fc, pl.U1p, U1p_exact);
    check(pl.converged, "tension capped solve converged");
    check(close(pl.N1, Fc, 1e-3), "anchor force capped at F_max,tens");
    check(close(pl.ux, ux_exact, 2e-3), "capped ux = (Fc - F)/k2 (closed form)");
    check(close(pl.U1p, U1p_exact, 2e-3), "permanent elongation U_p = -ux - Fc/k1");
    check(std::fabs(pl.ux) > std::fabs(el.ux), "yielding anchor -> larger displacement than elastic");
}

// Compression yield of A1 (symmetric closed form).
void test_anchor_compression_cap() {
    const double EA1 = 1.0e6, L1 = 5.0, k1 = EA1 / L1;
    const double EA2 = 1.0e6, L2 = 10.0, k2 = EA2 / L2;
    const double F = 300.0;
    const double Fcc = 120.0;

    const auto pl = run_two_spring(EA1, L1, EA2, L2, F, Fcc, /*push*/ true);
    const double ux_exact = (F - Fcc) / k2;         // 1.8e-3
    const double U1p_exact = -ux_exact + Fcc / k1;  // -1.2e-3
    std::printf("  compression: ux=%.6e (exact %.6e)  N1=%.3f (-Fcc=%.0f)  U_p=%.6e (exact %.6e)\n",
                pl.ux, ux_exact, pl.N1, -Fcc, pl.U1p, U1p_exact);
    check(pl.converged, "compression capped solve converged");
    check(close(pl.N1, -Fcc, 1e-3), "anchor force capped at -F_max,comp");
    check(close(pl.ux, ux_exact, 2e-3), "capped ux = (F - Fcc)/k2 (closed form)");
    check(close(pl.U1p, U1p_exact, 2e-3), "permanent elongation U_p = -ux + Fcc/k1");
}

} // namespace

int main() {
    test_anchor_tension_cap();
    test_anchor_compression_cap();
    if (g_failures == 0) {
        std::printf("OK: anchor elastoplasticity verified (F_max cap + permanent U_p, "
                    "two-spring closed form)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
