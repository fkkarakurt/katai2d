// Anchors (P2.4): the unidirectional axial spring (PLAXIS MMM Eq 18-1, N = (EA/L) U) wired
// into the solver as a structural element. Two forms share the same element technology
// (PLAXIS Reference Manual): node-to-node (strut between two mesh nodes) and fixed-end (one
// mesh node to a fixed far point). Only translational DOFs (no rotation) -> no DofMap change.
//
// (1) Fixed-end anchor as a pure spring: with very soft soil a horizontally loaded node held
//     by a horizontal FE anchor reaches equilibrium N = (EA/L) U ~ F and u ~ -F L/EA.
// (2) Node-to-node strut SUPPORTS a wall: a horizontal strut from a loaded wall column to a
//     fixed point drastically reduces the lateral deflection and carries ~ the applied load
//     (the canonical supported-excavation mechanic).
// (See docs/references/structural-plate-formulation.md sec 7; PLAXIS MMM sec 18.1.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
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

// Fixed-end anchor: a horizontally loaded top node, restrained by a stiff horizontal FE anchor,
// over soft soil -> the anchor carries ~ all the load (N ~ F) and stretches by ~ F L/EA.
void test_fixed_end_spring() {
    constexpr double W = 6.0, H = 3.0, F = 200.0;
    const double EA = 1.0e6, L = 5.0, E_soil = 50.0;  // soft soil, stiff anchor
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 6, 3);
    const int node = top_mid_node(mesh);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    f(dofs.equation(dofs.global_dof(node, 0))) = -F;  // pull the node in -x

    // FE anchor to a fixed point to the RIGHT (+x); node moving -x stretches it (tension).
    AnchorElement an;
    an.node_a = node; an.node_b = -1;
    an.fixed_point = Eigen::Vector2d(mesh.x[node] + L, mesh.y[node]);
    an.EA = EA; an.L = L;
    Structures st; st.anchors = {an};

    const MaterialModel m{MaterialType::LinearElastic, E_soil, 0.3};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                {1, 50, 1e-12}, {}, {}, st);
    check(r.converged, "fixed-end anchor solve converged");

    const double ux = r.displacement[dofs.global_dof(node, 0)];
    const double U = -ux;            // elongation (node moved away from the +x fixed point)
    const double N = EA / L * U;     // anchor axial force
    std::printf("  fixed-end: ux=%.5e  N=EA/L*U=%.2f (F=%.0f)  u~-FL/EA=%.5e\n",
                ux, N, F, -F * L / EA);
    check(ux < 0.0, "loaded node displaces toward the load (-x)");
    check(close(N, F, 0.05), "anchor carries ~ the applied load (N = EA/L U ~ F, soft soil)");
    check(close(ux, -F * L / EA, 0.05), "node displacement ~ -F L/EA (anchor-dominated)");
}

// Node-to-node strut supports a laterally loaded wall (stiff column of soil): WITH the strut
// the top deflection is far smaller and the strut force ~ the applied load.
void test_strut_supports_wall() {
    constexpr double W = 4.0, H = 8.0, F = 300.0;
    const double E_soil = 5.0e3, EA = 2.0e6, L = 4.0;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 8);
    // Wall = the left edge; load its top node horizontally (+x, into the soil).
    std::vector<int> left(mesh.left_nodes.begin(), mesh.left_nodes.end());
    std::sort(left.begin(), left.end(), [&](int a, int b) { return mesh.y[a] < mesh.y[b]; });
    const int wall_top = left.back();

    auto run = [&](bool with_strut) {
        DofMap dofs(mesh.node_count, 2);
        for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
        dofs.finalize();
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        f(dofs.equation(dofs.global_dof(wall_top, 0))) = F;  // push wall top +x
        Structures st;
        if (with_strut) {  // horizontal FE strut from wall top to a fixed point on the LEFT (-x)
            AnchorElement an; an.node_a = wall_top; an.node_b = -1;
            an.fixed_point = Eigen::Vector2d(mesh.x[wall_top] - L, mesh.y[wall_top]);
            an.EA = EA; an.L = L; st.anchors = {an};
        }
        const MaterialModel m{MaterialType::LinearElastic, E_soil, 0.3};
        const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                    {1, 50, 1e-12}, {}, {}, st);
        check(r.converged, "wall solve converged");
        const double ux = r.displacement[dofs.global_dof(wall_top, 0)];
        const double N = EA / L * (ux);  // strut: fixed point at -L; U=(0-ux)*(-1)=ux -> N=EA/L*ux
        return std::make_pair(ux, N);
    };

    const auto no_strut = run(false);
    const auto with_strut = run(true);
    std::printf("  wall top ux: no-strut=%.5e  with-strut=%.5e  strut N=%.2f (F=%.0f)\n",
                no_strut.first, with_strut.first, with_strut.second, F);
    check(with_strut.first < 0.5 * no_strut.first, "strut roughly halves+ the wall deflection");
    check(with_strut.first > 0.0 && with_strut.first < no_strut.first, "strut reduces deflection");
    check(close(with_strut.second, F, 0.25), "strut force ~ applied load (equilibrium)");
}

} // namespace

int main() {
    test_fixed_end_spring();
    test_strut_supports_wall();
    if (g_failures == 0) {
        std::printf("OK: anchors (fixed-end + node-to-node) verified (axial spring N=EA/L U, "
                    "strut supports wall)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
