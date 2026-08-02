// Mesh node-splitting for an embedded wall (barrier) -- the mesh surgery that lets a
// wall separate the soil into two sides (active behind / passive in front), connected
// only below the toe. split_mesh_at_wall duplicates the seam nodes above the toe and
// rebinds the left-side elements to the duplicates. Validation: load the LEFT region
// horizontally; a seam node on the right and its left duplicate (same location, above
// the toe) must move DIFFERENTLY (they are disconnected -> barrier), while below the toe
// the soil stays continuous. If the split failed (shared nodes) they would move alike.
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::core::SeamPair;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_wall_split() {
    constexpr double W = 4.0, H = 4.0, xw = 2.0, y_toe = 2.0, P = 100.0, E = 1.0e4, nu = 0.3;
    RectangularDomain domain{0.0, 0.0, W, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 8, 8);  // xw=2, y_toe=2 on grid lines

    const int n_before = mesh.node_count;
    const std::vector<SeamPair> seam =
        katai::core::split_mesh_at_wall(mesh, xw, y_toe, H);  // split the upper half (y in (2,4])
    std::printf("  split: nodes %d -> %d  (%zu seam pairs)\n", n_before, mesh.node_count, seam.size());
    check(mesh.node_count > n_before && !seam.empty(), "split added duplicate seam nodes");

    // Top-left node (0, H) -- push the left region toward +x (into the wall line).
    int top_left = -1;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n]) < 1e-9 && std::fabs(mesh.y[n] - H) < 1e-9) top_left = n;

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.right_nodes)  dofs.fix_node_component(n, 0);  // far-field roller (right region)
    dofs.finalize();

    const std::vector<LinearElastic> mats = {{E, nu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, mats, builder);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    f(dofs.equation(dofs.global_dof(top_left, 0))) = P;  // horizontal load on the left region

    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, s->solve(f));

    // Pick a seam pair near mid-height (y ~ 3) and compare the right vs left-duplicate ux.
    double best = 1e9; SeamPair sp{-1, -1, 0};
    for (const auto& p : seam) if (std::fabs(p.y - 3.0) < best) { best = std::fabs(p.y - 3.0); sp = p; }
    const double ux_right = u(dofs.global_dof(sp.right, 0));
    const double ux_left  = u(dofs.global_dof(sp.left, 0));
    std::printf("  seam @y=%.2f: ux_right=%.4e  ux_left(dup)=%.4e  gap=%.4e\n",
                sp.y, ux_right, ux_left, ux_left - ux_right);
    check(std::fabs(ux_left - ux_right) > 0.2 * std::fabs(ux_left),
          "seam right and left-duplicate move differently -> wall is a barrier (split works)");
    check(ux_left > ux_right, "loaded left region moves more than the shielded right region");
}

} // namespace

int main() {
    test_wall_split();
    if (g_failures == 0) {
        std::printf("OK: embedded-wall mesh node-splitting (barrier) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
