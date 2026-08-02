// tri15 surface-traction verification on a structured quartic mesh. A confined
// column under a uniform top pressure p has, by vertical equilibrium (no body
// force), sigma_yy = -p everywhere -- a field exactly representable by tri15, so a
// correct quartic-edge traction + assembly reproduces it to round-off. Exercises
// generate_structured_tri15 and the element-order-aware assemble_surface_traction
// (5-node quartic edges).
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri15.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace tri15 = katai::core::tri15;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_tri15_traction() {
    constexpr double W = 4.0, H = 6.0, p = 50.0, E = 10000.0, nu = 0.3;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri15(domain, 3, 4);
    check(mesh.nodes_per_element == 15, "structured tri15 mesh");

    // Base fixed, sides on vertical rollers, top loaded with pressure p.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{E, nu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0, -p, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, solver->solve(f));

    // Consistent traction: total vertical nodal load equals p * W.
    check(std::fabs(f.sum() + p * W) < 1e-9 * p * W,
          "tri15 traction: total load = p*W (consistent integral)");

    // sigma_yy = D * B * u at Gauss points must equal -p everywhere.
    const Eigen::Matrix3d D = materials[0].plane_strain_matrix();
    double max_err = 0.0, min_detJ = 1e300;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri15::NodeCoords coords;
        Eigen::Matrix<double, 30, 1> ue;
        for (int k = 0; k < 15; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            ue(2 * k) = u[dofs.global_dof(n, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(n, 1)];
        }
        for (const auto& gp : tri15::gauss_points()) {
            const auto g = tri15::strain_displacement(coords, gp.xi, gp.eta);
            min_detJ = std::fmin(min_detJ, g.det_jacobian);
            const Eigen::Vector3d sigma = D * (g.B * ue);
            max_err = std::fmax(max_err, std::fabs(sigma(1) + p));
        }
    }
    std::printf("  tri15 traction: elems=%d nodes=%d  max|sig_yy+p|=%.3e  min detJ=%.3e\n",
                mesh.element_count, mesh.node_count, max_err, min_detJ);
    check(min_detJ > 0.0, "tri15 traction: positive Jacobians");
    check(max_err < 1e-6 * p, "tri15 traction: sigma_yy = -p (uniform, round-off)");
}

} // namespace

int main() {
    test_tri15_traction();
    if (g_failures == 0) {
        std::printf("OK: tri15 surface traction verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
