// End-to-end tri15 pipeline verification: unstructured mesher -> tri15 enrichment
// -> generalized assembly (gravity + stiffness) -> solve. The exact geostatic field
// under self-weight with lateral confinement is sigma_yy(y) = -gamma*(H-y), linear
// in depth and hence in the tri15 (quartic) finite-element space, so a correct
// pipeline reproduces it to round-off at every Gauss point. This exercises the new
// tri15_from_triangulation enrichment and the element-generic assemble_* paths.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri15.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
using katai::mesh::Triangulation;
namespace tri15 = katai::core::tri15;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

constexpr double kW = 4.0, kH = 10.0;
constexpr double kGamma = 18.0;
constexpr double kE = 10000.0, kNu = 0.3;

void test_tri15_pipeline() {
    // 1. Mesh the rectangle, enrich to a 15-node quartic mesh.
    const std::vector<double> px = {0.0, kW, kW, 0.0};
    const std::vector<double> py = {0.0, 0.0, kH, kH};
    const std::vector<std::array<int, 2>> segs = {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}};
    const Triangulation tg = katai::mesh::quality_mesh(px, py, segs, 25.0, 3.0);
    const Mesh mesh = katai::mesh::tri15_from_triangulation(tg, 0);

    check(mesh.nodes_per_element == 15, "mesh is 15-node");
    check(mesh.element_count > 10, "non-trivial mesh");
    // Each element has 3 interior (private) nodes; shared corners/edges add the rest.
    check(mesh.node_count >= 3 * mesh.element_count, "tri15 nodes generated");

    // 2. BCs: base fixed, sides on vertical rollers.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.boundary_nodes) {
        const double x = mesh.x[n], y = mesh.y[n];
        if (std::fabs(y) < 1e-7) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        if (std::fabs(x) < 1e-7 || std::fabs(x - kW) < 1e-7) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();

    // 3. Assemble (element-generic) and solve elastically.
    const std::vector<LinearElastic> materials = {{kE, kNu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {kGamma}, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u =
        katai::core::expand_to_full(dofs, solver->solve(f));

    // 4. Recover stress directly at Gauss points: sigma = D * B * u_e. Compare to
    //    the exact geostatic field; also check positive Jacobians and u_x = 0.
    const Eigen::Matrix3d D = materials[0].plane_strain_matrix();
    double max_sigma_err = 0.0, max_ux = 0.0, min_detJ = 1e300;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri15::NodeCoords coords;
        Eigen::Matrix<double, 30, 1> ue;
        for (int k = 0; k < 15; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            ue(2 * k) = u[dofs.global_dof(n, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(n, 1)];
            max_ux = std::fmax(max_ux, std::fabs(u[dofs.global_dof(n, 0)]));
        }
        for (const auto& gp : tri15::gauss_points()) {
            const auto g = tri15::strain_displacement(coords, gp.xi, gp.eta);
            min_detJ = std::fmin(min_detJ, g.det_jacobian);
            const Eigen::Vector3d sigma = D * (g.B * ue);
            const tri15::ShapeValues N = tri15::shape_functions(gp.xi, gp.eta);
            double yq = 0.0;
            for (int k = 0; k < 15; ++k) yq += N(k) * coords(k, 1);
            const double exact = -kGamma * (kH - yq);
            max_sigma_err = std::fmax(max_sigma_err, std::fabs(sigma(1) - exact));
        }
    }
    std::printf("  tri15 pipeline: elems=%d nodes=%d  max|sig_yy err|=%.3e  "
                "max|ux|=%.3e  min detJ=%.3e\n",
                mesh.element_count, mesh.node_count, max_sigma_err, max_ux, min_detJ);
    check(min_detJ > 0.0, "all tri15 elements have positive Jacobian");
    check(max_sigma_err < 1e-6 * kGamma * kH, "tri15: sigma_yy matches -gamma*(H-y)");
    check(max_ux < 1e-9, "tri15: u_x is zero (lateral confinement)");
}

} // namespace

int main() {
    test_tri15_pipeline();
    if (g_failures == 0) {
        std::printf("OK: tri15 mesher->assembly->solve pipeline verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
