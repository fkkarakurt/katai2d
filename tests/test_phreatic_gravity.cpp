// Phreatic (water-table-aware) self-weight: assemble_gravity_phreatic picks the saturated
// unit weight below the phreatic surface and the moist unit weight above, per Gauss point.
// Combined with the hydrostatic pore-pressure load the recovered stress is EFFECTIVE.
//
// A confined column of height H with the water table partway up at y_w. By Terzaghi the
// effective vertical stress is piecewise linear:
//     above (y >= y_w):  sigma'_v = -gamma_unsat (H - y)
//     below (y <  y_w):  sigma'_v = -[gamma_unsat (H - y_w) + gamma' (y_w - y)],  gamma' = gamma_sat - gamma_w
// The water table is placed on a mesh node row so no element straddles it -> the field is
// piecewise linear in FE space and reproduced to round-off. Validates assemble_gravity_phreatic.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri6.hpp>
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
namespace tri6 = katai::core::tri6;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_partial_water_table() {
    constexpr double W = 4.0, H = 10.0, y_w = 5.0;   // water table on the k=4 node row (ny=8)
    constexpr double gamma_unsat = 17.0, gamma_sat = 20.0, gamma_w = 9.81;
    constexpr double gamma_eff = gamma_sat - gamma_w;
    constexpr double E = 1.0e4, nu = 0.3;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 3, 8);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{E, nu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto wt = [&](double) { return y_w; };
    katai::core::assemble_gravity_phreatic(mesh, dofs, {gamma_unsat}, {gamma_sat}, wt, f);
    auto pore = [&](double, double y) { return gamma_w * std::max(0.0, y_w - y); };
    katai::core::assemble_pore_pressure_load(mesh, dofs, pore, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, solver->solve(f));

    const Eigen::Matrix3d D = materials[0].plane_strain_matrix();
    double max_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri6::NodeCoords coords;
        Eigen::Matrix<double, 12, 1> ue;
        for (int k = 0; k < 6; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n]; coords(k, 1) = mesh.y[n];
            ue(2 * k) = u[dofs.global_dof(n, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(n, 1)];
        }
        for (const auto& gp : tri6::gauss_points()) {
            const auto g = tri6::strain_displacement(coords, gp.xi, gp.eta);
            const Eigen::Vector3d sigma_eff = D * (g.B * ue);
            const tri6::ShapeValues N = tri6::shape_functions(gp.xi, gp.eta);
            double yq = 0.0;
            for (int k = 0; k < 6; ++k) yq += N(k) * coords(k, 1);
            const double exact = yq >= y_w
                ? -gamma_unsat * (H - yq)
                : -(gamma_unsat * (H - y_w) + gamma_eff * (y_w - yq));
            max_err = std::fmax(max_err, std::fabs(sigma_eff(1) - exact));
        }
    }
    std::printf("  partial water table: max|sigma'_v err|=%.3e  (gamma'=%.2f)\n", max_err, gamma_eff);
    check(max_err < 1e-6 * gamma_sat * H, "effective sigma'_v matches piecewise buoyant profile");
}

} // namespace

int main() {
    test_partial_water_table();
    if (g_failures == 0) {
        std::printf("OK: phreatic self-weight effective stress verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
