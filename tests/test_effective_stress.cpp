// Effective stress with a prescribed (hydrostatic) pore pressure field (P2 / water).
// A saturated column with the water table at the top surface, under its saturated
// self-weight: by Terzaghi the effective vertical stress is the BUOYANT (submerged)
// geostatic stress
//     sigma'_v(y) = -gamma' (H - y),   gamma' = gamma_sat - gamma_w,
// and the pore pressure is u(y) = gamma_w (H - y). The skeleton (effective stress)
// responds to gravity(gamma_sat) + the pore-pressure load f_pore = integral B^T u m;
// the computed stress is then effective. The field is linear -> reproduced to
// round-off. Exercises assemble_pore_pressure_load and the effective-stress framing.
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

void test_submerged_column() {
    constexpr double W = 4.0, H = 10.0;
    constexpr double gamma_sat = 20.0, gamma_w = 9.81;
    constexpr double gamma_eff = gamma_sat - gamma_w;  // buoyant unit weight
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

    // Total self-weight (saturated) + hydrostatic pore-pressure load (water table
    // at the top, u = gamma_w (H - y)).
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {gamma_sat}, f);
    auto pore = [&](double, double y) { return gamma_w * std::max(0.0, H - y); };
    katai::core::assemble_pore_pressure_load(mesh, dofs, pore, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, solver->solve(f));

    // Gauss-point effective stress sigma' = D B u_e (the load includes the pore term,
    // so the computed stress is effective). Compare sigma'_v to -gamma'(H - y).
    const Eigen::Matrix3d D = materials[0].plane_strain_matrix();
    double max_err = 0.0, max_ux = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri6::NodeCoords coords;
        Eigen::Matrix<double, 12, 1> ue;
        for (int k = 0; k < 6; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            ue(2 * k) = u[dofs.global_dof(n, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(n, 1)];
            max_ux = std::fmax(max_ux, std::fabs(u[dofs.global_dof(n, 0)]));
        }
        for (const auto& gp : tri6::gauss_points()) {
            const auto g = tri6::strain_displacement(coords, gp.xi, gp.eta);
            const Eigen::Vector3d sigma_eff = D * (g.B * ue);
            const tri6::ShapeValues N = tri6::shape_functions(gp.xi, gp.eta);
            double yq = 0.0;
            for (int k = 0; k < 6; ++k) yq += N(k) * coords(k, 1);
            const double exact = -gamma_eff * (H - yq);  // buoyant effective stress
            max_err = std::fmax(max_err, std::fabs(sigma_eff(1) - exact));
        }
    }
    std::printf("  submerged column: max|sigma'_v err|=%.3e  max|ux|=%.3e  (gamma'=%.2f)\n",
                max_err, max_ux, gamma_eff);
    check(max_err < 1e-6 * gamma_eff * H,
          "effective sigma'_v matches buoyant -gamma'(H-y)");
    check(max_ux < 1e-9, "u_x is zero (lateral confinement)");
}

} // namespace

int main() {
    test_submerged_column();
    if (g_failures == 0) {
        std::printf("OK: effective stress with prescribed pore pressure verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
