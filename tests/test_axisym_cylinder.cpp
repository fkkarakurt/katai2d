// Axisymmetric global validation (P1.6): a thick-walled cylinder under internal
// pressure, against the closed-form Lamé solution. A long cylinder (eps_z = 0,
// plane strain in z) with inner radius a, outer b, internal pressure p_i:
//     sigma_r(r)     = (p_i a^2 / (b^2 - a^2)) (1 - b^2/r^2)
//     sigma_theta(r) = (p_i a^2 / (b^2 - a^2)) (1 + b^2/r^2)
// (independent of E, nu). We model the annulus a<=r<=b, 0<=z<=Lz with u_z = 0 on
// the z-faces (eps_z = 0), internal pressure on r=a, free outer wall, and compare
// the Gauss-point sigma_r, sigma_theta. Exercises the axisymmetric element + the
// r-weighted stiffness and traction assembly end-to-end.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::core::Tri15Element;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace axisym = katai::core::axisym;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_cylinder() {
    constexpr double a = 1.0, b = 2.0, Lz = 1.0, pi = 10.0, E = 10000.0, nu = 0.3;
    const RectangularDomain domain{a, 0.0, b - a, Lz, 0};  // (r, z) rectangle
    const Mesh mesh = katai::mesh::generate_structured_tri15(domain, 8, 2);

    // u_z = 0 on top and bottom (eps_z = 0, long cylinder); u_r free.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) dofs.fix_node_component(n, 1);
    for (int n : mesh.top_nodes) dofs.fix_node_component(n, 1);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{E, nu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_axisym_stiffness(mesh, dofs, materials, builder);

    // Internal pressure p_i on the inner wall (r = a, the left edge): traction in
    // +r pushes the solid outward, giving sigma_r(a) = -p_i.
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_axisym_traction(mesh, dofs, mesh.left_nodes, pi, 0.0, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, solver->solve(f));

    const Eigen::Matrix4d D = materials[0].axisymmetric_matrix();
    const double coef = pi * a * a / (b * b - a * a);
    double max_err = 0.0, scale = pi;
    for (int e = 0; e < mesh.element_count; ++e) {
        Tri15Element::NodeCoords coords;
        Eigen::Matrix<double, 30, 1> ue;
        for (int k = 0; k < 15; ++k) {
            const int nidx = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[nidx];
            coords(k, 1) = mesh.y[nidx];
            ue(2 * k) = u[dofs.global_dof(nidx, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(nidx, 1)];
        }
        for (const auto& gp : Tri15Element::gauss_points()) {
            const auto g = axisym::strain_displacement<Tri15Element>(coords, gp.xi, gp.eta);
            const Eigen::Vector4d sigma = D * (g.B * ue);  // [sr, sz, srz, stheta]
            const double r = g.radius;
            const double sr_exact = coef * (1.0 - b * b / (r * r));
            const double sth_exact = coef * (1.0 + b * b / (r * r));
            max_err = std::fmax(max_err, std::fabs(sigma(0) - sr_exact));
            max_err = std::fmax(max_err, std::fabs(sigma(3) - sth_exact));
        }
    }
    const double rel = max_err / scale;
    std::printf("  axisym cylinder: elems=%d  max|sigma err|/p_i = %.3e\n",
                mesh.element_count, rel);
    check(rel < 0.02, "axisym Lame cylinder: sigma_r, sigma_theta within 2% of p_i");
}

} // namespace

int main() {
    test_cylinder();
    if (g_failures == 0) {
        std::printf("OK: axisymmetric thick-cylinder (Lame) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
