// Axisymmetric Mohr-Coulomb global validation (P1.6): fully-plastic collapse of a
// thick-walled cylinder under internal pressure. For a Tresca material (phi = 0,
// cohesion c) the limit pressure has the closed form
//     p_collapse = 2 c ln(b/a)
// (plastic zone spreads from r=a to r=b; integrate dsigma_r/dr = 2c/r with
// sigma_theta - sigma_r = 2c). We ramp the internal pressure by axisymmetric limit
// analysis and recover p_collapse from the highest equilibrated load. Exercises the
// axisymmetric MC material point + the kinematics-generalized nonlinear solver.
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::Kinematics;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::NewtonOptions;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_cylinder_collapse() {
    constexpr double a = 1.0, b = 2.0, Lz = 0.5, c = 10.0;
    const RectangularDomain domain{a, 0.0, b - a, Lz, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri15(domain, 12, 1);

    // u_z = 0 on z-faces (eps_z = 0, long cylinder); u_r free.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) dofs.fix_node_component(n, 1);
    for (int n : mesh.top_nodes) dofs.fix_node_component(n, 1);
    dofs.finalize();

    // Weightless Tresca (phi = psi = 0): associated -> symmetric tangent.
    const std::vector<MaterialModel> materials = {
        {MaterialType::MohrCoulomb, 1.0e5, 0.3, c, 0.0, 0.0}};

    // Unit internal pressure on the inner wall (r = a, left edge): t_r = +1.
    Eigen::VectorXd f_unit = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_axisym_traction(mesh, dofs, mesh.left_nodes, 1.0, 0.0, f_unit);

    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricIndefinite);
        s->factorize(k);
        return s->solve(r);
    };

    const double p_exact = 2.0 * c * std::log(b / a);  // ~13.86
    const double p_max = 1.3 * p_exact;
    NewtonOptions opt;
    opt.kinematics = Kinematics::Axisymmetric;
    opt.load_steps = 30;
    opt.max_iterations = 100;
    opt.tolerance = 1e-5;

    const auto result = katai::core::solve_nonlinear(mesh, dofs, materials,
                                                     p_max * f_unit, solve, opt);
    const double p_fem = result.load_factor * p_max;
    const double rel_err = std::fabs(p_fem - p_exact) / p_exact;
    std::printf("  axisym cylinder collapse: p = %.3f vs 2c ln(b/a) = %.3f, hata = %.1f%%\n",
                p_fem, p_exact, 100.0 * rel_err);
    check(rel_err < 0.05, "axisym MC cylinder collapse within 5% of 2c ln(b/a)");
}

} // namespace

int main() {
    test_cylinder_collapse();
    if (g_failures == 0) {
        std::printf("OK: axisymmetric MC cylinder collapse verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
