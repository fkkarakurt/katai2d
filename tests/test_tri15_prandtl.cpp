// Prandtl bearing capacity with 15-node quartic elements -- the end-to-end tri15
// demonstration (footing surface traction + Mohr-Coulomb plasticity + limit
// analysis) and the payoff of the higher-order element: a COARSE tri15 mesh
// matches N_c = 2 + pi as well as a far finer tri6 mesh. Weightless, purely
// cohesive (Tresca) half-space; q_ult from incremental limit analysis.
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
constexpr double kPi = 3.14159265358979323846;

void test_tri15_prandtl() {
    // Coarse tri15 mesh: cell 0.6, footing B = 1.2 spans 2 cells (8 quartic edges).
    constexpr double W = 6.0, H = 4.0, x0 = 2.4, x1 = 3.6;
    constexpr int nx = 10, ny = 6;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri15(domain, nx, ny);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    constexpr double c = 10.0;
    const std::vector<MaterialModel> materials = {
        {MaterialType::MohrCoulomb, 10000.0, 0.3, c, 0.0, 0.0}};  // Tresca

    // Footing = ordered top nodes within [x0, x1].
    std::vector<int> footing;
    for (int n : mesh.top_nodes) {
        const double x = mesh.x[n];
        if (x >= x0 - 1e-9 && x <= x1 + 1e-9) footing.push_back(n);
    }
    Eigen::VectorXd f_unit = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, footing, 0.0, -1.0, f_unit);

    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricIndefinite);
        s->factorize(k);
        return s->solve(r);
    };
    constexpr double q_max = 6.0 * c;
    NewtonOptions opt;
    opt.load_steps = 30;
    opt.max_iterations = 100;
    opt.tolerance = 1e-5;

    const auto result = katai::core::solve_nonlinear(mesh, dofs, materials,
                                                     q_max * f_unit, solve, opt);
    const double nc = result.load_factor * q_max / c;
    const double nc_exact = 2.0 + kPi;
    const double rel_err = std::fabs(nc - nc_exact) / nc_exact;
    std::printf("  tri15 Prandtl: elems=%d eqs=%d  N_c = %.3f vs %.3f, hata = %.1f%%\n",
                mesh.element_count, dofs.equation_count(), nc, nc_exact, 100.0 * rel_err);
    check(rel_err < 0.05, "tri15 Prandtl N_c within 5% of 2+pi");
}

} // namespace

int main() {
    test_tri15_prandtl();
    if (g_failures == 0) {
        std::printf("OK: tri15 Prandtl bearing-capacity benchmark passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
