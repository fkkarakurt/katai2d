// Newton-Raphson skeleton verification (P1.1).
//
// The nonlinear solver, with a LINEAR ELASTIC material, must reproduce the direct linear
// PARDISO solution to round-off (constant tangent → ~2 iterations in a single step).
// Splitting the load into many steps (path-independent LE) must also give the same
// result. This verifies the internal-force/residual/tangent/load-step machinery
// independently of constitutive complexity.
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
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

constexpr double kE = 10000.0;
constexpr double kNu = 0.3;
constexpr double kPressure = 50.0;

// SPD linear solve (PARDISO) — handed to Newton as a callback.
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k,
                          const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(r);
}

// Laterally confined column: base clamped, side rollers (ux=0), pressure p on top.
struct Setup {
    Mesh mesh;
    DofMap dofs;
    Eigen::VectorXd f_ext;
    Setup(Mesh m, DofMap d, Eigen::VectorXd f)
        : mesh(std::move(m)), dofs(std::move(d)), f_ext(std::move(f)) {}
};

Setup build_column() {
    const RectangularDomain domain{0.0, 0.0, 4.0, 8.0, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 4);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) {
        dofs.fix_node_component(n, 0);
        dofs.fix_node_component(n, 1);
    }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0,
                                           -kPressure, f);
    return Setup(std::move(mesh), std::move(dofs), std::move(f));
}

void test_newton_matches_linear() {
    Setup s = build_column();

    // Reference: the direct linear solve.
    const std::vector<LinearElastic> le = {{kE, kNu}};
    katai::math::SparseMatrixBuilder builder(s.dofs.equation_count());
    katai::core::assemble_stiffness(s.mesh, s.dofs, le, builder);
    const Eigen::VectorXd u_ref =
        katai::core::expand_to_full(s.dofs, solve_spd(builder.build(), s.f_ext));

    // Newton (LE material model), a single load step.
    const std::vector<MaterialModel> mm = {
        {MaterialType::LinearElastic, kE, kNu}};
    const auto r1 = katai::core::solve_nonlinear(s.mesh, s.dofs, mm, s.f_ext,
                                                 solve_spd, {1, 50, 1e-12});

    check(r1.converged, "single-step Newton converged");
    check(r1.total_iterations == 2, "LE single step exactly 2 iterations");
    const double err1 = (r1.displacement - u_ref).cwiseAbs().maxCoeff();
    std::printf("  single-step: iter=%d, max|u-u_ref|=%.3e\n",
                r1.total_iterations, err1);
    check(err1 < 1e-10, "Newton(LE) matched the linear solve to round-off");

    // Multi-step load: path-independent LE → the same result.
    NewtonOptions opt;
    opt.load_steps = 4;
    opt.tolerance = 1e-12;
    const auto r4 =
        katai::core::solve_nonlinear(s.mesh, s.dofs, mm, s.f_ext, solve_spd, opt);
    check(r4.converged, "4-step Newton converged");
    const double err4 = (r4.displacement - u_ref).cwiseAbs().maxCoeff();
    std::printf("  4-step : iter=%d, max|u-u_ref|=%.3e\n",
                r4.total_iterations, err4);
    check(err4 < 1e-10, "the multi-step load gave the same LE result");
}

} // namespace

int main() {
    test_newton_matches_linear();
    if (g_failures == 0) {
        std::printf("OK: Newton-Raphson skeleton (LE) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
