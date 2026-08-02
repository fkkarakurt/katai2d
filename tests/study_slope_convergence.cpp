// Mesh convergence study for the slope SRM FoS -- NOT a regression test. Runs the
// Rocscience/Slide verification #1 slope at several mesh densities (decreasing
// target element area) and reports the factor of safety so we can see whether the
// reported FoS is mesh-converged. Results are recorded in docs/validation/.
#include <katai/analysis/strength_reduction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/mesh.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace katai;
using core::DofMap;
using core::MaterialModel;
using core::MaterialType;
using core::StrengthReductionOptions;
namespace linsolve = katai::linsolve;
using mesh::Mesh;
using mesh::Triangulation;

constexpr double kDeg = 3.14159265358979323846 / 180.0;

static double run_fos(double max_area, int load_steps, bool associated,
                      int element_order, int& nelem, double& secs) {
    const std::vector<double> px = {20.0, 70.0, 70.0, 50.0, 30.0, 20.0};
    const std::vector<double> py = {20.0, 20.0, 35.0, 35.0, 25.0, 25.0};
    const std::vector<std::array<int, 2>> segs = {
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 4}}, {{4, 5}}, {{5, 0}}};
    const Triangulation tg = mesh::quality_mesh(px, py, segs, 20.0, max_area);
    const Mesh mesh = element_order == 15 ? mesh::tri15_from_triangulation(tg, 0)
                                          : mesh::tri6_from_triangulation(tg, 0);
    nelem = mesh.element_count;

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.boundary_nodes) {
        const double x = mesh.x[n], y = mesh.y[n];
        if (std::fabs(y - 20.0) < 1e-7) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        if (std::fabs(x - 20.0) < 1e-7 || std::fabs(x - 70.0) < 1e-7) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();
    Eigen::VectorXd gravity = Eigen::VectorXd::Zero(dofs.equation_count());
    core::assemble_gravity(mesh, dofs, {20.2}, gravity);
    const double phi = 19.6 * kDeg;
    const MaterialModel base{MaterialType::MohrCoulomb, 1.0e5, 0.3, 3.0, phi,
                             associated ? phi : 0.0};
    auto solve = [](const math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(k); return s->solve(r);
    };
    StrengthReductionOptions opt;
    opt.srf_min = 0.5; opt.srf_max = 2.5; opt.bisection_iterations = 10;
    opt.newton.load_steps = load_steps; opt.newton.max_iterations = 500; opt.newton.tolerance = 1e-3;

    const auto t0 = std::chrono::steady_clock::now();
    const double fos = core::factor_of_safety(mesh, dofs, gravity, base, solve, opt);
    const auto t1 = std::chrono::steady_clock::now();
    secs = std::chrono::duration<double>(t1 - t0).count();
    return fos;
}

int main() {
    std::printf("Slope FoS (ref: Bishop 0.988, Spencer 0.987, Phase2 T6 0.997)\n");
    // Non-associated (psi=0, the benchmark) vs associated (psi=phi) mesh
    // convergence. Non-associated perfect plasticity lacks strict mesh objectivity
    // (shear-band localization without regularization); associated flow is far more
    // mesh-objective. This distinguishes a modeling property from a solver bug.
    std::printf("\n[tri6] %10s %10s %12s %12s\n", "max_area", "elements", "FoS(psi=0)", "FoS(psi=phi)");
    for (double a : {10.0, 5.0, 2.5, 1.25}) {
        int ne0 = 0, ne1 = 0; double s0 = 0.0, s1 = 0.0;
        const double f0 = run_fos(a, 8, false, 6, ne0, s0);
        const double f1 = run_fos(a, 8, true, 6, ne1, s1);
        std::printf("%10.2f %10d %12.4f %12.4f\n", a, ne0, f0, f1);
        std::fflush(stdout);
    }
    // tri15 at two densities (15-node quartic): exercises the full nonlinear MC
    // path (return mapping + consistent tangent + line search) on quartic elements.
    std::printf("\n[tri15] %9s %10s %12s %8s\n", "max_area", "elements", "FoS(psi=0)", "secs");
    for (double a : {10.0, 5.0}) {
        int ne = 0; double s = 0.0;
        const double f = run_fos(a, 8, false, 15, ne, s);
        std::printf("%10.2f %10d %12.4f %8.1f\n", a, ne, f, s);
        std::fflush(stdout);
    }
    return 0;
}
