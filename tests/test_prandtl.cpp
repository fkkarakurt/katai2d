// Prandtl bearing-capacity benchmark (P1.3) -- the first global validation of
// Mohr-Coulomb plasticity.
//
// A weightless, purely cohesive (phi = 0, Tresca) half-space under a flexible
// strip footing collapses at the bearing pressure
//     q_ult = (2 + pi) * c  =>  N_c = 2 + pi ~= 5.14   (Prandtl, 1921).
// We obtain q_ult by incremental limit analysis: the footing pressure is ramped
// in equal increments and the highest level the Newton-Raphson scheme can
// equilibrate (NewtonResult::load_factor) approaches the collapse load from
// below. The elastoplastic tangent loses positive-definiteness near collapse, so
// the linear solver uses the real-symmetric-indefinite path.
//
// N_c is approached from slightly above (a stiffer, under-resolved footing-edge
// singularity on the uniform mesh). Even so the structured mesh already attains
// ~1% error, comfortably inside the project's <5% accuracy target; the graded
// mesher (P1.4) will tighten it further.
// verify: KV-FND-005
//   oracle:   closed_form
//   source:   Prandtl (1921) bearing-capacity wedge solution for a weightless Tresca (phi = 0) half-space
//   locator:  q_ult = (2 + pi) c, i.e. N_c = 2 + pi = 5.14159 (stated in full above)
//   quantity: bearing-capacity factor N_c [-], incremental limit analysis approaching collapse from below
//   expected: 5.142
//   band:     5%, as asserted below -- the uniform structured mesh under-resolves the footing-edge singularity; tri15 at 120 elements measures +1.1% (test_tri15_prandtl)
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

void test_prandtl() {
    // --- Geometry and structured mesh ---------------------------------------
    constexpr double W = 6.0, H = 4.0;       // domain
    constexpr double x0 = 2.4, x1 = 3.6;     // footing extent (width B = 1.2)
    constexpr int nx = 30, ny = 20;          // cell size 0.2 (6 cells under footing)
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, ny);

    // --- Boundary conditions: base fixed, sides on vertical rollers ---------
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) {
        dofs.fix_node_component(n, 0);
        dofs.fix_node_component(n, 1);
    }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // --- Material: weightless, cohesive, frictionless (Tresca) --------------
    constexpr double c = 10.0;
    const std::vector<MaterialModel> materials = {
        {MaterialType::MohrCoulomb, 10000.0, 0.3, c, 0.0, 0.0}};

    // --- Unit footing pressure (flexible, uniform traction) -----------------
    std::vector<int> footing;
    for (int n : mesh.top_nodes) {
        const double x = mesh.x[n];
        if (x >= x0 - 1e-9 && x <= x1 + 1e-9) footing.push_back(n);
    }
    Eigen::VectorXd f_unit = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, footing, 0.0, -1.0, f_unit);

    // --- Incremental limit analysis -----------------------------------------
    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricIndefinite);
        s->factorize(k);
        return s->solve(r);
    };

    constexpr double q_max = 6.0 * c;  // ramp just past collapse (N_c ~ 5.14)
    NewtonOptions opt;
    // Adaptive load stepping refines automatically near collapse, so a modest
    // initial step count suffices (and is faster).
    opt.load_steps = 40;
    opt.max_iterations = 100;
    opt.tolerance = 1e-6;
    const Eigen::VectorXd f_ext = q_max * f_unit;

    const auto result =
        katai::core::solve_nonlinear(mesh, dofs, materials, f_ext, solve, opt);

    const double q_ult = result.load_factor * q_max;
    const double nc = q_ult / c;
    const double nc_exact = 2.0 + kPi;
    const double rel_err = std::fabs(nc - nc_exact) / nc_exact;
    std::printf("  Prandtl: N_c(FEM) = %.3f, exact = %.3f (2+pi), hata = %.1f%%\n",
                nc, nc_exact, 100.0 * rel_err);

    // Project accuracy target: < 5% deviation on critical results.
    check(rel_err < 0.05, "Prandtl N_c within 5% of 2+pi");
}

} // namespace

int main() {
    test_prandtl();
    if (g_failures == 0) {
        std::printf("OK: Prandtl bearing-capacity benchmark passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
