// Bearing capacity of a cohesive-frictional (phi > 0) soil -- the global
// validation of Mohr-Coulomb plasticity with friction, complementing the phi = 0
// Prandtl benchmark (P1.3) and exercising the analytic consistent tangent + line
// search on a second, independent closed-form solution.
//
// For a weightless soil with cohesion c and friction angle phi under a flexible
// strip footing, the ultimate bearing pressure is (Prandtl-Reissner, associated
// flow):
//     q_ult = c * Nc,   Nc = (Nq - 1) cot(phi),   Nq = e^(pi tan phi) tan^2(45 + phi/2).
// For phi = 20 deg this gives Nq = 6.400, Nc = 14.835 (standard tables). We obtain
// q_ult by incremental limit analysis (the highest equilibrated load factor) and
// require < 5% error. Associated flow (psi = phi) matches the exact slip-line
// solution and yields a symmetric tangent -> the symmetric-indefinite solver path.
//
// As in Prandtl, the footing-edge stress singularity on a uniform mesh makes the
// FE capacity approach the exact value from slightly above; the domain is sized to
// contain the (wider, deeper) phi > 0 collapse mechanism so it is not truncated.
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
constexpr double kDeg = kPi / 180.0;

void test_bearing_phi() {
    // --- Geometry sized to contain the phi=20 deg mechanism (W/B~6.7, H/B~3.3) --
    constexpr double W = 8.0, H = 4.0;        // domain
    constexpr double B = 1.2;                 // footing width
    constexpr double x0 = 0.5 * (W - B), x1 = 0.5 * (W + B);  // centered footing
    constexpr int nx = 40, ny = 20;           // cell size 0.2 (6 cells under footing)
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, ny);

    // --- Base fixed, sides on vertical rollers -------------------------------
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // --- Weightless cohesive-frictional soil, associated flow (psi = phi) -----
    constexpr double c = 10.0;
    constexpr double phi = 20.0 * kDeg;
    const std::vector<MaterialModel> materials = {
        {MaterialType::MohrCoulomb, 1.0e5, 0.3, c, phi, phi}};

    // --- Unit footing pressure (flexible, uniform traction) -------------------
    std::vector<int> footing;
    for (int n : mesh.top_nodes) {
        const double x = mesh.x[n];
        if (x >= x0 - 1e-9 && x <= x1 + 1e-9) footing.push_back(n);
    }
    Eigen::VectorXd f_unit = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, footing, 0.0, -1.0, f_unit);

    // --- Incremental limit analysis ------------------------------------------
    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricIndefinite);
        s->factorize(k);
        return s->solve(r);
    };

    const double nq = std::exp(kPi * std::tan(phi)) *
                      std::pow(std::tan(0.25 * kPi + 0.5 * phi), 2.0);
    const double nc_exact = (nq - 1.0) / std::tan(phi);  // (Nq-1) cot(phi)

    const double q_max = 1.1 * c * 14.835;  // ramp just past collapse
    NewtonOptions opt;
    // Adaptive load stepping refines automatically near collapse, so a modest
    // initial step count suffices (and is faster -- less past-collapse grinding).
    opt.load_steps = 30;
    opt.max_iterations = 150;
    opt.tolerance = 1e-5;
    const Eigen::VectorXd f_ext = q_max * f_unit;

    const auto result =
        katai::core::solve_nonlinear(mesh, dofs, materials, f_ext, solve, opt);

    const double q_ult = result.load_factor * q_max;
    const double nc = q_ult / c;
    const double rel_err = std::fabs(nc - nc_exact) / nc_exact;
    std::printf("  Bearing phi=20deg: Nc(FEM) = %.3f, exact = %.3f, hata = %.1f%%\n",
                nc, nc_exact, 100.0 * rel_err);

    check(rel_err < 0.05, "phi>0 bearing capacity Nc within 5% of (Nq-1)cot(phi)");
}

} // namespace

int main() {
    test_bearing_phi();
    if (g_failures == 0) {
        std::printf("OK: phi>0 bearing-capacity benchmark passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
