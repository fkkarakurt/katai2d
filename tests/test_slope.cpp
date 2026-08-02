// Slope stability by phi-c reduction (P1.7) -- the first sellable vertical.
//
// Rocscience / Slide verification #1 (after Griffiths & Lane): a homogeneous 1:2
// slope on a foundation layer, gamma = 20.2, c = 3 kPa, phi = 19.6 deg, psi = 0.
// Reference factor of safety: Bishop 0.988, Spencer/GLE 0.987, Phase2 T6 0.997.
// This foundation geometry has all corners obtuse/right; sharp (< 60 deg) corners
// are now handled too (concentric-shell subsegment splitting in the mesher, see
// test_delaunay::test_acute_corner_domain). We mesh with the unstructured quality
// mesher, build a tri6 mesh, load under gravity, and bisect the strength reduction
// factor to collapse.
#include <katai/analysis/strength_reduction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::StrengthReductionOptions;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
using katai::mesh::Triangulation;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
constexpr double kDeg = 3.14159265358979323846 / 180.0;

void test_slope() {
    // Foundation slope (CCW): base (20,20)-(70,20), right (70,20)-(70,35), top
    // (70,35)-(50,35), slope face (50,35)-(30,25), foundation top (30,25)-(20,25),
    // back (20,25)-(20,20). All corners >= 90 deg.
    const std::vector<double> px = {20.0, 70.0, 70.0, 50.0, 30.0, 20.0};
    const std::vector<double> py = {20.0, 20.0, 35.0, 35.0, 25.0, 25.0};
    const std::vector<std::array<int, 2>> segs = {
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 4}}, {{4, 5}}, {{5, 0}}};
    const Triangulation tg = katai::mesh::quality_mesh(px, py, segs, 20.0, 5.0);
    const Mesh mesh = katai::mesh::tri6_from_triangulation(tg, 0);
    std::printf("  slope mesh: %d elements, %d nodes\n", mesh.element_count,
                mesh.node_count);

    // Base fixed, back/toe vertical sides on vertical rollers.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.boundary_nodes) {
        const double x = mesh.x[n], y = mesh.y[n];
        if (std::fabs(y - 20.0) < 1e-7) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        if (std::fabs(x - 20.0) < 1e-7 || std::fabs(x - 70.0) < 1e-7) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();

    // Gravity body force (free-DOF vector).
    Eigen::VectorXd gravity = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {20.2}, gravity);

    const MaterialModel base{MaterialType::MohrCoulomb, 1.0e5, 0.3,
                             3.0, 19.6 * kDeg, 0.0};

    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(k);
        return s->solve(r);
    };

    StrengthReductionOptions opt;
    opt.srf_min = 0.5;
    opt.srf_max = 2.5;
    opt.bisection_iterations = 9;  // resolution ~2/2^9 ~= 4e-3
    opt.newton.load_steps = 4;
    opt.newton.max_iterations = 150;
    opt.newton.tolerance = 1e-3;

    const double fos =
        katai::core::factor_of_safety(mesh, dofs, gravity, base, solve, opt);
    const double reference = 0.99;  // LEM 0.987-0.988, Phase2 T6 0.997
    const double rel_err = std::fabs(fos - reference) / reference;
    std::printf("  phi-c reduction: FoS = %.3f  (ref ~0.99, T6 0.997)  hata = %.1f%%\n",
                fos, 100.0 * rel_err);
    check(rel_err < 0.05, "slope factor of safety within 5% of the benchmark");
}

} // namespace

int main() {
    test_slope();
    if (g_failures == 0) {
        std::printf("OK: slope stability (phi-c reduction) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
