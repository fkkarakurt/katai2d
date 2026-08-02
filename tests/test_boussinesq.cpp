// Boussinesq/Flamant strip-load verification (P0.8b).
//
// A uniform strip load of width 2a and intensity q on an elastic half-plane. The vertical
// stress at depth z under the center (analytic, independent of E and v):
//     σz = (q/π)(α + sin α),   α = 2·arctan(a/z)   [the angle subtended by the strip]
// The half-plane (plane strain) FEM approximates this on a large finite domain; a few %
// deviation is expected from the finite boundary + discretization (< 5% target).
// verify: KV-FND-006
//   oracle:   closed_form
//   source:   classical elasticity: Boussinesq point solution integrated over a uniform strip (independent of E and nu)
//   locator:  sigma_z = (q/pi)(alpha + sin alpha), alpha = 2 atan(a/z), beneath the strip centre (stated in full above)
//   quantity: vertical stress sigma_z along the centreline under a uniform strip load [kPa]
//   expected: the closed form above, evaluated per sampled depth
//   band:     5%, as asserted below -- a half-plane approximated by a large finite domain plus discretisation
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
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

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

constexpr double kWidth = 40.0;   // domain width
constexpr double kHeight = 20.0;  // domain depth
constexpr double kHalfWidth = 2.0;  // strip half-width a
constexpr double kPressure = 100.0; // q
constexpr double kCenterX = kWidth / 2.0;
constexpr double kPi = 3.14159265358979323846;

double boussinesq_sigma_z(double z) {
    const double alpha = 2.0 * std::atan(kHalfWidth / z);
    return (kPressure / kPi) * (alpha + std::sin(alpha));
}

Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k,
                          const Eigen::VectorXd& f) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(f);
}

// The node nearest to (x, y).
int nearest_node(const Mesh& mesh, double x, double y) {
    int best = 0;
    double best_d2 = 1e300;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double dx = mesh.x[n] - x, dy = mesh.y[n] - y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; best = n; }
    }
    return best;
}

void test_strip_load() {
    const RectangularDomain domain{0.0, 0.0, kWidth, kHeight, 0};
    const int nx = 40, ny = 20;  // 1x1 cells, 0.5 fine grid
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, ny);

    // BC: taban ankastre, yanlar rulo (ux=0).
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) {
        dofs.fix_node_component(n, 0);
        dofs.fix_node_component(n, 1);
    }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{10000.0, 0.3}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    const katai::math::CsrMatrix k = builder.build();

    // Pressure q (downward) on the central strip [xc-a, xc+a] of the top surface.
    std::vector<int> strip_nodes;
    for (int n : mesh.top_nodes) {
        const double x = mesh.x[n];
        if (x >= kCenterX - kHalfWidth - 1e-9 && x <= kCenterX + kHalfWidth + 1e-9)
            strip_nodes.push_back(n);
    }
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, strip_nodes, 0.0,
                                           -kPressure, f);

    const Eigen::VectorXd full =
        katai::core::expand_to_full(dofs, solve_spd(k, f));
    const auto field = katai::core::recover_nodal_stresses(mesh, materials, full);

    // Compare the vertical stress at several depths under the center.
    const double depths[] = {2.0, 4.0, 6.0};
    double max_rel_err = 0.0;
    for (double z : depths) {
        const int n = nearest_node(mesh, kCenterX, kHeight - z);
        const double sigma_z_fem = -field.stress[n](1);  // as compression (+)
        const double sigma_z_exact = boussinesq_sigma_z(z);
        const double rel = std::fabs(sigma_z_fem - sigma_z_exact) / sigma_z_exact;
        max_rel_err = std::fmax(max_rel_err, rel);
        std::printf("  z=%.1f: FEM σz=%.3f  analytic=%.3f  rel err=%.2f%%\n",
                    z, sigma_z_fem, sigma_z_exact, 100.0 * rel);
    }
    check(max_rel_err < 0.05, "Boussinesq vertical stress < 5% deviation");
}

} // namespace

int main() {
    test_strip_load();
    if (g_failures == 0) {
        std::printf("OK: Boussinesq strip-load verification passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
