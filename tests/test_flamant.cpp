// Point / line load on an elastic half-space -- Flamant verification (Faz A.3). A point load in
// plane strain FEM is a concentrated nodal force = a line load P (force per unit out-of-plane
// length). Flamant's solution gives the vertical stress under it; on the centerline (x=0):
//     sigma_z = 2 P / (pi z)
// and the load must be transmitted in full (vertical equilibrium: sum of base reactions = P).
// Concentrated loads are singular at the load point, so we compare at depth (away from it).
// Reference: Flamant (1892); any elasticity text. (Distributed load is covered by test_boussinesq.)
//
// verify: KV-FND-007
//   oracle:   closed_form
//   source:   Flamant (1892) line load on an elastic half-plane
//   locator:  sigma_z = 2 P / (pi z) on the centreline (stated in full above)
//   quantity: vertical stress under a concentrated line load, compared at depth away from the singular point [kPa]
//   expected: the closed form above, evaluated per sampled depth
//   band:     6%, as asserted below -- the concentrated nodal load is singular at its point of application
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
constexpr double kPi = 3.14159265358979323846;
constexpr double kW = 40.0, kH = 20.0, kCx = 20.0, kP = 100.0;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& f) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(k);
    return s->solve(f);
}
int nearest_node(const Mesh& mesh, double x, double y) {
    int best = 0; double bd = 1e300;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double dx = mesh.x[n] - x, dy = mesh.y[n] - y, d2 = dx * dx + dy * dy;
        if (d2 < bd) { bd = d2; best = n; }
    }
    return best;
}

void test_flamant_line_load() {
    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 40, 20);  // fine 0.5

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{10000.0, 0.3}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    const katai::math::CsrMatrix k = builder.build();

    // Concentrated vertical nodal force P (downward) at the top-centre node = a line load.
    const int load_node = nearest_node(mesh, kCx, kH);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const int eq = dofs.equation(dofs.global_dof(load_node, 1));
    f(eq) = -kP;

    const Eigen::VectorXd full = katai::core::expand_to_full(dofs, solve_spd(k, f));
    const auto field = katai::core::recover_nodal_stresses(mesh, materials, full);

    // Centerline vertical stress vs Flamant 2P/(pi z), at depth (away from the singular point).
    const double depths[] = {4.0, 6.0, 8.0, 10.0};
    double max_rel = 0.0;
    for (double z : depths) {
        const int n = nearest_node(mesh, kCx, kH - z);
        const double sz_fem = -field.stress[n](1);
        const double sz_exact = 2.0 * kP / (kPi * z);
        const double rel = std::fabs(sz_fem - sz_exact) / sz_exact;
        max_rel = std::fmax(max_rel, rel);
        std::printf("  z=%.1f: FEM sz=%.4f  Flamant=%.4f  rel=%.2f%%\n", z, sz_fem, sz_exact, 100.0 * rel);
    }
    check(max_rel < 0.06, "Flamant centerline sigma_z = 2P/(pi z) within 6% (concentrated load)");
}

} // namespace

int main() {
    test_flamant_line_load();
    if (g_failures == 0) {
        std::printf("OK: Flamant point/line load verified (centerline sigma_z = 2P/pi z)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
