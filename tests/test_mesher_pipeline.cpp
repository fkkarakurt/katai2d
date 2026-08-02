// End-to-end mesher -> solver verification (P1.4e).
//
// A rectangular domain is meshed by the unstructured quality mesher, converted to
// a quadratic (tri6) FEM mesh, and analysed under self-weight with lateral
// confinement (base fixed, sides on vertical rollers, free top). The exact
// geostatic solution is linear in depth,
//     sigma_yy(y) = -gamma * (H - y),   u_x = 0,
// and lies in the tri6 finite-element space, so a correct pipeline reproduces it
// to round-off on ANY valid mesh -- a strong check that the mesher feeds the
// solver a consistent, well-formed mesh.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/geom_predicates.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
using katai::mesh::Triangulation;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

constexpr double kW = 4.0, kH = 10.0;
constexpr double kGamma = 18.0;
constexpr double kE = 10000.0, kNu = 0.3;

void test_pipeline() {
    // 1. Mesh the rectangle with the unstructured quality mesher.
    const std::vector<double> px = {0.0, kW, kW, 0.0};
    const std::vector<double> py = {0.0, 0.0, kH, kH};
    const std::vector<std::array<int, 2>> segs = {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}};
    const Triangulation tg = katai::mesh::quality_mesh(px, py, segs, 25.0, 1.5);
    const Mesh mesh = katai::mesh::tri6_from_triangulation(tg, 0);

    check(mesh.element_count > 20, "mesher produced a non-trivial mesh");
    check(!mesh.boundary_nodes.empty(), "boundary nodes identified");

    // Every element must be positively oriented (valid Jacobian for tri6).
    bool ccw = true;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        if (katai::mesh::orient2d(mesh.x[a], mesh.y[a], mesh.x[b], mesh.y[b],
                                  mesh.x[c], mesh.y[c]) <= 0.0)
            ccw = false;
    }
    check(ccw, "all tri6 elements positively oriented");

    // 2. Boundary conditions by coordinate: base fixed, sides roller (u_x = 0).
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.boundary_nodes) {
        const double x = mesh.x[n], y = mesh.y[n];
        if (std::fabs(y) < 1e-7) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        if (std::fabs(x) < 1e-7 || std::fabs(x - kW) < 1e-7) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();

    // 3. Assemble K and the gravity load, solve, recover nodal stresses.
    const std::vector<LinearElastic> materials = {{kE, kNu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {kGamma}, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u =
        katai::core::expand_to_full(dofs, solver->solve(f));
    const auto field = katai::core::recover_nodal_stresses(mesh, materials, u);

    // 4. Compare against the exact geostatic field.
    double max_sigma_err = 0.0, max_ux = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double sigma_yy = field.stress[n](1);
        const double exact = -kGamma * (kH - mesh.y[n]);
        max_sigma_err = std::fmax(max_sigma_err, std::fabs(sigma_yy - exact));
        max_ux = std::fmax(max_ux, std::fabs(u[2 * n]));
    }
    std::printf("  pipeline: elems=%d nodes=%d  max|sigma_yy err|=%.3e  max|ux|=%.3e\n",
                mesh.element_count, mesh.node_count, max_sigma_err, max_ux);
    check(max_sigma_err < 1e-6 * kGamma * kH, "sigma_yy matches -gamma*(H-y)");
    check(max_ux < 1e-9, "u_x is zero (lateral confinement)");
}

} // namespace

int main() {
    test_pipeline();
    if (g_failures == 0) {
        std::printf("OK: mesher -> solver pipeline verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
