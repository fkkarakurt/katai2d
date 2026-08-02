// Staged construction -- excavation (P2.2). A laterally confined (oedometer) column
// under a K0 initial stress state is excavated: the top layer (y > h) is deactivated.
// Removing that overburden unloads the remaining column, which rebounds (heave).
//
// The excavation mechanics fall out automatically from (i) the carried-over initial
// stress sigma0 (solve_nonlinear initial_state), (ii) the active-element mask, and
// (iii) self-weight over active elements: the first residual = f_grav,active -
// integral_active B^T sigma0 equals exactly the release of the removed overburden's
// support (the surcharge -gamma(H-h) that the removed soil applied at y=h). The
// confined column then rebounds with a uniform vertical strain.
//
// Analytic (linear elastic, confined constrained modulus M' = D'(1,1)):
//   unloading    Delta sigma'_v = gamma (H - h)         (uniform)
//   heave at y=h u_y = Delta sigma'_v / M' * h = gamma (H-h) h / M'
//   final stress sigma'_v(y) = -gamma (h - y)           (new surface stress-free)
// Linear fields -> reproduced to round-off.
// (See docs/references/initial-stress-k0.md and staged_construction.hpp.)
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/staged_construction.hpp>
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
using katai::core::GaussState;
using katai::core::K0Options;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(r);
}

void test_excavation_heave() {
    constexpr double W = 4.0, H = 10.0, h = 6.0;  // excavate top H-h = 4 m
    constexpr double gamma = 18.0, K0 = 0.5;
    constexpr double E = 1.0e4, nu = 0.3;

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 5);  // rows align at y=6

    // Active set: elements below the excavation level (centroid y < h).
    std::vector<char> active(mesh.element_count, 1);
    for (int e = 0; e < mesh.element_count; ++e) {
        double yc = 0.0;
        for (int k = 0; k < 3; ++k) yc += mesh.y[mesh.node_of(e, k)] / 3.0;
        active[e] = yc < h ? 1 : 0;
    }

    // BCs: bottom fixed, side rollers; orphaned (excavated) nodes fully fixed.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    katai::core::fix_inactive_nodes(mesh, active, dofs);
    dofs.finalize();

    // Initial stress (K0) over the full column, carried as prestress; self-weight over
    // the ACTIVE elements only.
    const std::vector<GaussState> init =
        katai::core::compute_k0_initial_stress(mesh, K0Options{H, gamma, K0});
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {gamma}, f, active);

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_spd,
                                                {1, 50, 1e-12}, init, active);
    check(r.converged, "excavation phase converged");

    const double Mp = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double heave_exact = gamma * (H - h) * h / Mp;

    // Heave = max upward displacement (occurs at the excavation surface y = h).
    double heave = 0.0;
    for (int n = 0; n < mesh.node_count; ++n)
        heave = std::fmax(heave, r.displacement[dofs.global_dof(n, 1)]);
    std::printf("  excavation: heave(y=%.0f)=%.4e  exact=%.4e  iters=%d\n",
                h, heave, heave_exact, r.total_iterations);
    check(close(heave, heave_exact, 1e-6), "rebound heave = gamma(H-h)h/M'");

    // Final effective vertical stress in the remaining column: sigma'_v = -gamma(h-y).
    const auto gps = katai::core::tri6::gauss_points();
    double max_sv_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active[e]) continue;
        for (int g = 0; g < 3; ++g) {
            const auto N = katai::core::tri6::shape_functions(gps[g].xi, gps[g].eta);
            double y = 0.0;
            for (int k = 0; k < 6; ++k) y += N(k) * mesh.y[mesh.node_of(e, k)];
            const double sv_exact = -gamma * (h - y);  // new equilibrium
            max_sv_err = std::fmax(
                max_sv_err, std::fabs(r.gauss_states[e * 3 + g].stress(1) - sv_exact));
        }
    }
    std::printf("  remaining column: max|sigma'_v - (-gamma(h-y))|=%.3e\n", max_sv_err);
    check(max_sv_err < 1e-6 * gamma * H, "post-excavation sigma'_v = -gamma(h-y)");
}

} // namespace

int main() {
    test_excavation_heave();
    if (g_failures == 0) {
        std::printf("OK: staged construction (excavation) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
