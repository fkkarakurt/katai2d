// Staged construction -- fill (element activation) + multi-phase chaining (P2.2).
// A confined (oedometer) column is built in two lifts; each phase carries the previous
// phase's committed stress as its prestress (initial_state) and resets displacement
// (incremental per phase) -- exactly PLAXIS staged construction.
//
//   Phase 1: activate lift 1 [0, H/2], zero initial stress, self-weight.
//            -> sigma'_v = -gamma(H/2 - y),  top(y=H/2) settle = -gamma(H/2)^2/(2M')
//   Phase 2: activate BOTH lifts; lift 1 carries its phase-1 stress, lift 2 is fresh
//            (zero); self-weight over both. The residual is exactly lift 2's weight.
//
// Analytic (linear elastic, confined constrained modulus M' = D'(1,1)):
//   final stress (path-independent):  sigma'_v(y) = -gamma(H - y)      [equilibrium]
//   phase-2 incremental top settle:   u_y(H) = -3 gamma H^2 / (8 M')
// The staged settlement (3/8) is LESS than a single-pour fill (-gamma H^2/(2M') = 4/8)
// because lift 1 already settled under its own weight in phase 1 and that is not
// re-counted in phase 2. Linear fields -> round-off.
// (See docs/references/initial-stress-k0.md and staged_construction.hpp.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::GaussState;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::NewtonResult;
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

constexpr double W = 4.0, H = 8.0;     // two lifts of H/2 = 4 (rows align at y=4)
constexpr double gamma = 18.0;
constexpr double E = 1.0e4, nu = 0.3;

// One staged phase: build the active DofMap (bottom fixed, side rollers, orphaned nodes
// fixed), assemble self-weight over the active set, and solve carrying the prior stress.
NewtonResult run_phase(const Mesh& mesh, const std::vector<char>& active,
                       const std::vector<GaussState>& prestress) {
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    katai::core::fix_inactive_nodes(mesh, active, dofs);
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {gamma}, f, active);

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    NewtonResult r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_spd,
                                                  {1, 50, 1e-12}, prestress, active);
    // Attach the active DofMap result via full displacement (already total_dofs sized).
    return r;
}

// Vertical displacement at nodes on the line y = y_target (global_dof = node*2 + 1,
// structural and independent of constraints; r.displacement is full total_dofs sized).
double settle_at_y(const Mesh& mesh, const NewtonResult& r, double y_target) {
    double s = 0.0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.y[n] - y_target) < 1e-9)
            s = std::fmin(s, r.displacement[n * 2 + 1]);
    return s;
}

void test_staged_fill() {
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 4);

    auto lift_mask = [&](double y_lo, double y_hi) {
        std::vector<char> a(mesh.element_count, 0);
        for (int e = 0; e < mesh.element_count; ++e) {
            double yc = 0.0;
            for (int k = 0; k < 3; ++k) yc += mesh.y[mesh.node_of(e, k)] / 3.0;
            a[e] = (yc > y_lo && yc < y_hi) ? 1 : 0;
        }
        return a;
    };
    const std::vector<char> lift1 = lift_mask(0.0, H / 2);
    const std::vector<char> both = lift_mask(0.0, H);

    const double Mp = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Phase 1: place lift 1 (fresh, zero stress) under self-weight.
    const NewtonResult r1 = run_phase(mesh, lift1, {});
    check(r1.converged, "phase 1 (lift 1) converged");
    const double settle1 = settle_at_y(mesh, r1, H / 2);
    const double settle1_exact = -gamma * (H / 2) * (H / 2) / (2.0 * Mp);
    std::printf("  phase 1: lift-1 top settle=%.4e  exact=%.4e\n", settle1, settle1_exact);
    check(close(settle1, settle1_exact, 1e-6), "phase 1 self-weight settlement");

    // Phase 2: place lift 2; carry lift-1 stress as prestress; self-weight over both.
    const NewtonResult r2 = run_phase(mesh, both, r1.gauss_states);
    check(r2.converged, "phase 2 (lift 2) converged");

    // (1) Final stress is path-independent: sigma'_v = -gamma(H - y) everywhere.
    const auto gps = katai::core::tri6::gauss_points();
    double max_sv_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        for (int g = 0; g < 3; ++g) {
            const auto N = katai::core::tri6::shape_functions(gps[g].xi, gps[g].eta);
            double y = 0.0;
            for (int k = 0; k < 6; ++k) y += N(k) * mesh.y[mesh.node_of(e, k)];
            const double sv_exact = -gamma * (H - y);
            max_sv_err = std::fmax(
                max_sv_err, std::fabs(r2.gauss_states[e * 3 + g].stress(1) - sv_exact));
        }
    }
    std::printf("  phase 2: max|sigma'_v - (-gamma(H-y))|=%.3e\n", max_sv_err);
    check(max_sv_err < 1e-6 * gamma * H, "final sigma'_v = -gamma(H-y) (path-independent)");

    // (2) Phase-2 incremental top settlement = -3 gamma H^2 / (8 M').
    const double settle2 = settle_at_y(mesh, r2, H);
    const double settle2_exact = -3.0 * gamma * H * H / (8.0 * Mp);
    const double single_pour = -gamma * H * H / (2.0 * Mp);
    std::printf("  phase 2: top settle=%.4e  exact(staged)=%.4e  single-pour=%.4e\n",
                settle2, settle2_exact, single_pour);
    check(close(settle2, settle2_exact, 1e-6),
          "phase-2 incremental top settlement = -3 gamma H^2/(8 M') (staged < single-pour)");
}

} // namespace

int main() {
    test_staged_fill();
    if (g_failures == 0) {
        std::printf("OK: staged construction (fill + multi-phase chaining) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
