// K0 procedure: geostatic initial stress + solver prestress seeding (P2.2 / staged
// construction foundation). Two PLAXIS initial-stress methods are exercised:
//
// (a) K0 PROCEDURE -- prescribe sigma'_v = -gamma'(z_surf - y), sigma'_h = K0 sigma'_v
//     directly at every Gauss point (K0 = Jaky 1-sin(phi'), a free value, NOT tied to
//     nu'). Seed it into the solver as the committed prestress and apply self-weight in
//     ONE load step. Because horizontal layering makes this field self-equilibrated for
//     ANY K0, the body is already in equilibrium with its weight -> the initial phase
//     produces ZERO displacement (round-off) and the committed stress stays the K0
//     field. This is exactly PLAXIS's "K0 initial phase yields no displacement" rule and
//     the prestress-carrying hook that staged construction needs.
//
// (b) GRAVITY LOADING -- no prestress; apply self-weight and solve elastically. The
//     laterally confined (oedometer) response forces sigma'_h / sigma'_v = nu'/(1-nu'),
//     the "natural" elastic K0, which differs from the free K0 of method (a).
//
// (See docs/references/initial-stress-k0.md.)
#include <katai/analysis/initial_stress.hpp>
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

Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(r);
}

constexpr double W = 4.0, H = 10.0;
constexpr double gamma = 18.0;     // dry unit weight
constexpr double E = 1.0e4, nu = 0.3;
constexpr double phi = 30.0 * 3.14159265358979323846 / 180.0;

// Laterally confined column: bottom fixed, side rollers (u_x = 0), top free.
struct Setup { Mesh mesh; DofMap dofs; };
Setup build_column() {
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 5);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    return {std::move(mesh), std::move(dofs)};
}

// Gauss-point y for tri6 (matches the solver / initial-stress ordering).
double gauss_y(const Mesh& mesh, int e, const katai::core::tri6::GaussPoint& gp) {
    const auto N = katai::core::tri6::shape_functions(gp.xi, gp.eta);
    double y = 0.0;
    for (int k = 0; k < 6; ++k) y += N(k) * mesh.y[mesh.node_of(e, k)];
    return y;
}

void test_k0_procedure() {
    Setup s = build_column();
    const double K0 = 1.0 - std::sin(phi);  // Jaky -> 0.5 for phi = 30 deg

    // (a) K0 procedure: prescribe the geostatic stress field as the prestress.
    const K0Options k0opt{H, gamma, K0};
    const std::vector<GaussState> init =
        katai::core::compute_k0_initial_stress(s.mesh, k0opt);

    // Self-weight as the external load; one load step (the prestress already balances
    // the FULL weight, so it must not be ramped).
    Eigen::VectorXd f = Eigen::VectorXd::Zero(s.dofs.equation_count());
    katai::core::assemble_gravity(s.mesh, s.dofs, {gamma}, f);

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const auto r = katai::core::solve_nonlinear(s.mesh, s.dofs, mm, f, solve_spd,
                                                {1, 50, 1e-10}, init);
    check(r.converged, "K0 initial phase converged");

    // The prestressed body is in equilibrium with its weight -> zero displacement.
    const double max_u = r.displacement.cwiseAbs().maxCoeff();
    std::printf("  K0 procedure (K0=%.3f): max|u|=%.3e  iters=%d\n",
                K0, max_u, r.total_iterations);
    check(max_u < 1e-9 * gamma * H / E * H, "K0 initial phase produces zero displacement");

    // Committed stress is unchanged = the prescribed K0 field (sigma'_h = K0 sigma'_v).
    const auto gps = katai::core::tri6::gauss_points();
    double max_sv_err = 0.0, max_ratio_err = 0.0;
    for (int e = 0; e < s.mesh.element_count; ++e) {
        for (int g = 0; g < 3; ++g) {
            const int gi = e * 3 + g;
            const GaussState& gs = r.gauss_states[gi];
            const double y = gauss_y(s.mesh, e, gps[g]);
            const double sv_exact = -gamma * (H - y);
            max_sv_err = std::fmax(max_sv_err, std::fabs(gs.stress(1) - sv_exact));
            if (std::fabs(gs.stress(1)) > 1e-6 * gamma * H)
                max_ratio_err = std::fmax(
                    max_ratio_err, std::fabs(gs.stress(0) / gs.stress(1) - K0));
        }
    }
    check(max_sv_err < 1e-6 * gamma * H, "K0 sigma'_v = -gamma(H-y)");
    check(max_ratio_err < 1e-9, "K0 sigma'_h / sigma'_v = K0 (free Jaky value)");
}

void test_gravity_loading_contrast() {
    Setup s = build_column();
    Eigen::VectorXd f = Eigen::VectorXd::Zero(s.dofs.equation_count());
    katai::core::assemble_gravity(s.mesh, s.dofs, {gamma}, f);

    // No prestress: gravity loading from zero stress -> elastic confined response.
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const auto r = katai::core::solve_nonlinear(s.mesh, s.dofs, mm, f, solve_spd,
                                                {1, 50, 1e-12});
    check(r.converged, "gravity loading converged");

    const double k0_elastic = nu / (1.0 - nu);  // 0.4286 for nu = 0.3
    double max_ratio_err = 0.0;
    for (int e = 0; e < s.mesh.element_count; ++e) {
        for (int g = 0; g < 3; ++g) {
            const int gi = e * 3 + g;
            const GaussState& gs = r.gauss_states[gi];
            if (std::fabs(gs.stress(1)) > 1e-3 * gamma * H)
                max_ratio_err = std::fmax(
                    max_ratio_err, std::fabs(gs.stress(0) / gs.stress(1) - k0_elastic));
        }
    }
    std::printf("  gravity loading: K0_elastic=nu/(1-nu)=%.4f  max ratio err=%.3e\n",
                k0_elastic, max_ratio_err);
    check(max_ratio_err < 1e-9, "gravity loading forces K0 = nu/(1-nu)");
}

} // namespace

int main() {
    test_k0_procedure();
    test_gravity_loading_contrast();
    if (g_failures == 0) {
        std::printf("OK: K0 procedure + gravity loading initial stress verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
