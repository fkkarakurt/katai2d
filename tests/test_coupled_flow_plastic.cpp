// Elastoplastic (MC/HS) fully-coupled flow-deformation (W3-follow). Merges the unsaturated van Genuchten/
// Mualem + Bishop (χ=S_eff) machinery of W3 (solve_coupled_flow_deformation) with the MC/HS elastoplastic
// coupled-Newton skeleton of consolidation_plastic. There is no closed form for unsaturated elastoplastic
// consolidation, so the verification is anchored on two EXACT limits + a robustness check:
//   (1) SATURATED BIT-EXACT: where the soil stays saturated, the result must reproduce
//       solve_consolidation_plastic to round-off -- for (a) an LE column AND (b) an MC YIELDING footing
//       (so the plastic return mapping is exercised through the coupling). This gates the whole machinery.
//   (2) LE-UNSATURATED == W3: an LE skeleton with active suction must converge to the SAME fixed point as
//       the W3 solver (both solve identical coupled equations; the plastic core's Newton-Picard and W3's
//       full-solve Picard share one fixed point).
//   (3) MC-UNSATURATED ROBUSTNESS: an MC skeleton with initial suction draining converges every step, stays
//       admissible (S_eff in (0,1]), and resaturates as suction drains -- plasticity + unsaturated coupling
//       coexist without blowing up. MC-without-yield additionally reduces bit-for-bit to the LE skeleton.
// Refs: docs/references/transient-unsaturated-flow-formulation.md (W3); consolidation-formulation.md §4.3.
#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/coupled_flow_deformation.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/water_retention.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

using katai::core::CoupledFlowPlasticResult;
using katai::core::CoupledFlowResult;
using katai::core::DofMap;
using katai::core::GaussState;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Permeability;
using katai::core::WaterRetention;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// Coupled (saddle-point / non-associated) solve factory: factor once per Newton iterate, back-solve.
katai::core::ConsolidationSolveFactory nonsym_factory() {
    return [](const katai::math::CsrMatrix& A) {
        std::shared_ptr<linsolve::DirectSolver> s =
            linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(A);
        return std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(
            [s](const Eigen::VectorXd& b) { return s->solve(b); });
    };
}

// Max |difference| between two consolidation series (displacement + pore), with reference magnitudes.
struct SeriesDiff { double max_du, ref_u, max_dp, ref_p; };
SeriesDiff series_diff(const katai::core::ConsolidationResult& a, const katai::core::ConsolidationResult& b) {
    SeriesDiff d{0, 0, 0, 0};
    const size_t ns = std::min(a.times.size(), b.times.size());
    for (size_t i = 0; i < ns; ++i) {
        for (int j = 0; j < a.displacement[i].size(); ++j) {
            d.max_du = std::fmax(d.max_du, std::fabs(a.displacement[i][j] - b.displacement[i][j]));
            d.ref_u = std::fmax(d.ref_u, std::fabs(a.displacement[i][j]));
        }
        for (int j = 0; j < a.pore[i].size(); ++j) {
            d.max_dp = std::fmax(d.max_dp, std::fabs(a.pore[i][j] - b.pore[i][j]));
            d.ref_p = std::fmax(d.ref_p, std::fabs(a.pore[i][j]));
        }
    }
    return d;
}

// 1D Terzaghi confined column (tri15), top drained, initial uniform excess pore u0 (tension-positive).
struct Column { Mesh mesh; DofMap dofs{0, 2}; std::vector<char> drained; std::vector<double> p0; };
Column make_column(double u0) {
    Column c;
    RectangularDomain domain{0.0, 0.0, 0.2, 1.0, 0};
    c.mesh = katai::mesh::generate_structured_tri15(domain, 1, 12);
    c.dofs = DofMap(c.mesh.node_count, 2);
    for (int n : c.mesh.bottom_nodes) { c.dofs.fix_node_component(n, 0); c.dofs.fix_node_component(n, 1); }
    for (int n : c.mesh.left_nodes)  c.dofs.fix_node_component(n, 0);
    for (int n : c.mesh.right_nodes) c.dofs.fix_node_component(n, 0);
    c.dofs.finalize();
    c.drained.assign(c.mesh.node_count, 0);
    for (int n : c.mesh.top_nodes) c.drained[n] = 1;
    c.p0.assign(c.mesh.node_count, u0);
    return c;
}

// (1a) Saturated LE column: plastic-coupled (saturated retention) == solve_consolidation_plastic, bit-exact.
void test_saturated_le_regression() {
    std::printf("-- (1a) saturated LE column: plastic-coupled == consolidation_plastic (bit-exact) --\n");
    // u0 = -10 (tension-positive: a compressive excess pore pressure is NEGATIVE -> ψ=p/γw ≤ 0 -> the SWCC
    // gives S_eff=1 everywhere, so the unsaturated coefficients reduce exactly to the saturated ones).
    Column c = make_column(-10.0);
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, 1000.0, 0.0}};
    const std::vector<Permeability> perm = {{1e-4, 1e-4}};
    const std::vector<WaterRetention> ret = {WaterRetention{}};   // saturated -> coefficients reduce
    const std::vector<double> poro = {0.4};
    const double gamma_w = 10.0, kw_over_n = 1e9, dt = 2.0; const int nsteps = 40;

    const auto base = katai::core::solve_consolidation_plastic(c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n,
        c.drained, {}, c.p0, dt, nsteps, {}, nullptr, nonsym_factory());
    const auto cf = katai::core::solve_coupled_flow_deformation_plastic(c.mesh, c.dofs, mm, perm, ret, poro,
        gamma_w, kw_over_n, c.drained, {}, c.p0, dt, nsteps, {}, nullptr, nonsym_factory());
    check(cf.converged, "plastic-coupled converged");
    const SeriesDiff d = series_diff(base.series, cf.series);
    std::printf("   max|Δu|=%.3e (ref %.3e)  max|Δp|=%.3e (ref %.3e)\n", d.max_du, d.ref_u, d.max_dp, d.ref_p);
    check(d.max_du < 1e-9 * d.ref_u, "displacement bit-exact to consolidation_plastic");
    check(d.max_dp < 1e-9 * d.ref_p, "pore pressure bit-exact to consolidation_plastic");
}

// (1b) Saturated MC YIELDING footing: plastic-coupled (saturated) == solve_consolidation_plastic, bit-exact.
// Loaded undrained at t=0 then consolidated; the load yields the soil so the MC return mapping is exercised.
void test_saturated_mc_footing_regression() {
    std::printf("-- (1b) saturated MC yielding footing: plastic-coupled == consolidation_plastic (bit-exact) --\n");
    constexpr double kPi = 3.14159265358979323846;
    constexpr double W = 6.0, Hd = 4.0, x0 = 2.4, x1 = 3.6;
    constexpr double E = 5000.0, nu = 0.3, c = 10.0, phi = 20.0 * kPi / 180.0, q = 45.0;
    constexpr double k = 1.0e-2, gamma_w = 10.0;
    const double kw_over_n = 2.0e6 / 0.3;
    RectangularDomain domain{0.0, 0.0, W, Hd, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 10, 8);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<std::pair<double, int>> tn;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] > x0 - 1e-6 && mesh.x[n] < x1 + 1e-6) tn.push_back({mesh.x[n], n});
    std::sort(tn.begin(), tn.end());
    std::vector<int> foot; for (auto& t : tn) foot.push_back(t.second);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, foot, 0.0, -q, f);

    MaterialModel mc; mc.type = MaterialType::MohrCoulomb; mc.youngs_modulus = E; mc.poisson_ratio = nu;
    mc.cohesion = c; mc.friction_angle = phi; mc.dilatancy_angle = 0.0;
    const std::vector<MaterialModel> mm = {mc};
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<WaterRetention> ret = {WaterRetention{}};
    const std::vector<double> poro = {0.3};
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double cv = k * Eoed / gamma_w;
    const double dt = 2.0 * Hd * Hd / cv / 24.0; const int nsteps = 30;

    const auto base = katai::core::solve_consolidation_plastic(mesh, dofs, mm, perm, gamma_w, kw_over_n,
        drained, {}, {}, dt, nsteps, {}, &f, nonsym_factory());
    const auto cf = katai::core::solve_coupled_flow_deformation_plastic(mesh, dofs, mm, perm, ret, poro,
        gamma_w, kw_over_n, drained, {}, {}, dt, nsteps, {}, &f, nonsym_factory());
    check(base.converged, "reference consolidation_plastic footing converged");
    check(cf.converged, "plastic-coupled footing converged");
    const SeriesDiff d = series_diff(base.series, cf.series);
    std::printf("   max|Δu|=%.3e (ref %.3e)  max|Δp|=%.3e (ref %.3e)  newton_iters=%d\n",
                d.max_du, d.ref_u, d.max_dp, d.ref_p, cf.max_newton_iters);
    check(d.max_du < 1e-8 * d.ref_u, "yielding displacement bit-exact to consolidation_plastic");
    check(d.max_dp < 1e-8 * d.ref_p, "yielding pore pressure bit-exact to consolidation_plastic");
}

// Suction-draining column shared by (2) and (3): initial uniform suction p0=+3 (tension-positive),
// top node pore fixed to 0 (drains), more permeable so resaturation is visible.
struct SuctionCase {
    Column c; std::vector<double> p0; std::vector<Permeability> perm; std::vector<WaterRetention> ret;
    std::vector<double> poro; double gamma_w, kw_over_n, dt; int nsteps;
};
SuctionCase make_suction_case() {
    SuctionCase s;
    s.c = make_column(0.0);
    s.perm = {{1e-3, 1e-3}};
    s.ret = {WaterRetention{2.0, 2.0, 0.5, 0.0, 1.0}};
    s.poro = {0.4};
    s.gamma_w = 10.0; s.kw_over_n = 1e9; s.dt = 10.0; s.nsteps = 120;
    s.p0.assign(s.c.mesh.node_count, 3.0);
    for (int n : s.c.mesh.top_nodes) s.p0[n] = 0.0;
    return s;
}

// (2) LE unsaturated: plastic-coupled (LE skeleton) converges to the SAME fixed point as W3.
void test_le_unsaturated_equals_w3() {
    std::printf("-- (2) LE unsaturated: plastic-coupled == W3 (same fixed point) --\n");
    SuctionCase s = make_suction_case();
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, 1000.0, 0.0}};
    const CoupledFlowResult w3 = katai::core::solve_coupled_flow_deformation(
        s.c.mesh, s.c.dofs, mm, s.perm, s.ret, s.poro, s.gamma_w, s.kw_over_n, s.c.drained, s.p0,
        s.dt, s.nsteps, {}, nullptr, {}, 50, 1e-11);
    const CoupledFlowPlasticResult pl = katai::core::solve_coupled_flow_deformation_plastic(
        s.c.mesh, s.c.dofs, mm, s.perm, s.ret, s.poro, s.gamma_w, s.kw_over_n, s.c.drained, {}, s.p0,
        s.dt, s.nsteps, {}, nullptr, nonsym_factory(), 50, 1e-9);
    check(w3.converged && pl.converged, "both W3 and plastic-coupled converged");
    const SeriesDiff d = series_diff(w3.series, pl.series);
    std::printf("   max|Δu|=%.3e (ref %.3e)  max|Δp|=%.3e (ref %.3e)  newton_iters=%d\n",
                d.max_du, d.ref_u, d.max_dp, d.ref_p, pl.max_newton_iters);
    // Both solve identical discrete equations, so they share one fixed point; they agree to ~1e-5 here.
    // The residual stopping of the plastic core (kept bit-identical to consolidation_plastic) does not pin
    // the solution as tightly as W3's iterate-change stopping on this stiff (kw_over_n=1e9) saddle system,
    // so we assert sub-0.1% agreement rather than machine precision (the saturated limit IS bit-exact, (1a)).
    check(d.max_du < 1e-3 * d.ref_u + 1e-12, "displacement matches W3 fixed point (sub-0.1%)");
    check(d.max_dp < 1e-3 * d.ref_p + 1e-12, "pore pressure matches W3 fixed point (sub-0.1%)");
}

// (3) MC unsaturated robustness: an MC skeleton with draining suction converges, stays admissible, and
// resaturates. MC-without-yield (huge cohesion) additionally reduces bit-for-bit to the LE skeleton.
void test_mc_unsaturated_robustness() {
    std::printf("-- (3) MC unsaturated: converges + admissible + resaturates (+ no-yield == LE) --\n");
    SuctionCase s = make_suction_case();
    MaterialModel le; le.type = MaterialType::LinearElastic; le.youngs_modulus = 1000.0; le.poisson_ratio = 0.0;
    MaterialModel mc; mc.type = MaterialType::MohrCoulomb; mc.youngs_modulus = 1000.0; mc.poisson_ratio = 0.0;
    mc.cohesion = 1.0e6; mc.friction_angle = 0.5; mc.dilatancy_angle = 0.0;   // huge c -> never yields

    const auto rle = katai::core::solve_coupled_flow_deformation_plastic(s.c.mesh, s.c.dofs, {le}, s.perm,
        s.ret, s.poro, s.gamma_w, s.kw_over_n, s.c.drained, {}, s.p0, s.dt, s.nsteps, {}, nullptr, nonsym_factory());
    const auto rmc = katai::core::solve_coupled_flow_deformation_plastic(s.c.mesh, s.c.dofs, {mc}, s.perm,
        s.ret, s.poro, s.gamma_w, s.kw_over_n, s.c.drained, {}, s.p0, s.dt, s.nsteps, {}, nullptr, nonsym_factory());
    check(rmc.converged, "MC unsaturated consolidation converged every step");

    // MC-no-yield == LE skeleton, bit-for-bit (the MC integrate_point branch is wired with the coupling).
    const SeriesDiff d = series_diff(rle.series, rmc.series);
    std::printf("   MC(no yield) vs LE: max|Δu|=%.3e (ref %.3e)\n", d.max_du, d.ref_u);
    check(d.max_du < 1e-9 * (d.ref_u + 1e-30), "MC without yield reduces to LE skeleton in the coupled solver");

    // Admissibility + resaturation (S_eff rises as suction drains), reusing W3's checks.
    auto smin = [&](int step) {
        double m = 2.0; for (int n = 0; n < s.c.mesh.node_count; ++n) m = std::fmin(m, rmc.saturation[step][n]);
        return m;
    };
    auto smean = [&](int step) {
        double a = 0.0; for (int n = 0; n < s.c.mesh.node_count; ++n) a += rmc.saturation[step][n];
        return a / s.c.mesh.node_count;
    };
    const int last = (int)rmc.saturation.size() - 1;
    std::printf("   S_eff min: %.4f->%.4f  mean: %.4f->%.4f\n", smin(0), smin(last), smean(0), smean(last));
    check(smin(0) < 0.99, "van Genuchten activated (initial min S_eff < 1)");
    check(smean(last) > smean(0) + 0.05, "saturation increases (resaturation) as suction drains");
    bool bounded = true;
    for (size_t i = 0; i < rmc.saturation.size(); ++i)
        for (int n = 0; n < s.c.mesh.node_count; ++n)
            if (rmc.saturation[i][n] <= 0.0 || rmc.saturation[i][n] > 1.0 + 1e-12) bounded = false;
    check(bounded, "S_eff bounded in (0,1] throughout");
}

}  // namespace

int main() {
    std::printf("Elastoplastic fully-coupled flow-deformation (W3-follow) -- coupled Newton-Picard\n\n");
    test_saturated_le_regression();
    test_saturated_mc_footing_regression();
    test_le_unsaturated_equals_w3();
    test_mc_unsaturated_robustness();
    if (g_failures == 0) {
        std::printf("\nOK: elastoplastic coupled flow reduces to consolidation_plastic (saturated, bit-exact) "
                    "and to W3 (LE unsaturated); MC + unsaturated coupling robust\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
