// Elastoplastic (MC/HS) Biot consolidation -- monolithic coupled Newton (PLAXIS Sci.Man 4.3).
// There is no simple closed form for elastoplastic consolidation, so the verification is anchored two
// ways: (1) LE REDUCTION -- with a linear-elastic skeleton the coupled Newton solver must reproduce
// the validated LE consolidation core (and hence the classical Terzaghi U-Tv) to round-off; (2) an MC
// case that genuinely yields must converge every step, stay admissible, and (consolidation being
// drained at t->inf) approach the directly-DRAINED elastoplastic solution (solve_nonlinear), the
// gold-standard oracle (with a tolerance acknowledging plastic path-dependence).
#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Permeability;
using katai::core::GaussState;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 80; ++j) { const double M = (2 * j + 1) * kPi / 2.0; s += (2.0 / (M * M)) * std::exp(-M * M * Tv); }
    return 1.0 - s;
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

// Confined column shared by the cases. Returns mesh + dofs (oedometer constraint, top drained).
struct Column {
    Mesh mesh; DofMap dofs; std::vector<char> drained; int topn; double H;
    Column(Mesh m, DofMap d, std::vector<char> dr, int t, double h)
        : mesh(std::move(m)), dofs(std::move(d)), drained(std::move(dr)), topn(t), H(h) {}
};
Column make_column(double W, double H, int nx, int ny) {
    RectangularDomain domain{0.0, 0.0, W, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri15(domain, nx, ny);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;
    const int topn = mesh.top_nodes[mesh.top_nodes.size() / 2];
    return Column(std::move(mesh), std::move(dofs), std::move(drained), topn, H);
}

// (1) LE reduction: the SAME 1D Terzaghi column as test_consolidation, dissipating an initial uniform
// excess pore u0, solved through the elastoplastic coupled-Newton path with a LINEAR-ELASTIC skeleton.
// Must reproduce both the LE core solver and the Terzaghi series.
void test_le_reduction() {
    std::printf("-- (1) LE reduction: elastoplastic solver with an elastic skeleton = Terzaghi --\n");
    constexpr double W = 0.2, H = 1.0, E = 1000.0, nu = 0.0, u0 = 10.0;
    constexpr double k = 1.0e-4, gamma_w = 10.0, kw_over_n = 1.0e9;
    const double Eoed = E;  // nu = 0
    const double cv = (k / gamma_w) / (1.0 / Eoed + 1.0 / kw_over_n);
    const double s_inf = u0 * H / Eoed;
    const double dt = 2.0; const int nsteps = 50;

    Column c = make_column(W, H, 1, 16);
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> p0(c.mesh.node_count, u0);

    const auto le = katai::core::solve_consolidation(c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n,
                                                     c.drained, p0, dt, nsteps);
    const auto pl = katai::core::solve_consolidation_plastic(c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n,
                                                             c.drained, {}, p0, dt, nsteps, {}, nullptr,
                                                             nonsym_factory());
    check(pl.converged, "elastoplastic(LE) consolidation converged every step");
    check(pl.series.times.size() == le.times.size(), "same number of recorded steps");

    auto settle = [&](const katai::core::ConsolidationResult& r, int step) {
        return std::fabs(r.displacement[step][c.dofs.global_dof(c.topn, 1)]);
    };
    double max_diff = 0.0, max_err = 0.0; int checked = 0;
    for (size_t i = 1; i < pl.series.times.size(); ++i) {
        max_diff = std::fmax(max_diff, std::fabs(settle(pl.series, (int)i) - settle(le, (int)i)));
        const double Tv = cv * pl.series.times[i] / (H * H);
        if (std::fabs(Tv - 0.2) < 0.011 || std::fabs(Tv - 0.4) < 0.011 ||
            std::fabs(Tv - 0.6) < 0.011 || std::fabs(Tv - 0.9) < 0.011) {
            const double Ufe = settle(pl.series, (int)i) / s_inf, Uth = terzaghi_U(Tv);
            std::printf("   Tv=%.3f  U_plastic=%.4f  U_Terzaghi=%.4f  (%.1f%%)\n",
                        Tv, Ufe, Uth, 100.0 * (Ufe - Uth) / Uth);
            max_err = std::fmax(max_err, std::fabs(Ufe - Uth)); ++checked;
        }
    }
    std::printf("   max |settlement_plastic - settlement_LE| = %.3e (s_inf=%.3e)\n", max_diff, s_inf);
    check(checked >= 3, "checked several Tv points");
    check(max_diff < 1e-6 * s_inf, "elastoplastic(LE) == LE core solver (reduces exactly)");
    check(max_err < 0.03, "U(Tv) matches Terzaghi within 3%");
}

// (2) MC that stays ELASTIC (huge cohesion): must match the LE reduction exactly -> confirms the MC
// integrate_point branch is wired correctly in the coupled assembly.
void test_mc_elastic_equals_le() {
    std::printf("-- (2) MC (no yield) == LE in the coupled solver --\n");
    constexpr double W = 0.2, H = 1.0, E = 1000.0, nu = 0.0, u0 = 10.0;
    constexpr double k = 1.0e-4, gamma_w = 10.0, kw_over_n = 1.0e9;
    const double dt = 2.0; const int nsteps = 30;
    Column c = make_column(W, H, 1, 16);
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> p0(c.mesh.node_count, u0);

    MaterialModel le; le.type = MaterialType::LinearElastic; le.youngs_modulus = E; le.poisson_ratio = nu;
    MaterialModel mc; mc.type = MaterialType::MohrCoulomb; mc.youngs_modulus = E; mc.poisson_ratio = nu;
    mc.cohesion = 1.0e6; mc.friction_angle = 0.6; mc.dilatancy_angle = 0.0;   // c huge -> never yields
    const auto rle = katai::core::solve_consolidation_plastic(c.mesh, c.dofs, {le}, perm, gamma_w,
        kw_over_n, c.drained, {}, p0, dt, nsteps, {}, nullptr, nonsym_factory());
    const auto rmc = katai::core::solve_consolidation_plastic(c.mesh, c.dofs, {mc}, perm, gamma_w,
        kw_over_n, c.drained, {}, p0, dt, nsteps, {}, nullptr, nonsym_factory());
    check(rmc.converged, "MC(no yield) consolidation converged");
    double md = 0.0;
    for (size_t i = 0; i < rmc.series.times.size(); ++i)
        md = std::fmax(md, std::fabs(rmc.series.displacement[i][c.dofs.global_dof(c.topn, 1)] -
                                     rle.series.displacement[i][c.dofs.global_dof(c.topn, 1)]));
    std::printf("   max |u_MC - u_LE| = %.3e\n", md);
    check(md < 1e-8, "MC without yield reduces to LE in the coupled solver");
}

// (3) YIELDING drained-limit oracle: an MC footing strip, loaded UNDRAINED then CONSOLIDATED to
// t->inf (p->0). At full dissipation the effective state must approach the directly-DRAINED
// elastoplastic solution (solve_nonlinear) -- the gold-standard oracle (both share integrate_point).
// A modest plastic path-dependence (undrained-shear vs drained path) is tolerated.
void test_mc_footing_drained_limit() {
    std::printf("-- (3) MC footing: undrained -> consolidate ~ directly drained (yielding) --\n");
    constexpr double W = 6.0, Hd = 4.0, x0 = 2.4, x1 = 3.6;   // footing B = 1.2, aligned to columns
    constexpr double E = 5000.0, nu = 0.3, c = 10.0, phi = 20.0 * kPi / 180.0, q = 45.0;
    constexpr double k = 1.0e-2, gamma_w = 10.0;
    const double kw_over_n = 2.0e6 / 0.3;                     // near-incompressible water
    RectangularDomain domain{0.0, 0.0, W, Hd, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 10, 8);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // Footing chain: top nodes in [x0,x1], ordered by x (corner,mid,corner,... edge triples).
    std::vector<std::pair<double, int>> tn;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] > x0 - 1e-6 && mesh.x[n] < x1 + 1e-6) tn.push_back({mesh.x[n], n});
    std::sort(tn.begin(), tn.end());
    std::vector<int> foot; for (auto& t : tn) foot.push_back(t.second);
    std::vector<int> footing_top = foot;   // for settlement readout

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, foot, 0.0, -q, f);

    MaterialModel mc; mc.type = MaterialType::MohrCoulomb; mc.youngs_modulus = E; mc.poisson_ratio = nu;
    mc.cohesion = c; mc.friction_angle = phi; mc.dilatancy_angle = 0.0;
    const std::vector<MaterialModel> mm = {mc};
    const std::vector<Permeability> perm = {{k, k}};

    auto solve_ns = [](const katai::math::CsrMatrix& K, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric); s->factorize(K); return s->solve(r);
    };
    // Drained reference (directly drained elastoplastic; ramped load steps for plasticity).
    const auto dr = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_ns, {20, 60, 1e-7});
    check(dr.converged, "drained MC footing converged (load below drained capacity)");
    auto foot_settle = [&](const Eigen::VectorXd& d) {
        double s = 0.0; for (int n : footing_top) s = std::fmax(s, std::fabs(d[dofs.global_dof(n, 1)])); return s;
    };
    const double s_drained = foot_settle(dr.displacement);

    // Consolidation: load applied undrained at t=0, dissipate to large Tv.
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double cv = k * Eoed / gamma_w;
    const double dt = 2.0 * Hd * Hd / cv / 24.0; const int nsteps = 30;   // reach Tv ~ 2.5
    const auto co = katai::core::solve_consolidation_plastic(mesh, dofs, mm, perm, gamma_w, kw_over_n,
        drained, {}, {}, dt, nsteps, {}, &f, nonsym_factory());
    check(co.converged, "elastoplastic footing consolidation converged every step");
    if (!co.converged) return;
    const double s_consol = foot_settle(co.series.displacement.back());
    double pore_max = 0.0; for (double pv : co.series.pore.back()) pore_max = std::fmax(pore_max, std::fabs(pv));
    const double Tvf = cv * co.series.times.back() / (Hd * Hd);

    std::printf("   footing settlement: drained=%.5e  consol(t->inf)=%.5e  (%.1f%%)  Tv_f=%.2f  pore_f=%.3f\n",
                s_drained, s_consol, 100.0 * (s_consol - s_drained) / s_drained, Tvf, pore_max);
    check(pore_max < 0.05 * q, "excess pore dissipated to ~0 (fully consolidated)");
    check(s_consol > 0.0, "footing settles");
    check(std::fabs(s_consol - s_drained) < 0.05 * s_drained,
          "consolidated (t->inf) settlement = directly drained elastoplastic (within 5%)");
}

}  // namespace

int main() {
    std::printf("Elastoplastic (Biot) consolidation -- coupled Newton\n\n");
    test_le_reduction();
    test_mc_elastic_equals_le();
    test_mc_footing_drained_limit();
    if (g_failures == 0) {
        std::printf("\nOK: elastoplastic consolidation reduces to LE/Terzaghi (coupled Newton verified)\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
