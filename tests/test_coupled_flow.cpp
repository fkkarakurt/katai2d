// Fully-coupled flow-deformation (W3) verification. PLAXIS 2D's most general analysis; unsaturated van Genuchten/
// Mualem + Bishop χ=S_eff — the generalization of Biot consolidation.
//   (1) SATURATED-LIMIT BIT-FOR-BIT REGRESSION: on a problem that stays saturated,
//       W3 == solve_consolidation (a core already Terzaghi-verified to <2.3%) → the whole
//       coupling machinery (saddle point, time integration, L/H/S blocks) is right.
//   (2) Terzaghi U(Tv) closed form (via the W3 path) <3%.
//   (3) UNSATURATED ACTIVATION + BISHOP: initial suction (negative pore) → van Genuchten
//       S_eff<1; with drainage pore→0 → resaturation (S_eff→1 monotone); the nonlinear
//       Picard converges. The Bishop χ=S_eff coupling is engaged.
// See docs/references/transient-unsaturated-flow-formulation.md (W3); PLAXIS Sci.Man §3 Eq 3-8 + §4 (Biot).
#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/coupled_flow_deformation.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/water_retention.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::CoupledFlowResult;
using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Permeability;
using katai::core::solve_consolidation;
using katai::core::solve_coupled_flow_deformation;
using katai::core::WaterRetention;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int m = 0; m < 80; ++m) { const double M = (2 * m + 1) * kPi / 2.0; s += (2.0 / (M * M)) * std::exp(-M * M * Tv); }
    return 1.0 - s;
}

// Shared Terzaghi 1D configuration (same as test_consolidation).
struct Setup {
    Mesh mesh; DofMap dofs{0, 2};
    std::vector<MaterialModel> mm; std::vector<Permeability> perm;
    std::vector<char> drained; std::vector<double> p0;
    double E, nu, Eoed, k, gamma_w, kw_over_n, u0, Hc, dt; int nsteps;
};
Setup make_terzaghi() {
    Setup s;
    s.E = 1000.0; s.nu = 0.0; s.Eoed = s.E; s.k = 1e-4; s.gamma_w = 10.0; s.kw_over_n = 1e9;
    s.u0 = -10.0; s.Hc = 1.0; s.dt = 2.0; s.nsteps = 40;  // p TENSION-POSITIVE: compressive pore pressure<0 → saturated
    RectangularDomain domain{0.0, 0.0, 0.2, s.Hc, 0};
    s.mesh = katai::mesh::generate_structured_tri6(domain, 1, 12);
    s.dofs = DofMap(s.mesh.node_count, 2);
    for (int n : s.mesh.bottom_nodes) { s.dofs.fix_node_component(n, 0); s.dofs.fix_node_component(n, 1); }
    for (int n : s.mesh.left_nodes)  s.dofs.fix_node_component(n, 0);
    for (int n : s.mesh.right_nodes) s.dofs.fix_node_component(n, 0);
    s.dofs.finalize();
    s.drained.assign(s.mesh.node_count, 0);
    for (int n : s.mesh.top_nodes) s.drained[n] = 1;
    s.p0.assign(s.mesh.node_count, s.u0);
    s.mm = {{MaterialType::LinearElastic, s.E, s.nu}};
    s.perm = {{s.k, s.k}};
    return s;
}

// (1) Saturated-limit bit-for-bit regression: W3 (stays saturated) == solve_consolidation.
void test_saturated_regression() {
    Setup s = make_terzaghi();
    const auto rc = solve_consolidation(s.mesh, s.dofs, s.mm, s.perm, s.gamma_w, s.kw_over_n,
                                        s.drained, s.p0, s.dt, s.nsteps);
    const std::vector<WaterRetention> ret = {WaterRetention{}};  // S_sat=1,S_res=0 → saturated coefficients exact
    const std::vector<double> porosity = {0.4};
    const CoupledFlowResult rw = solve_coupled_flow_deformation(
        s.mesh, s.dofs, s.mm, s.perm, ret, porosity, s.gamma_w, s.kw_over_n, s.drained, s.p0, s.dt, s.nsteps);
    check(rw.converged, "(1) W3 converged");
    double maxd = 0.0, maxp = 0.0, refd = 0.0, refp = 0.0;
    for (size_t i = 0; i < rc.times.size(); ++i) {
        for (int j = 0; j < rc.displacement[i].size(); ++j) {
            maxd = std::fmax(maxd, std::fabs(rc.displacement[i][j] - rw.series.displacement[i][j]));
            refd = std::fmax(refd, std::fabs(rc.displacement[i][j]));
        }
        for (int j = 0; j < rc.pore[i].size(); ++j) {
            maxp = std::fmax(maxp, std::fabs(rc.pore[i][j] - rw.series.pore[i][j]));
            refp = std::fmax(refp, std::fabs(rc.pore[i][j]));
        }
    }
    std::printf("  (1) saturated regression vs consolidation: max|Δu|=%.3e (ref %.3e), max|Δp|=%.3e (ref %.3e)\n",
                maxd, refd, maxp, refp);
    check(maxd < 1e-9 * refd, "(1) displacement bit-exact to consolidation");
    check(maxp < 1e-9 * refp, "(1) pore pressure bit-exact to consolidation");
    check(rw.max_picard_iters <= 2, "(1) saturated → Picard converges in 1-2 iters");
}

// (2) Terzaghi U(Tv) via the W3 path (closed form).
void test_terzaghi_via_w3() {
    Setup s = make_terzaghi();
    const double cv = (s.k / s.gamma_w) / (1.0 / s.Eoed + 1.0 / s.kw_over_n);
    const double s_inf = std::fabs(s.u0) * s.Hc / s.Eoed;
    const std::vector<WaterRetention> ret = {WaterRetention{}};
    const std::vector<double> porosity = {0.4};
    const CoupledFlowResult rw = solve_coupled_flow_deformation(
        s.mesh, s.dofs, s.mm, s.perm, ret, porosity, s.gamma_w, s.kw_over_n, s.drained, s.p0, s.dt, s.nsteps);
    const int topn = s.mesh.top_nodes[s.mesh.top_nodes.size() / 2];
    int checked = 0;
    for (size_t i = 1; i < rw.series.times.size(); ++i) {
        const double Tv = cv * rw.series.times[i] / (s.Hc * s.Hc);
        if (std::fabs(Tv - 0.2) < 0.011 || std::fabs(Tv - 0.5) < 0.011 || std::fabs(Tv - 0.8) < 0.011) {
            const double U = std::fabs(rw.series.displacement[i][s.dofs.global_dof(topn, 1)]) / s_inf;
            std::printf("    Tv=%.3f  U_W3=%.4f  U_Terzaghi=%.4f\n", Tv, U, terzaghi_U(Tv));
            check(std::fabs(U - terzaghi_U(Tv)) < 0.03, "(2) W3 settlement matches Terzaghi U(Tv)");
            ++checked;
        }
    }
    check(checked >= 2, "(2) checked Terzaghi points");
}

// (3) Unsaturated activation + Bishop: initial suction → S_eff<1; resaturation with drainage (S_eff↑ trend).
void test_unsaturated_activation() {
    Setup s = make_terzaghi();
    s.perm = {{1e-3, 1e-3}};   // more permeable → resaturation visible in reasonable time
    const double dt = 10.0; const int nsteps = 120;
    // Initial UNIFORM SUCTION p0=+3 (tension-pos: p>0 = suction; except the drained top p=0) → ψ=p/γw=0.3 → S_eff<1.
    std::vector<double> p0(s.mesh.node_count, 3.0);
    for (int n : s.mesh.top_nodes) p0[n] = 0.0;
    const std::vector<WaterRetention> ret = {WaterRetention{2.0, 2.0, 0.5, 0.0, 1.0}};
    const std::vector<double> porosity = {0.4};
    const CoupledFlowResult rw = solve_coupled_flow_deformation(
        s.mesh, s.dofs, s.mm, s.perm, ret, porosity, s.gamma_w, s.kw_over_n, s.drained, p0, dt, nsteps);
    check(rw.converged, "(3) unsaturated W3 converged");
    std::printf("  (3) unsaturated activation: max Picard iters=%d\n", rw.max_picard_iters);

    auto stat = [&](int step, bool mean) {
        double acc = mean ? 0.0 : 2.0; int cnt = 0;
        for (int n = 0; n < s.mesh.node_count; ++n) { const double S = rw.saturation[step][n];
            if (mean) { acc += S; ++cnt; } else acc = std::fmin(acc, S); }
        return mean ? acc / cnt : acc;
    };
    const int last = (int)rw.saturation.size() - 1;
    const double Smin_i = stat(0, false), Smin_f = stat(last, false);
    const double Smean_i = stat(0, true), Smean_f = stat(last, true);
    const double S_e_hand = 1.0 / std::sqrt(1.0 + std::pow(2.0 * 0.3, 2.0));  // van Genuchten n=2, ψ=0.3
    std::printf("    S_eff min: %.4f→%.4f  mean: %.4f→%.4f  (initial bulk hand %.4f)\n",
                Smin_i, Smin_f, Smean_i, Smean_f, S_e_hand);
    check(Smin_i < 0.99, "(3) van Genuchten activated (initial min S_eff < 1)");
    check(std::fabs(Smin_i - S_e_hand) < 1e-3, "(3) initial S_eff matches van Genuchten closed form");
    check(Smean_f > Smean_i + 0.05, "(3) saturation increases (resaturation) as suction drains");
    check(Smin_f > Smin_i, "(3) min saturation rises (deepest node also resaturates)");
    // S_eff ∈ (0,1] at all nodes.
    bool bounded = true;
    for (size_t i = 0; i < rw.saturation.size(); ++i)
        for (int n = 0; n < s.mesh.node_count; ++n)
            if (rw.saturation[i][n] <= 0.0 || rw.saturation[i][n] > 1.0 + 1e-12) bounded = false;
    check(bounded, "(3) S_eff bounded in (0,1] throughout");
    check(rw.max_picard_iters >= 2, "(3) unsaturated nonlinearity exercised Picard (>1 iter)");
}

}  // namespace

int main() {
    test_saturated_regression();
    test_terzaghi_via_w3();
    test_unsaturated_activation();
    if (g_failures == 0) {
        std::printf("OK: coupled flow-deformation verified (saturated==consolidation bit-exact + "
                    "Terzaghi U(Tv) + unsaturated van Genuchten/Bishop activation)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
