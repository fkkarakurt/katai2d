// Biot coupled consolidation -- Terzaghi 1D verification (Faz A.2). A laterally confined column
// with an initial uniform excess pore pressure u0, drained at the top, impermeable elsewhere,
// consolidates: the pore pressure dissipates and the column settles. The FE degree of
// consolidation U(t) = settlement(t)/settlement_inf must follow the classical Terzaghi series
//   U = 1 - sum (2/M^2) exp(-M^2 Tv),  M=(2m+1)pi/2,  Tv = cv t / H_dr^2,  cv = k Eoed/gamma_w.
// Final settlement = u0 H / Eoed. Closed-form analytic reference (gold standard).
// (See docs/references/consolidation-formulation.md; PLAXIS 2D Scientific Manual ch.4.)
//
// verify: KV-CON-001
//   oracle:   closed_form
//   source:   Terzaghi one-dimensional consolidation theory; derivation recorded in docs/references/consolidation-formulation.md (cites PLAXIS 2D Scientific Manual ch. 4)
//   locator:  U(Tv) = 1 - sum (2/M^2) exp(-M^2 Tv), M = (2m+1) pi/2; final settlement u0 H / Eoed
//   quantity: degree of consolidation U(Tv) at sampled time factors, and the final settlement [-]
//   expected: the Terzaghi series above at Tv = 0.2, 0.4, 0.6, 0.9
//   band:     0.03 absolute on U, as asserted below -- first-order time stepping with Tv steps of 0.02
#include <katai/analysis/consolidation.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Permeability;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// Terzaghi average degree of consolidation U(Tv) via the Fourier series.
double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int m = 0; m < 80; ++m) {
        const double M = (2 * m + 1) * kPi / 2.0;
        s += (2.0 / (M * M)) * std::exp(-M * M * Tv);
    }
    return 1.0 - s;
}

void test_terzaghi_1d() {
    constexpr double W = 0.2, Hc = 1.0;          // column: H_dr = Hc (single drainage, top)
    constexpr double E = 1000.0, nu = 0.0;        // -> Eoed = E = 1000 (nu=0)
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    constexpr double u0 = 10.0;                    // initial uniform excess pore pressure
    constexpr double k = 1.0e-4, gamma_w = 10.0;   // permeability, unit weight of water
    const double kw_over_n = 1.0e9;                // ~incompressible pore fluid -> cv = k Eoed/gamma_w
    const double cv = (k / gamma_w) / (1.0 / Eoed + 1.0 / kw_over_n);
    const double s_inf = u0 * Hc / Eoed;           // final settlement

    RectangularDomain domain{0.0, 0.0, W, Hc, 0};
    Mesh mesh = katai::mesh::generate_structured_tri15(domain, 1, 16);  // tri15: high-order pore field

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);   // laterally confined
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // Drainage: top nodes p=0; elsewhere impermeable (natural). Initial pore = u0 everywhere.
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;
    std::vector<double> p0(mesh.node_count, u0);

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    const std::vector<Permeability> perm = {{k, k}};

    // Time steps: Tv = cv t / Hc^2 = 0.01 t. dt=2 -> Tv step 0.02; 50 steps -> Tv up to 1.0.
    const double dt = 2.0;
    const int nsteps = 50;
    const auto r = katai::core::solve_consolidation(mesh, dofs, mm, perm, gamma_w, kw_over_n,
                                                    drained, p0, dt, nsteps);

    // Top settlement = |uy| at a top node (compression -> uy<0 in tension-positive).
    const int topn = mesh.top_nodes[mesh.top_nodes.size() / 2];
    auto settlement = [&](int step) {
        return std::fabs(r.displacement[step][dofs.global_dof(topn, 1)]);
    };
    std::printf("  Terzaghi 1D: cv=%.4f  s_inf(analytic)=%.5e  Eoed=%.0f\n", cv, s_inf, Eoed);
    std::printf("   Tv     U_FE     U_Terzaghi   (settlement/s_inf)\n");
    int checked = 0;
    for (size_t i = 1; i < r.times.size(); ++i) {
        const double Tv = cv * r.times[i] / (Hc * Hc);
        const double Ufe = settlement(i) / s_inf;
        const double Uth = terzaghi_U(Tv);
        // Check representative Tv across the curve (avoid the steep early-time boundary layer Tv<0.15).
        if ((std::fabs(Tv - 0.2) < 0.011 || std::fabs(Tv - 0.4) < 0.011 ||
             std::fabs(Tv - 0.6) < 0.011 || std::fabs(Tv - 0.9) < 0.011)) {
            std::printf("   %.3f  %.4f   %.4f       err=%.1f%%\n", Tv, Ufe, Uth,
                        100.0 * (Ufe - Uth) / Uth);
            check(std::fabs(Ufe - Uth) < 0.03, "U_FE matches Terzaghi U(Tv) within 3%");
            ++checked;
        }
    }
    check(checked >= 3, "checked several Tv points");

    // Consistency at the final step: settlement / s_inf should equal Terzaghi U at that Tv.
    const double Tvf = cv * r.times.back() / (Hc * Hc);
    const double Uf = settlement(r.times.size() - 1) / s_inf;
    std::printf("  final: Tv=%.2f  U_FE=%.4f  U_Terzaghi=%.4f  (settlement %.5e of s_inf %.5e)\n",
                Tvf, Uf, terzaghi_U(Tvf), settlement(r.times.size() - 1), s_inf);
    check(std::fabs(Uf - terzaghi_U(Tvf)) < 0.03, "final settlement consistent with Terzaghi");
}

// Vermeer & Verruijt (1981) critical time step -- calibrated against the LBB checkerboard study
// (tests/study_lbb_undrained.cpp; docs/validation/lbb-undrained-checkerboard.md). That study's 12x8 m
// block meshed 24x16 has h~0.5 m elements, E=1e4, nu=0.3, ky=1e-3, gamma_w=10, Kw/n=6.67e6 (n=0.3,
// Kw=2e6). The measured checkerboard is severe at dt=1e-8 and gone by dt~1, so dt_crit must land
// between them; the closed formula gives ~4.65e-3 day for tri6 (eta=40) and half that for tri15.
void test_critical_dt() {
    constexpr double h = 0.5, E = 1.0e4, nu = 0.3, ky = 1.0e-3, gamma_w = 10.0, n = 0.3, Kw = 2.0e6;
    const double Eoed = katai::core::oedometer_modulus(E, nu);
    check(std::fabs(Eoed - 13461.538) < 1e-2, "oedometer_modulus = E(1-nu)/((1+nu)(1-2nu))");

    const double dt6 = katai::core::consolidation_critical_dt(h, 40.0, Eoed, ky, n, Kw, gamma_w);
    const double dt15 = katai::core::consolidation_critical_dt(h, 80.0, Eoed, ky, n, Kw, gamma_w);
    std::printf("  dt_crit: tri6=%.4e day  tri15=%.4e day  (Eoed=%.1f)\n", dt6, dt15, Eoed);
    check(std::fabs(dt6 - 4.6522e-3) < 1e-5, "dt_crit(tri6) matches Vermeer-Verruijt closed form");
    check(std::fabs(dt15 - 0.5 * dt6) < 1e-12, "tri15 (eta=80) threshold is half of tri6 (eta=40)");
    // It must bracket the study's measured regime: checkerboard severe at dt=1e-8, clean by dt=1.
    check(1e-8 < dt6 && dt6 < 1.0, "dt_crit brackets the checkerboard onset (1e-8 < dt_crit < 1)");
    // Degenerate inputs -> 0 ("no constraint"), so the GUI warning treats them as adequate.
    check(katai::core::consolidation_critical_dt(h, 40.0, Eoed, 0.0, n, Kw, gamma_w) == 0.0,
          "zero permeability -> no dt_crit constraint");
    check(katai::core::oedometer_modulus(E, 0.5) == 0.0, "degenerate nu=0.5 -> oedometer modulus 0");
}

} // namespace

int main() {
    test_terzaghi_1d();
    test_critical_dt();
    if (g_failures == 0) {
        std::printf("OK: Biot consolidation verified (Terzaghi 1D U-Tv + final settlement + dt_crit)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
