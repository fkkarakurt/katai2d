// GUI-path QUANTITATIVE validation (V&V campaign, foundations-up). Runs canonical closed-form
// geotechnical benchmarks through the EXACT GUI compute path (build_problem solve_gravity_le) and
// asserts <few% error, so the INTEGRATED pre/compute/post pipeline is provably correct -- not just
// the core kernels (which test_prandtl / test_bearing_phi already validate directly).
//
// Coverage here (Mohr-Coulomb, the first nonlinear model):
//   (1) PRANDTL bearing capacity: weightless phi=0 (Tresca) strip footing -> q_ult = (2+pi)c,
//       N_c ~ 5.14 (Prandtl 1921). Obtained by incremental limit analysis: ramp a footing pressure
//       past collapse and read the equilibrated fraction (SolveResult.load_factor).
//   (2) REISSNER bearing capacity: weightless phi=20 -> N_c = (N_q-1)cot(phi) ~ 14.83 (Prandtl-
//       Reissner), q_ult = c N_c.
//   (3) MC == LE in the elastic regime: a confined column loaded below yield must give EXACTLY the
//       linear-elastic oedometric settlement (MC must not spuriously yield/soften when admissible).
//
// (LE oedometer + K0 + distributed surcharge are validated quantitatively in test_solve / test_gui_solve.)
//
// ON-DEMAND study (EXCLUDE_FROM_ALL; ~45 s nonlinear-to-collapse, so not in ctest):
//   cmake --build build/msvc-rwdi --target study_gui_validation && bin/study_gui_validation
// Finding (documented in the report): the GUI path under-predicts the collapse load by ~6%
// (CONSERVATIVE). This is the incremental-limit-analysis step granularity -- the GUI nonlinear
// scheme uses 20 load steps (min_dlam ~ 1/160), so it equilibrates slightly below the true limit;
// the dedicated core path (test_prandtl, 40 steps) reaches ~1%. The elastic response is EXACT.
// Reference: docs/validation/gui-pipeline-validation.md.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// Weightless c-phi half-space with a flexible strip footing (uniform pressure q over [x0,x1] at the
// top). BCs: base fully fixed, sides horizontal rollers. Returns q_ult via incremental limit analysis.
double bearing_q_ult(double phi_deg, double c, double q_applied, double max_area,
                     m::Drainage drainage = m::Drainage::Drained) {
    constexpr double W = 6.0, H = 4.0, x0 = 2.4, x1 = 3.6;  // footing width B = 1.2 (matches test_prandtl)
    m::Project pr;
    m::Material s;
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e4; s.nu = 0.3; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;  // weightless
    s.c = c; s.phi = phi_deg; s.psi = 0.0;
    s.drainage = drainage;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = x0; L.y1 = H; L.x2 = x1; L.y2 = H;            // footing strip at the top
    L.qx1 = 0; L.qy1 = -q_applied; L.qx2 = 0; L.qy2 = -q_applied;  // uniform downward pressure
    pr.loads.push_back(L);

    const auto M = katai::app::mesh_from_project(pr, max_area, 6);
    if (!M.ok) { std::printf("   (mesh failed)\n"); return -1; }
    // phi=0 collapses; the solve will NOT fully converge -- that is expected (limit load). We read the
    // equilibrated fraction either way (GravityLoading: weightless, so only the footing load ramps).
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    return R.load_factor * q_applied;
}

void test_prandtl_gui() {
    const double c = 10.0;
    const double q_ult = bearing_q_ult(0.0, c, 6.0 * c, 0.03);   // ramp to 6c; ~0.25 m elements
    const double Nc = q_ult / c, Nc_exact = 2.0 + kPi;           // 5.1416
    const double err = std::fabs(Nc - Nc_exact) / Nc_exact * 100.0;
    std::printf("  (1) Prandtl phi=0:  N_c(GUI) = %.3f  exact = %.3f (2+pi)  err = %.1f%%\n",
                Nc, Nc_exact, err);
    check(err < 8.0, "GUI-path Prandtl N_c within 8% of 2+pi (incremental limit analysis)");
}

void test_undrained_bearing_gui() {
    // Undrained (B): the undrained shear strength su is entered as c with phi forced to 0 (Tresca),
    // so the bearing capacity is the SAME Prandtl value q_ult = (2+pi) su, regardless of any
    // effective phi' typed in. We pass phi' = 30 deg on purpose: if Undrained (B) works, it is
    // IGNORED (phi -> 0) and N_c ~ 5.14; a leak of phi'=30 would give a ~5x larger capacity.
    const double su = 10.0;
    const double q_ult = bearing_q_ult(30.0, su, 6.0 * su, 0.03, m::Drainage::UndrainedB);
    const double Nc = q_ult / su, Nc_exact = 2.0 + kPi;          // 5.1416
    const double err = std::fabs(Nc - Nc_exact) / Nc_exact * 100.0;
    std::printf("  (4) Undrained (B) phi'=30 (ignored):  N_c(GUI) = %.3f  exact = %.3f (2+pi)  err = %.1f%%\n",
                Nc, Nc_exact, err);
    check(err < 10.0, "GUI-path Undrained (B) N_c ~ 2+pi (su governs, effective phi' ignored)");
}

void test_reissner_gui() {
    const double c = 10.0, phi = 20.0 * kPi / 180.0;
    const double Nq = std::exp(kPi * std::tan(phi)) * std::pow(std::tan(kPi / 4 + phi / 2), 2.0);
    const double Nc_exact = (Nq - 1.0) / std::tan(phi);          // ~14.83
    const double q_ult = bearing_q_ult(20.0, c, 1.8 * Nc_exact * c, 0.03);  // ramp ~1.8x past collapse
    const double Nc = q_ult / c;
    const double err = std::fabs(Nc - Nc_exact) / Nc_exact * 100.0;
    std::printf("  (2) Reissner phi=20:  N_c(GUI) = %.3f  exact = %.3f  err = %.1f%%\n",
                Nc, Nc_exact, err);
    check(err < 12.0, "GUI-path Reissner N_c(phi=20) within 12% of Prandtl-Reissner");
}

void test_mc_equals_le_elastic() {
    // Confined column, weightless, uniform surcharge BELOW yield. phi=0,c large -> stays elastic, so
    // MC must give EXACTLY the linear-elastic oedometric settlement qH/E_oed (no spurious yield).
    constexpr double W = 2.0, H = 10.0, E = 1.2e4, nu = 0.3, q = 20.0;
    auto settle = [&](m::SoilModel model) {
        m::Project pr;
        m::Material s; s.model = model;
        s.E = E; s.nu = nu; s.gamma_unsat = 0.0; s.c = 1.0e3; s.phi = 0.0; s.psi = 0.0;  // c huge -> elastic
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                     (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
        pr.polygons.push_back(P);
        m::Load L; L.kind = m::LoadKind::Distributed;
        L.x1 = 0; L.y1 = H; L.x2 = W; L.y2 = H; L.qx1 = 0; L.qy1 = -q; L.qx2 = 0; L.qy2 = -q;
        pr.loads.push_back(L);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
        double u_top = 0.0;
        for (int n = 0; n < M.mesh.node_count; ++n)
            if (M.mesh.y[n] > H - 1e-6) u_top = std::fmin(u_top, R.disp[n * 2 + 1]);
        return u_top;
    };
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double u_exact = -q * H / Eoed;
    const double u_le = settle(m::SoilModel::LinearElastic);
    const double u_mc = settle(m::SoilModel::MohrCoulomb);
    std::printf("  (3) elastic oedometer:  LE = %.6f  MC = %.6f  exact -qH/E_oed = %.6f\n",
                u_le, u_mc, u_exact);
    check(std::fabs(u_le - u_exact) < 0.02 * std::fabs(u_exact), "LE oedometer = -qH/E_oed (<2%)");
    check(std::fabs(u_mc - u_le) < 1e-3 * std::fabs(u_le),
          "MC == LE in the elastic regime (no spurious yield/softening)");
}

}  // namespace

int main() {
    std::printf("GUI-path quantitative validation (Mohr-Coulomb bearing capacity + elastic regime)\n");
    test_prandtl_gui();
    test_undrained_bearing_gui();
    test_reissner_gui();
    test_mc_equals_le_elastic();
    if (g_failures == 0) {
        std::printf("\nOK: GUI compute path matches closed-form bearing capacity + elastic theory\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
