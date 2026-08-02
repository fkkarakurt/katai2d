// Consolidation through the GUI compute path (build_problem PhaseType::Consolidation + solve_phases).
// The Biot CORE is validated quantitatively in test_consolidation (Terzaghi 1D U-Tv, dense LU);
// here the INTEGRATED path -- project + materials + flow drainage BC -> mesh -> initial K0 phase ->
// consolidation phase (surcharge applied at t=0+, excess pore dissipates) -- must reproduce the
// same closed form. A laterally-confined column, drained at the top, loaded by a top surcharge q:
//   final settlement  s_inf = q H / Eoed        (drained oedometric)
//   degree of consol.  U(Tv) = 1 - sum (2/M^2) exp(-M^2 Tv),  Tv = cv t / H_dr^2,  cv = k Eoed/gamma_w.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 80; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += (2.0 / (M * M)) * std::exp(-M * M * Tv);
    }
    return 1.0 - s;
}

void test_terzaghi_gui(m::SoilModel model, const char* label) {
    std::printf("-- 1D Terzaghi consolidation through the GUI phase path (%s) --\n", label);
    constexpr double W = 1.0, H = 12.0;          // column: single (top) drainage -> H_dr = H
    constexpr double E = 1000.0, nu = 0.0;        // nu = 0 -> Eoed = E
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    constexpr double k = 0.1, q = 10.0;           // permeability [m/day], surcharge [kPa]
    const double gw = katai::app::kGammaWater;
    const double cv = k * Eoed / gw;              // near-incompressible pore fluid
    const double s_inf = q * H / Eoed;            // final (drained) settlement

    m::Project pr;
    m::Material s; s.model = model;
    s.E = E; s.nu = nu; s.gamma_unsat = 16.0; s.gamma_sat = 18.0; s.e_init = 0.5;
    s.c = 50.0; s.phi = 25.0; s.psi = 0.0;   // MC: K0 confined compression stays elastic -> = Terzaghi
    s.kx = k; s.ky = k;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0};
    P.y = {0, 0, H, H};
    //            bottom                  right                       top               left
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};   // top drains (p = 0)
    P.edge_head = {0.0, 0.0, H, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;   // pore = excess only (no hydrostatic background)

    // Surcharge on the top edge -- installed in the consolidation phase (not the initial phase).
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 0; L.y1 = H; L.x2 = W; L.y2 = H;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -q;
    pr.loads.push_back(L);

    pr.initial.load_active = {0};                 // surcharge OFF in the initial K0 phase
    // Duration chosen so Tv = cv t / H^2 reaches ~2 (near-complete consolidation, U ~ 0.99).
    m::Phase consol; consol.name = "Consolidation";
    consol.type = m::PhaseType::Consolidation;
    consol.duration = 2.0 * H * H / cv; consol.time_steps = 120;
    consol.load_active = {1};                      // surcharge ON -> applied at t = 0+
    pr.phases.push_back(consol);

    const auto M = katai::app::mesh_from_project(pr, 0.08, 6);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("  (%s)\n", M.message.c_str()); return; }

    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2, "initial + consolidation phases ran");
    if (res.size() != 2) return;
    check(res[0].ok, "initial K0 phase converged");
    check(res[1].ok, "consolidation phase solved");
    if (!res[1].ok) { std::printf("  (%s)\n", res[1].message.c_str()); return; }
    const auto& C = res[1];
    check(!C.consol_time.empty(), "consolidation produced a time series");
    if (C.consol_time.empty()) return;

    std::printf("  cv=%.4f  s_inf(analytic)=%.5e m  Eoed=%.0f  steps=%zu\n",
                cv, s_inf, Eoed, C.consol_time.size());

    // Peak excess pore (undrained generation, early time, far from drainage) ~ surcharge q; by the
    // end (Tv ~ 2) it has dissipated to ~0. The peak is index-robust (the load is applied in the
    // first time step, so the t=0 record is the pre-load zero state).
    double pore_peak = 0.0;
    for (double p : C.consol_excess_pore) pore_peak = std::fmax(pore_peak, p);
    std::printf("   excess pore: peak = %.3f kPa (~q=%.1f),  final = %.4f kPa\n",
                pore_peak, q, C.consol_excess_pore.back());
    check(pore_peak > 0.85 * q && pore_peak < 1.10 * q,
          "undrained excess pore generation ~ surcharge q (Skempton B~1, within ~10%)");
    check(C.consol_excess_pore.back() < 0.05 * q, "excess pore dissipated to ~0 by Tv ~ 2");

    // Settlement-time curve vs Terzaghi U(Tv) at representative points (avoid the steep early
    // boundary layer Tv < 0.15). U_FE = settlement(t) / s_inf.
    std::printf("   Tv     U_FE     U_Terzaghi   err\n");
    int checked = 0;
    for (size_t i = 1; i < C.consol_time.size(); ++i) {
        const double Tv = cv * C.consol_time[i] / (H * H);
        const double Ufe = C.consol_settlement[i] / s_inf;
        const double Uth = terzaghi_U(Tv);
        if (std::fabs(Tv - 0.2) < 0.013 || std::fabs(Tv - 0.4) < 0.013 ||
            std::fabs(Tv - 0.6) < 0.013 || std::fabs(Tv - 0.9) < 0.013) {
            std::printf("   %.3f  %.4f   %.4f       %.1f%%\n", Tv, Ufe, Uth, 100.0 * (Ufe - Uth) / Uth);
            check(std::fabs(Ufe - Uth) < 0.05, "U_FE matches Terzaghi U(Tv) within 5%");
            ++checked;
        }
    }
    check(checked >= 3, "checked several Tv points");

    // Final settlement consistent with Terzaghi U at the final Tv (the run does not reach Tv=inf).
    const double Tvf = cv * C.consol_time.back() / (H * H);
    const double Uf_expected = terzaghi_U(Tvf);
    const double Uf = C.consol_settlement.back() / s_inf;
    std::printf("   final: Tv=%.2f  U_FE=%.4f  U_Terzaghi=%.4f  (settlement %.5e of s_inf %.5e)\n",
                Tvf, Uf, Uf_expected, C.consol_settlement.back(), s_inf);
    check(std::fabs(Uf - Uf_expected) < 0.05, "final settlement consistent with Terzaghi U(Tv_final)");
}

}  // namespace

int main() {
    std::printf("Consolidation through the GUI compute path\n\n");
    test_terzaghi_gui(m::SoilModel::LinearElastic, "linear-elastic");
    test_terzaghi_gui(m::SoilModel::MohrCoulomb, "Mohr-Coulomb -> elastoplastic coupled-Newton path");
    if (g_failures == 0) {
        std::printf("\nOK: GUI consolidation path = Terzaghi 1D U-Tv + final settlement\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
