// "Water story" through the GUI compute path (build_problem + solve_phases). The cores are validated
// quantitatively elsewhere (test_transient_flow / test_unsaturated_flow / test_coupled_flow); here the
// INTEGRATED path -- project + materials (van Genuchten) + flow BCs -> mesh -> phases -- must work.
//   (1) FULLY-COUPLED flow-deformation phase on a saturated Terzaghi column reproduces U(Tv) (it
//       reduces to consolidation in the saturated regime) + reports a saturation field (~1).
//   (2) TRANSIENT flow phase (no deformation) on a column with prescribed-head boundaries runs,
//       converges, fills the pore/saturation fields, and resaturates an initially-dry interior.
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
    for (int j = 0; j < 80; ++j) { const double M = (2 * j + 1) * kPi / 2.0; s += (2.0 / (M * M)) * std::exp(-M * M * Tv); }
    return 1.0 - s;
}

// (1) Fully-coupled flow-deformation through the GUI = Terzaghi (saturated reduces to consolidation).
void test_fully_coupled_gui() {
    std::printf("-- Fully-coupled flow-deformation through the GUI phase path --\n");
    constexpr double W = 1.0, H = 12.0;
    constexpr double E = 1000.0, nu = 0.0;
    const double Eoed = E;
    constexpr double k = 0.1, q = 10.0;
    const double gw = katai::app::kGammaWater;
    const double cv = k * Eoed / gw;
    const double s_inf = q * H / Eoed;

    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = nu; s.gamma_unsat = 16.0; s.gamma_sat = 18.0; s.e_init = 0.5;
    s.kx = k; s.ky = k;   // van Genuchten defaults (g_a=2,g_n=2): saturated zone -> S_e=1
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};   // top drains
    P.edge_head = {0.0, 0.0, H, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 0; L.y1 = H; L.x2 = W; L.y2 = H; L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -q;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};

    m::Phase fc; fc.name = "Fully coupled"; fc.type = m::PhaseType::FullyCoupled;
    fc.duration = 2.0 * H * H / cv; fc.time_steps = 120; fc.load_active = {1};
    pr.phases.push_back(fc);

    const auto M = katai::app::mesh_from_project(pr, 0.08, 6);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("  (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + fully-coupled phases ran");
    if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); return; }
    const auto& C = res[1];
    check(!C.consol_time.empty(), "fully-coupled produced a settlement-time series");
    if (C.consol_time.empty()) return;
    check(!C.saturation.empty(), "fully-coupled reports a saturation field");

    // Saturated regime: S_eff ~ 1 everywhere (excess pore >= 0 -> psi <= 0).
    double smin = 2.0; for (double Sv : C.saturation) smin = std::fmin(smin, Sv);
    std::printf("  min saturation = %.4f (saturated ~1)  s_inf=%.4e  cv=%.4f\n", smin, s_inf, cv);
    check(smin > 0.98, "saturated zone S_eff ~ 1");

    std::printf("   Tv     U_FE     U_Terzaghi   err\n");
    int checked = 0;
    for (size_t i = 1; i < C.consol_time.size(); ++i) {
        const double Tv = cv * C.consol_time[i] / (H * H);
        const double Ufe = C.consol_settlement[i] / s_inf, Uth = terzaghi_U(Tv);
        if (std::fabs(Tv - 0.2) < 0.013 || std::fabs(Tv - 0.4) < 0.013 ||
            std::fabs(Tv - 0.6) < 0.013 || std::fabs(Tv - 0.9) < 0.013) {
            std::printf("   %.3f  %.4f   %.4f       %.1f%%\n", Tv, Ufe, Uth, 100.0 * (Ufe - Uth) / Uth);
            check(std::fabs(Ufe - Uth) < 0.05, "U_FE matches Terzaghi U(Tv) within 5%");
            ++checked;
        }
    }
    check(checked >= 3, "checked several Tv points");
}

// (2) Transient groundwater flow through the GUI (no deformation): pressurize an initially-dry column.
void test_transient_flow_gui() {
    std::printf("-- Transient groundwater flow through the GUI phase path --\n");
    constexpr double W = 1.0, H = 10.0;
    constexpr double k = 1.0, head_bc = H + 3.0;   // prescribed head above surface (artesian -> pore>0)

    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.0e4; s.nu = 0.3; s.e_init = 0.5; s.kx = k; s.ky = k;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    // Top + bottom prescribed head = H+3 (pressurize); sides closed. Interior fills from dry.
    P.edge_flow = {(int)m::FlowBCType::Head, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};
    P.edge_head = {head_bc, 0.0, head_bc, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;   // interior initial head = elevation (pressure head 0) -> partly unsaturated

    m::Phase tf; tf.name = "Transient flow"; tf.type = m::PhaseType::TransientFlow;
    tf.duration = 5.0; tf.time_steps = 40;
    pr.phases.push_back(tf);

    const auto M = katai::app::mesh_from_project(pr, 0.1, 6);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("  (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + transient-flow phases ran");
    if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); return; }
    const auto& T = res[1];
    check(!T.consol_time.empty(), "transient flow produced a time series");
    check(!T.saturation.empty(), "transient flow reports a saturation field");
    check((int)T.pore.size() == T.mesh.node_count, "transient flow reports a pore field");
    if (T.consol_time.empty() || T.saturation.empty()) return;

    // No deformation.
    check(T.max_disp == 0.0, "transient flow has no deformation");
    // Saturation bounded; pressurization -> pore develops; interior resaturates toward ~1.
    double smin = 2.0, smax = 0.0, pmax = 0.0;
    for (double Sv : T.saturation) { smin = std::fmin(smin, Sv); smax = std::fmax(smax, Sv); }
    for (double pv : T.pore) pmax = std::fmax(pmax, pv);
    double pore_peak = 0.0; for (double pv : T.consol_excess_pore) pore_peak = std::fmax(pore_peak, pv);
    std::printf("  saturation in [%.3f, %.3f]  final max pore = %.2f kPa  series peak = %.2f kPa\n",
                smin, smax, pmax, pore_peak);
    check(smin > 0.0 && smax <= 1.0 + 1e-9, "saturation bounded in (0,1]");
    check(smin > 0.98, "interior resaturated under pressurization (final S_eff ~ 1)");
    check(pmax > 0.0, "transient flow develops a pore-pressure field");
}

}  // namespace

int main() {
    std::printf("Water story through the GUI compute path\n\n");
    test_fully_coupled_gui();
    std::printf("\n");
    test_transient_flow_gui();
    if (g_failures == 0) {
        std::printf("\nOK: GUI water-story path (fully-coupled = Terzaghi + transient flow) verified\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
