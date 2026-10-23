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

// (3) A PRESCRIBED INFLOW through the GUI path. Until 2026-10-22 a Flux edge reached every
// time-dependent phase as a Closed one: a rainfall infiltration was silently zero. Oracle (closed
// form, Darcy): a saturated column of length L, permeability k, head h_b held at its base and an
// inflow q entering at its top reaches the steady state h(y) = h_b + q (y - y_b) / k, so the top
// stands q L / k above the base and its pore pressure is gamma_w (h_b + q L / k - y_top). Its
// diffusion time L^2 S_s / k is a ten-thousandth of a day, so a one-day phase ends at the
// steady state; the head field is linear, which 6-noded elements hold exactly.
void test_transient_inflow_gui() {
    std::printf("-- A prescribed inflow in a transient flow phase --\n");
    constexpr double W = 1.0, H = 10.0, k = 1.0, q = 0.05, hb = H + 5.0;
    const double gw = katai::app::kGammaWater;
    auto column = [&](m::FlowBCType top) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 1.0e4; s.nu = 0.3; s.e_init = 0.5; s.kx = k; s.ky = k;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
        P.edge_flow = {(int)m::FlowBCType::Head, (int)m::FlowBCType::Closed, (int)top, (int)m::FlowBCType::Closed};
        P.edge_head = {hb, 0.0, 0.0, 0.0};
        P.edge_flux = {0.0, 0.0, top == m::FlowBCType::Flux ? q : 0.0, 0.0};
        pr.polygons.push_back(P);
        pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {hb, hb};   // saturated from the start
        return pr;
    };
    auto top_pore = [&](const katai::app::SolveResult& R) {
        double p = 0.0; int n_top = 0;
        for (int n = 0; n < R.mesh.node_count; ++n)
            if (std::fabs(R.mesh.y[n] - H) < 1e-9) { p += R.pore[n]; ++n_top; }
        return n_top ? p / n_top : -1.0;
    };
    {
        m::Project pr = column(m::FlowBCType::Flux);
        m::Phase tf; tf.name = "Rain"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 1.0; tf.time_steps = 20;
        pr.phases.push_back(tf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        check(M.ok, "column meshed");
        if (!M.ok) return;
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        check(res.size() == 2 && res[1].ok, "transient flow with an inflow edge solved");
        if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); return; }
        const double want = gw * (hb + q * H / k - H);
        const double got = top_pore(res[1]);
        std::printf("  top pore pressure %.6f kPa, closed form %.6f kPa (inflow ignored: %.6f)\n",
                    got, want, gw * (hb - H));
        check(std::fabs(got - want) < 1e-6 * want, "the inflow is applied: top pore pressure = gamma_w (h_b + qL/k - H)");
    }
    // The same edge in a consolidation phase is said, not dropped (K2D-A033).
    {
        m::Project pr = column(m::FlowBCType::Flux);
        m::Phase c; c.name = "Consolidation"; c.type = m::PhaseType::Consolidation;
        c.duration = 1.0; c.time_steps = 5;
        pr.phases.push_back(c);
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool said = false;
        if (res.size() == 2)
            for (const auto& d : res[1].diagnostics) said = said || d.code == "K2D-A033";
        check(res.size() == 2 && said, "a consolidation phase says its flux edge is not applied (K2D-A033)");
    }
    // A SEEPAGE FACE in a transient flow phase is open. Until 2026-10-23 it was closed there (and
    // said so, K2D-A034). Oracle (Darcy): the base held at h_b = H + 5 m, the top a seepage face,
    // so water rises through the column and leaves at its top at zero pressure; the head is then
    // linear from h_b at the base to H at the top. Closed, the top would stand at h_b (49 kPa).
    {
        m::Project pr = column(m::FlowBCType::Seepage);
        m::Phase tf; tf.name = "Transient"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 1.0; tf.time_steps = 10;
        pr.phases.push_back(tf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        check(res.size() == 2 && res[1].ok, "transient flow with a seepage face solved");
        if (res.size() == 2 && res[1].ok) {
            const auto& T = res[1];
            double mid = 0.0; int n_mid = 0;
            for (int n = 0; n < T.mesh.node_count; ++n)
                if (std::fabs(T.mesh.y[n] - 0.5 * H) < 1e-9) { mid += T.pore[n]; ++n_mid; }
            mid = n_mid ? mid / n_mid : -1.0;
            const double want_mid = gw * (0.5 * (hb + H) - 0.5 * H);
            std::printf("  top pore %.3e kPa (closed: %.3f), mid-height %.6f kPa, closed form %.6f kPa\n",
                        top_pore(T), gw * (hb - H), mid, want_mid);
            check(std::fabs(top_pore(T)) < 1e-6 * gw * (hb - H), "the face is open: zero pressure at the top");
            check(std::fabs(mid - want_mid) < 1e-6 * want_mid, "the head is linear from the base to the face (Darcy)");
            bool said = false;
            for (const auto& d : T.diagnostics) said = said || d.code == "K2D-A034";
            check(!said, "K2D-A034 is retired: not raised");
        }
    }
}

// (4) A WALL'S SEAM IN A TRANSIENT FLOW PHASE. A wall with interfaces splits the mesh; the
// coupled phases tie the two sides of a permeable seam to one pore equation, but the transient
// flow solver was handed the split mesh as it was, so every seam -- fully permeable included --
// stopped the water. Oracle (closed form): a confined box with the head held at h_L on its left
// side and h_R on its right, base and top closed, reaches h(x) = h_L + (h_R - h_L) x / L, which
// 6-noded elements hold exactly; a fully permeable wall inside it changes nothing. An impermeable
// wall must change the field (the water goes under it), or the flag is not read either way.
void test_transient_wall_seam_gui() {
    std::printf("-- A wall's seam in a transient flow phase --\n");
    constexpr double L = 20.0, H = 10.0, hL = 14.0, hR = 12.0, k = 1.0;
    auto box = [&](int barrier, bool with_wall) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 1.0e4; s.nu = 0.3; s.e_init = 0.5; s.kx = k; s.ky = k;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, L, L, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
        P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head,
                       (int)m::FlowBCType::Closed, (int)m::FlowBCType::Head};
        P.edge_head = {0.0, hR, 0.0, hL};
        pr.polygons.push_back(P);
        pr.has_water = true; pr.wx = {0.0, L}; pr.wy = {13.0, 13.0};   // saturated everywhere
        if (with_wall) {
            m::PlateMaterial pm; pm.EA = 1.0e6; pm.EI = 1.0e4;
            pr.plates.push_back(pm);
            m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall"; w.material = 0;
            w.x1 = 10.0; w.y1 = H; w.x2 = 10.0; w.y2 = 4.0;
            w.iface_pos = w.iface_neg = true;
            w.flow_barrier = barrier;
            pr.structs.push_back(w);
        }
        m::Phase tf; tf.name = "Flow"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 10.0; tf.time_steps = 10;
        pr.phases.push_back(tf);
        return pr;
    };
    // The largest departure from the linear field, and the largest head jump between the two
    // sides of the wall (nodes at the same position).
    struct Out { bool ok = false; double dev = 0.0, jump = 0.0; std::string msg; };
    auto run = [&](const m::Project& pr) {
        Out o;
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) { o.msg = M.message; return o; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.size() != 2 || !res[1].ok) { o.msg = res.empty() ? "" : res.back().message; return o; }
        const auto& R = res[1];
        o.ok = true;
        const double gw = katai::app::kGammaWater;
        for (int n = 0; n < R.mesh.node_count; ++n) {
            const double h = R.pore[n] / gw + R.mesh.y[n];
            o.dev = std::fmax(o.dev, std::fabs(h - (hL + (hR - hL) * R.mesh.x[n] / L)));
            for (int j = n + 1; j < R.mesh.node_count; ++j)
                if (std::fabs(R.mesh.x[j] - R.mesh.x[n]) < 1e-9 && std::fabs(R.mesh.y[j] - R.mesh.y[n]) < 1e-9)
                    o.jump = std::fmax(o.jump, std::fabs(R.pore[j] - R.pore[n]));
        }
        return o;
    };
    const Out plain = run(box(0, false));
    const Out perm = run(box(0, true));
    const Out imp = run(box(1, true));
    std::printf("  no wall: max |h - linear| %.2e m | permeable wall: %.2e m, seam jump %.2e kPa | "
                "impermeable wall: %.2e m, seam jump %.3f kPa\n",
                plain.dev, perm.dev, perm.jump, imp.dev, imp.jump);
    if (!perm.ok) std::printf("  (%s)\n", perm.msg.c_str());
    check(plain.ok && plain.dev < 1e-6, "no wall: the linear head field of the closed form");
    check(perm.ok && perm.dev < 1e-6 && perm.jump < 1e-9,
          "a FULLY PERMEABLE wall leaves the field as it is: its seam is one pore pressure");
    check(imp.ok && imp.dev > 0.1 && imp.jump > 1.0,
          "an IMPERMEABLE wall changes it: the head jumps across the seam");
}

}  // namespace

int main() {
    std::printf("Water story through the GUI compute path\n\n");
    test_fully_coupled_gui();
    std::printf("\n");
    test_transient_flow_gui();
    std::printf("\n");
    test_transient_inflow_gui();
    std::printf("\n");
    test_transient_wall_seam_gui();
    if (g_failures == 0) {
        std::printf("\nOK: GUI water-story path (fully-coupled = Terzaghi + transient flow) verified\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
