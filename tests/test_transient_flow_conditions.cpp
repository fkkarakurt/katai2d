// The one-sided conditions and the moving water of a TRANSIENT groundwater flow phase -- against
// closed forms.
//
//   (A) A water boundary that rises at a constant rate. The phase starts from the head the
//       ground stands at and moves its Head edge linearly to the head of the phase. Oracle
//       (Carslaw & Jaeger 1959, sec. 2.5): a semi-infinite saturated column at a uniform head h0
//       whose end is raised as h0 + A t has
//           h(d, t) - h0 = 4 A t i2erfc(d / (2 sqrt(D t))),  D = k / S_s,
//           i2erfc(z) = [(1 + 2 z^2) erfc(z) - (2 / sqrt(pi)) z exp(-z^2)] / 4.
//       The same column with the change made at the first instant is the erfc front instead,
//       which is far from it: the ramp is applied, not the end value.
//   (B) A drawdown below the ground surface. The water over a saturated column falls below its
//       top, so the top becomes a seepage face: open while water leaves through it, at zero
//       pressure. The excess head over the face dissipates as Terzaghi's isochrone,
//           (h - y_top) / (h0 - y_top) = sum (2/M) sin(M Z) exp(-M^2 Tv),  M = (2m+1) pi / 2,
//       Z = depth / L, Tv = D t / L^2, and the ground ends at the face's level -- not at the head
//       of the water below the surface, which the edge, held at that head, would have pulled the
//       ground down to by suction through a face open to the air.
//   (C) A well with h_min. A column under a head h_top pumped at its base: by Darcy the base
//       stands Q L / (k W) below the top while the ground can supply the discharge Q, and at h_min,
//       with the discharge that the ground gives there, k W (h_top - h_min) / L, when it cannot.
//   (D) A drain only takes water out. Below the water it holds its head; in ground that is drier
//       than it, it does nothing.
//
// verify: KV-FLW-008
//   oracle:   closed_form
//   source:   Carslaw, H.S. & Jaeger, J.C. (1959). Conduction of Heat in Solids, 2nd ed. Oxford, sec. 2.5 (a surface value rising linearly in time, the i2erfc solution) and sec. 3.3 (the Fourier series of a slab with one face held, Terzaghi's isochrone); Darcy's law for the well and the drain
//   locator:  (A) the head at depths 0.5..4 m of a 30 m saturated column, k = 1e-4 m/day, its top raised from 40 to 45 m over 0.25 day; (B) the head at four depths of a 10 m column whose water falls from 15 m to 5 m (below its top at 10 m), at Tv = 0.31; (C) the head at the base and at mid-height of a 10 m column, top at 20 m, pumped 0.5 and 2 m3/day/m with h_min = 5 m; (D) the head above and below a drain at mid-height, ground water above and below the drain's head
//   quantity: hydraulic head [m]
//   expected: (A) within 3% of the rise A t of the i2erfc solution; (B) within 3% of the initial excess of the isochrone, the face at zero pressure; (C) the Darcy heads to 1e-6; (D) the Darcy head below the drain to 1e-6 and the hydrostatic head above it to 1 mm (that ground still drains at 100 days), and in drier ground no change at all
//   band:     3% for the diffusion fronts (time and space discretisation), 1e-6 for the steady states (the head field is linear, which 6-noded elements hold exactly)
#include <katai/analysis/phase_solver/transient_flow.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
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
// The specific storage the phase gives a saturated soil, S_s = n gamma_w / K_w.
double storage(double e) { return e / (1.0 + e) * katai::app::kGammaWater / 2.0e6; }

m::Project column(double W, double L, double k, double e) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.0e4; s.nu = 0.3; s.e_init = e; s.kx = k; s.ky = k;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, L, L};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed};
    P.edge_head = {0.0, 0.0, 0.0, 0.0};
    pr.polygons.push_back(P);
    pr.x_min = 0; pr.x_max = W; pr.y_min = 0; pr.y_max = L;
    return pr;
}
// The head at elevation y on the column's left side, linear between the nodes there (the column
// is one-dimensional, so a side is as good as any vertical).
double mean_head_at(const katai::app::SolveResult& R, double y) {
    std::vector<std::pair<double, double>> side;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.x[n]) < 1e-9) side.push_back({R.mesh.y[n], R.flow_head[n]});
    std::sort(side.begin(), side.end());
    for (std::size_t i = 1; i < side.size(); ++i)
        if (y <= side[i].first + 1e-12) {
            const double w = (y - side[i - 1].first) / (side[i].first - side[i - 1].first);
            return (1.0 - w) * side[i - 1].second + w * side[i].second;
        }
    return -1e300;
}

// (A) -------------------------------------------------------------------------------------------
void rising_water() {
    std::printf("-- (A) a water boundary rising at a constant rate --\n");
    constexpr double W = 1.0, L = 30.0, k = 1.0e-4, e = 0.5, h0 = 40.0, h1 = 45.0, T = 0.25;
    const double D = k / storage(e), A = (h1 - h0) / T;
    const m::Project pr = column(W, L, k, e);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "column meshed");
    if (!M.ok) return;
    const auto run = [&](bool ramp) {
        katai::core::TransientFlowPhase tin;
        katai::core::TransientFlowMaterial fm;
        fm.kx = fm.ky = k; fm.porosity = e / (1.0 + e);
        tin.materials = {fm};
        // Every edge, the closed ones included: a boundary node takes the kind of the edge
        // nearest to it, so an edge left out would hand its nodes to the next one.
        katai::core::FlowEdge base{0, 0, W, 0}, right{W, 0, W, L}, top{W, L, 0, L}, left{0, L, 0, 0};
        top.kind = katai::core::FlowEdgeKind::Head; top.head = h1;
        tin.flow_edges = {base, right, top, left};
        tin.initial_head.assign(M.mesh.node_count, h0);
        tin.fallback_head.assign(M.mesh.node_count, h0);
        tin.ramp_heads = ramp;
        tin.duration_day = T; tin.time_steps = 100;
        katai::app::SolveResult R;
        R.mesh = M.mesh;
        const bool ok = katai::core::solve_transient_flow_phase(M.mesh, tin, R);
        return std::make_pair(ok, R);
    };
    const auto [ok, R] = run(true);
    check(ok, "the phase ran");
    if (!ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    const auto [ok_step, S] = run(false);
    double worst = 0.0, gap_to_step = 0.0;
    for (const double d : {0.5, 1.0, 2.0, 3.0, 4.0}) {
        const double z = d / (2.0 * std::sqrt(D * T));
        const double i2 = ((1.0 + 2.0 * z * z) * std::erfc(z) - 2.0 / std::sqrt(kPi) * z * std::exp(-z * z)) / 4.0;
        const double want = h0 + 4.0 * A * T * i2;
        const double got = mean_head_at(R, L - d);
        const double step = mean_head_at(S, L - d);
        std::printf("  d = %.1f m: h = %.4f m, closed form %.4f m (the change at once: %.4f m)\n", d, got, want, step);
        worst = std::fmax(worst, std::fabs(got - want) / (A * T));
        gap_to_step = std::fmax(gap_to_step, std::fabs(step - want) / (A * T));
    }
    std::printf("  worst %.2f%% of the rise; the change made at once is %.1f%% away\n", 100 * worst, 100 * gap_to_step);
    check(worst < 0.03, "the head follows the i2erfc solution within 3% of the rise");
    check(ok_step && gap_to_step > 0.10, "and the change made at once does not (the ramp is applied)");
}

// (B) -------------------------------------------------------------------------------------------
void drawdown_below_the_surface() {
    std::printf("-- (B) a drawdown below the ground surface: the top becomes a seepage face --\n");
    constexpr double W = 1.0, L = 10.0, k = 1.0e-4, e = 0.5, h0 = 15.0, hlow = 5.0, T = 0.5;
    const double D = k / storage(e), Tv = D * T / (L * L);
    m::Project pr = column(W, L, k, e);
    pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {h0, h0};
    pr.polygons[0].edge_flow[2] = (int)m::FlowBCType::Head;
    pr.polygons[0].edge_head[2] = hlow;
    m::Phase tf; tf.name = "Drawdown"; tf.type = m::PhaseType::TransientFlow;
    tf.duration = T; tf.time_steps = 200;
    pr.phases.push_back(tf);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2 && res[1].ok, "the drawdown phase ran");
    if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); return; }
    const auto& R = res[1];
    double worst = 0.0;
    for (const double Z : {0.25, 0.5, 0.75, 1.0}) {
        double s = 0.0;
        for (int j = 0; j < 200; ++j) {
            const double Mj = (2 * j + 1) * kPi / 2.0;
            s += 2.0 / Mj * std::sin(Mj * Z) * std::exp(-Mj * Mj * Tv);
        }
        const double want = L + (h0 - L) * s;
        const double got = mean_head_at(R, L - Z * L);
        std::printf("  Z = %.2f: h = %.4f m, isochrone %.4f m (Tv = %.3f)\n", Z, got, want, Tv);
        worst = std::fmax(worst, std::fabs(got - want) / (h0 - L));
    }
    const double top = mean_head_at(R, L);
    std::printf("  the face stands at h = %.6f m (its level %.1f m; the water below it %.1f m)\n", top, L, hlow);
    check(worst < 0.03, "the excess over the face dissipates as Terzaghi's isochrone (3%)");
    check(std::fabs(top - L) < 1e-6, "the face is at zero pressure, not at the head of the water below it");
}

// (C) -------------------------------------------------------------------------------------------
void well_with_h_min() {
    std::printf("-- (C) a well pumped at its rated discharge, and limited by its h_min --\n");
    constexpr double W = 1.0, L = 10.0, k = 1.0, e = 0.5, htop = 20.0, hmin = 5.0;
    for (const double Q : {0.5, 2.0}) {
        m::Project pr = column(W, L, k, e);
        pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {htop, htop};
        pr.polygons[0].edge_flow[2] = (int)m::FlowBCType::Head;
        pr.polygons[0].edge_head[2] = htop;
        m::HydroLine w; w.name = "Well"; w.kind = m::HydroKind::Well;
        w.behaviour = (int)m::WellBehaviour::Extraction;
        w.x1 = 0.0; w.y1 = 0.0; w.x2 = W; w.y2 = 0.0; w.q = Q; w.h_min = hmin;
        pr.hydros.push_back(w);
        pr.initial.hydro_active = {0};
        m::Phase tf; tf.name = "Pumping"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 1.0; tf.time_steps = 10;
        pr.phases.push_back(tf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        check(res.size() == 2 && res[1].ok, "the pumping phase ran");
        if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); continue; }
        const double free_base = htop - Q * L / (k * W);
        const bool limited = free_base < hmin;
        const double base = limited ? hmin : free_base;
        const double want_mid = 0.5 * (htop + base);
        const double got_base = mean_head_at(res[1], 0.0), got_mid = mean_head_at(res[1], 0.5 * L);
        std::printf("  Q = %.1f: base %.6f m (want %.6f, %s), mid %.6f m (want %.6f)\n", Q, got_base, base,
                    limited ? "held at h_min" : "free", got_mid, want_mid);
        check(std::fabs(got_base - base) < 1e-6 * htop && std::fabs(got_mid - want_mid) < 1e-6 * htop,
              limited ? "beyond what the ground can give: the base at h_min, the head linear to it"
                      : "within what the ground can give: the base Q L / (k W) below the top");
    }
}

// (D) -------------------------------------------------------------------------------------------
void drain_only_takes_water_out() {
    std::printf("-- (D) a drain at mid-height --\n");
    constexpr double W = 1.0, L = 10.0, k = 1.0, e = 0.5, hd = 6.0;
    for (const double hb : {12.0, 4.0}) {
        m::Project pr = column(W, L, k, e);
        pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {hb, hb};
        pr.polygons[0].edge_flow[0] = (int)m::FlowBCType::Head;
        pr.polygons[0].edge_head[0] = hb;
        m::HydroLine d; d.name = "Drain"; d.kind = m::HydroKind::Drain;
        d.x1 = 0.0; d.y1 = 0.5 * L; d.x2 = W; d.y2 = 0.5 * L; d.head = hd;
        pr.hydros.push_back(d);
        pr.initial.hydro_active = {0};
        // A silt's retention curve (alpha 1 1/m, n 2). The phase is long enough to be steady above
        // the drain too, where the drier ground releases its water slowly.
        pr.materials[0].gw_ga = 1.0; pr.materials[0].gw_gn = 2.0;
        m::Phase tf; tf.name = "Drained"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 100.0; tf.time_steps = 50;
        pr.phases.push_back(tf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        check(res.size() == 2 && res[1].ok, "the drain phase ran");
        if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); continue; }
        const double above = mean_head_at(res[1], 0.75 * L), below = mean_head_at(res[1], 0.25 * L);
        if (hb > hd) {
            std::printf("  ground water at %.0f m, drain at %.0f m: above %.6f (want %.6f), below %.6f (want %.6f)\n",
                        hb, hd, above, hd, below, 0.5 * (hb + hd));
            // Below the drain the ground is saturated and the head is steady at once; above it the
            // drier ground is still releasing its last water at 100 days, so it is held to a
            // millimetre of the hydrostatic head it tends to.
            check(std::fabs(below - 0.5 * (hb + hd)) < 1e-6 * hb && std::fabs(above - hd) < 1e-3,
                  "below the water the drain holds its head: linear below it, hydrostatic above (1 mm)");
        } else {
            std::printf("  ground water at %.0f m, drain at %.0f m: above %.6f, below %.6f (want %.6f both)\n",
                        hb, hd, above, below, hb);
            check(std::fabs(above - hb) < 1e-6 * hd && std::fabs(below - hb) < 1e-6 * hd,
                  "in drier ground the drain does nothing: it feeds no water in");
        }
    }
    // THE DEFAULT SAND (alpha 14.5 1/m, n 2.68). The ground the drain dries is almost dry four
    // metres above the water, and there the Picard iteration did not converge within its budget
    // at all -- a 6 m drop asked for a 1142 m update. Damped, capped at a metre an update, and
    // given four times the iterations once a step has been split, it does. Its unsaturated
    // conductivity is then a few ten-thousandths of the saturated one, so the ground above the
    // drain is still draining at 100 days: the check there is that it has drained, not that it
    // has finished.
    {
        constexpr double hb = 12.0;
        m::Project pr = column(W, L, k, e);
        pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {hb, hb};
        pr.polygons[0].edge_flow[0] = (int)m::FlowBCType::Head;
        pr.polygons[0].edge_head[0] = hb;
        m::HydroLine d; d.name = "Drain"; d.kind = m::HydroKind::Drain;
        d.x1 = 0.0; d.y1 = 0.5 * L; d.x2 = W; d.y2 = 0.5 * L; d.head = hd;
        pr.hydros.push_back(d);
        pr.initial.hydro_active = {0};
        m::Phase tf; tf.name = "Drained"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = 100.0; tf.time_steps = 50;
        pr.phases.push_back(tf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        check(res.size() == 2 && res[1].ok, "on the default sand the drain phase runs (it did not converge)");
        if (res.size() == 2 && res[1].ok) {
            const double above = mean_head_at(res[1], 0.75 * L), below = mean_head_at(res[1], 0.25 * L);
            std::printf("  sand: above %.6f (draining from %.0f towards %.0f), below %.6f (want %.6f)\n",
                        above, hb, hd, below, 0.5 * (hb + hd));
            check(std::fabs(below - 0.5 * (hb + hd)) < 1e-6 * hb && above > hd && above < hb,
                  "below the drain the Darcy head; above it the ground has drained part of the way");
        }
    }
}
// (E) -------------------------------------------------------------------------------------------
// The same ramp through the file contract and the phase sequence. Each phase states its water
// level, which its Head edges take as their head: the first raises the column to 20 m while the
// project's water table stays at 15 m, so the second can only match the closed form if it starts
// from the head the first one ended with; the second lowers its level to 10 m with `headramp`.
void ramp_through_the_phase_sequence() {
    std::printf("-- (E) a phase sequence: the level raised, then lowered at a rate --\n");
    constexpr double W = 1.0, L = 10.0, k = 1.0e-4, e = 0.5, hw = 15.0, hup = 20.0, hdown = 10.0, T = 0.02;
    const double D = k / storage(e), A = (hup - hdown) / T;
    m::Project pr = column(W, L, k, e);
    pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {hw, hw};
    pr.polygons[0].edge_flow[2] = (int)m::FlowBCType::Head;
    pr.polygons[0].edge_head[2] = hw;
    m::Phase rise; rise.name = "Rise"; rise.type = m::PhaseType::TransientFlow;
    rise.duration = 50.0; rise.time_steps = 20;
    rise.water_override = true; rise.wx = {0.0, W}; rise.wy = {hup, hup};
    pr.phases.push_back(rise);
    m::Phase fall = rise; fall.name = "Fall"; fall.duration = T; fall.time_steps = 80;
    fall.wy = {hdown, hdown}; fall.head_ramp = true;
    pr.phases.push_back(fall);
    const auto M = katai::app::mesh_from_project(pr, 0.25, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[1].ok && res[2].ok, "the rise and the fall ran");
    if (res.size() != 3 || !res[2].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); return; }
    std::printf("  after the rise: %.4f m at the base (the project's water table is at %.0f m)\n",
                mean_head_at(res[1], 0.0), hw);
    check(std::fabs(mean_head_at(res[1], 0.0) - hup) < 1e-6, "the rise takes its Head edge to the phase's level");
    double worst = 0.0;
    for (const double d : {0.25, 0.5, 1.0}) {
        const double z = d / (2.0 * std::sqrt(D * T));
        const double i2 = ((1.0 + 2.0 * z * z) * std::erfc(z) - 2.0 / std::sqrt(kPi) * z * std::exp(-z * z)) / 4.0;
        const double want = hup - 4.0 * A * T * i2;
        const double got = mean_head_at(res[2], L - d);
        std::printf("  d = %.2f m: h = %.4f m, closed form %.4f m\n", d, got, want);
        worst = std::fmax(worst, std::fabs(got - want) / (A * T));
    }
    std::printf("  worst %.2f%% of the fall\n", 100 * worst);
    check(worst < 0.03, "the fall follows the i2erfc solution from the head the rise ended with (3%)");
}
// (F) -------------------------------------------------------------------------------------------
// An excavated region holds no water. A column of two stacked regions, the upper one excavated, is
// drained from its base; the same drainage on the lower region meshed alone is the independent
// path. Until 2026-10-23 every element took part in a transient flow, and the excavated soil went
// on storing water that the drainage then had to remove.
void excavated_ground_holds_no_water() {
    std::printf("-- (F) a transient flow after an excavation: the excavated soil takes no part --\n");
    constexpr double W = 1.0, Lb = 5.0, k = 1.0e-4, e = 0.5, h0 = 5.0, hdown = 3.0, T = 0.1;
    const auto lower = [&](double top) {
        m::Project pr = column(W, top, k, e);
        pr.materials[0].gw_ga = 1.0; pr.materials[0].gw_gn = 2.0;   // a silt: the soil above the water holds some
        pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {h0, h0};
        pr.polygons[0].y = {0, 0, Lb, Lb};
        pr.polygons[0].edge_flow[0] = (int)m::FlowBCType::Head;
        pr.polygons[0].edge_head[0] = hdown;
        return pr;
    };
    m::Project alone = lower(Lb);
    m::Project stacked = lower(2.0 * Lb);
    m::SoilPolygon up = stacked.polygons[0];
    up.name = "Excavated"; up.y = {Lb, Lb, 2.0 * Lb, 2.0 * Lb};
    up.edge_flow = {0, 0, 0, 0};
    stacked.polygons.push_back(up);
    // The phase lowers its water to the base's head, below the floor the excavation leaves: that
    // floor is then a seepage face (the faces an excavation exposes take the phase's water), which
    // stays closed while the ground drains -- the closed top of the column meshed alone.
    for (m::Project* pr : {&alone, &stacked}) {
        m::Phase tf; tf.name = "Drain"; tf.type = m::PhaseType::TransientFlow;
        tf.duration = T; tf.time_steps = 50;
        tf.water_override = true; tf.wx = {0.0, W}; tf.wy = {hdown, hdown};
        if (pr == &stacked) { stacked.initial.poly_active = {1, 1}; tf.poly_active = {1, 0}; }
        pr->phases.push_back(tf);
    }
    const auto Ma = katai::app::mesh_from_project(alone, 0.5, 6);
    const auto Ms = katai::app::mesh_from_project(stacked, 0.5, 6);
    const auto ra = katai::app::solve_phases(alone, Ma.mesh, InitialPhase::K0Procedure);
    const auto rs = katai::app::solve_phases(stacked, Ms.mesh, InitialPhase::K0Procedure);
    check(ra.size() == 2 && ra[1].ok && rs.size() == 2 && rs[1].ok, "both drainages ran");
    if (ra.size() != 2 || !ra[1].ok || rs.size() != 2 || !rs[1].ok) {
        if (!rs.empty()) std::printf("  (%s)\n", rs.back().message.c_str());
        return;
    }
    double worst = 0.0;
    for (const double y : {0.5, 1.5, 2.5, 3.5, 4.5}) {
        const double a = mean_head_at(ra[1], y), b = mean_head_at(rs[1], y);
        std::printf("  y = %.1f m: alone %.6f m, under the excavation %.6f m\n", y, a, b);
        worst = std::fmax(worst, std::fabs(a - b));
    }
    check(worst < 1e-4, "the ground under the excavation drains as the same ground alone (0.1 mm)");
    // The control: with the upper region left in place its soil does take part, holding water the
    // drainage has to remove, and the heads under it are not those of the ground alone.
    m::Project kept = stacked;
    kept.phases[0].poly_active = {1, 1};
    const auto rk = katai::app::solve_phases(kept, Ms.mesh, InitialPhase::K0Procedure);
    double gap = 0.0;
    if (rk.size() == 2 && rk[1].ok)
        for (const double y : {0.5, 1.5, 2.5, 3.5, 4.5})
            gap = std::fmax(gap, std::fabs(mean_head_at(ra[1], y) - mean_head_at(rk[1], y)));
    std::printf("  with the upper region left in place: %.2e m from the ground alone\n", gap);
    check(gap > 10.0 * 1e-4, "and the same model with the region kept differs by far more (the check discriminates)");
}
// (G) -------------------------------------------------------------------------------------------
// A pit pumped down: its faces take the phase's water, and a transient flow run to its steady state
// is the steady flow of a staged phase -- a separate solver (the free-surface seepage solve of
// `flowpore`) on the same ground and the same water. The block is 20 x 10 m, the pit 4 m wide and
// 4 m deep; the phase's water stands at 9 m at the model's sides and at the pit floor (6 m) in
// the pit, so water flows under the pit's walls into its floor.
void pit_faces_reach_the_steady_flow(bool sand) {
    std::printf("-- (G) a pit's faces in a transient flow, run to the steady flow of a staged phase (%s) --\n",
                sand ? "the default sand: 50-day steps that have to be split" : "a silt");
    constexpr double Wb = 20.0, Hb = 10.0, x0 = 8.0, x1 = 12.0, yf = 6.0, hw = 9.0, k = 1.0, e = 0.5;
    const auto base = [&]() {
        m::Project pr = column(Wb, Hb, k, e);
        if (!sand) { pr.materials[0].gw_ga = 1.0; pr.materials[0].gw_gn = 2.0; }
        m::SoilPolygon& P = pr.polygons[0];
        P.x = {0, Wb, Wb, x1, x1, x0, x0, 0};
        P.y = {0, 0, Hb, Hb, yf, yf, Hb, Hb};
        P.edge_bc.assign(8, 0); P.edge_bc[0] = (int)m::BCType::FullyFixed;
        P.edge_bc[1] = P.edge_bc[7] = (int)m::BCType::HorizontallyFixed;
        // The model's sides hold the water; its ground surface is a seepage face; the base is closed.
        P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head, (int)m::FlowBCType::Seepage,
                       (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                       (int)m::FlowBCType::Seepage, (int)m::FlowBCType::Head};
        P.edge_head.assign(8, hw);
        m::SoilPolygon pit; pit.name = "Pit"; pit.material = 0;
        pit.x = {x0, x1, x1, x0}; pit.y = {yf, yf, Hb, Hb};
        pit.edge_bc.assign(4, 0); pit.edge_flow.assign(4, 0); pit.edge_head.assign(4, 0.0);
        pr.polygons.push_back(pit);
        pr.has_water = true; pr.wx = {0.0, Wb}; pr.wy = {hw, hw};
        return pr;
    };
    m::Phase dig; dig.name = "Pit"; dig.poly_active = {1, 0};
    dig.water_override = true;
    dig.wx = {0.0, x0 - 1e-3, x0, x1, x1 + 1e-3, Wb};   // the step just outside the walls: they drain
    dig.wy = {hw, hw, yf, yf, hw, hw};
    m::Project tr = base();
    m::Phase t = dig; t.type = m::PhaseType::TransientFlow; t.duration = 2000.0; t.time_steps = 40;
    tr.phases.push_back(t);
    m::Project st = base();
    m::Phase p = dig; p.type = m::PhaseType::Plastic; p.pore_from_flow = true;
    st.phases.push_back(p);
    const auto M = katai::app::mesh_from_project(tr, 0.5, 6);
    check(M.ok, "pit meshed");
    if (!M.ok) return;
    const auto rt = katai::app::solve_phases(tr, M.mesh, InitialPhase::K0Procedure);
    const auto rs = katai::app::solve_phases(st, M.mesh, InitialPhase::K0Procedure);
    check(rt.size() == 2 && rt[1].ok && rs.size() == 2 && rs[1].ok, "the transient and the steady flow ran");
    if (rt.size() != 2 || !rt[1].ok || rs.size() != 2 || !rs[1].ok) {
        if (!rt.empty()) std::printf("  transient: %s\n", rt.back().message.c_str());
        if (!rs.empty()) std::printf("  steady: %s\n", rs.back().message.c_str());
        return;
    }
    const auto& A = rt[1];
    const auto& B = rs[1];
    check(A.flow_head.size() == B.flow_head.size(), "both carry a head on the same nodes");
    if (A.flow_head.size() != B.flow_head.size()) return;
    // Below the pit floor (saturated ground, where the pore pressure lives): the heads agree.
    double worst = 0.0, at_x = 0.0, at_y = 0.0;
    for (int n = 0; n < A.mesh.node_count; ++n) {
        if (A.mesh.y[n] > yf + 1e-9) continue;
        const double d = std::fabs(A.flow_head[n] - B.flow_head[n]);
        if (d > worst) { worst = d; at_x = A.mesh.x[n]; at_y = A.mesh.y[n]; }
    }
    double floor_head = 0.0; int nf = 0;
    for (int n = 0; n < A.mesh.node_count; ++n)
        if (std::fabs(A.mesh.y[n] - yf) < 1e-9 && A.mesh.x[n] > x0 + 1e-9 && A.mesh.x[n] < x1 - 1e-9) {
            floor_head += A.flow_head[n]; ++nf;
        }
    floor_head = nf ? floor_head / nf : 0.0;
    std::printf("  pit floor head %.6f m (its water %.1f m); below the floor the transient and steady heads "
                "differ by at most %.2e m, at (%.2f, %.2f)\n", floor_head, yf, worst, at_x, at_y);
    check(std::fabs(floor_head - yf) < 1e-6, "the pit floor holds the pit's water");
    check(worst < 0.01 * (hw - yf), "the transient flow reaches the steady flow of the staged phase (1% of the head drop)");
}
} // namespace

int main() {
    rising_water();
    drawdown_below_the_surface();
    well_with_h_min();
    drain_only_takes_water_out();
    ramp_through_the_phase_sequence();
    excavated_ground_holds_no_water();
    pit_faces_reach_the_steady_flow(false);
    pit_faces_reach_the_steady_flow(true);
    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all transient flow condition checks passed\n");
    return 0;
}
