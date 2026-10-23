// A drawdown's pore pressures reach the deformation: a phase that takes its pore pressures from a
// groundwater flow, right after a transient flow phase, takes the head that flow ended with.
//
// A transient flow phase moves the water and nothing else; the ground feels the drawdown in the
// phase after it. That phase used to solve a steady flow of its own -- or take the water table --
// whatever the drawdown before it had left, so the effective stress of a lowered reservoir's slope
// was that of the water table, not of the drawdown. It now takes the drawdown's head, and records
// the head it used, so a Safety phase after it reads the same field.
//
// Oracle (one-dimensional, closed form): a saturated column of height L, its water table at the
// top, is drained from its base to a head h_d. When the flow has settled the water table stands at
// h_d: below it the pore pressure has fallen by gamma_w (L - h_d), above it to zero from
// gamma_w (L - y). With a single unit weight above and below the water table the total stress does
// not change, so the effective stress rises by exactly the pore pressure lost, and with nu = 0 the
// top settles
//     s = (gamma_w / E) [ (L - h_d) h_d + (L - h_d)^2 / 2 ].
// The effective stress varies linearly in y, so 6-noded elements hold the displacement exactly.
//
// verify: KV-FLW-009
//   oracle:   closed_form
//   source:   Terzaghi's effective stress principle with one-dimensional (oedometric, nu = 0) linear elasticity: a pore pressure change du under constant total stress changes the vertical effective stress by -du, and the settlement is the integral of -du / E over the height
//   locator:  a 1 x 10 m linear elastic column (E = 1e4 kPa, nu = 0, one unit weight above and below the water), its water table at the top; a transient flow phase lowers the base to 5 m over 400 days; the Plastic phase after it takes its pore pressures from that flow
//   quantity: nodal pore pressure [kPa]; settlement of the top [m]
//   expected: the Plastic phase's pore pressures equal those the transient flow ended with, node by node, and the top settles (gamma_w / E)(25 + 12.5) = 0.0368 m
//   band:     1e-9 relative on the pore pressures handed over (the same field); 1% on the settlement (the unsaturated zone above the new water table still drains at 400 days)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
}  // namespace

int main() {
    constexpr double W = 1.0, L = 10.0, E = 1.0e4, hd = 5.0;
    const double gw = katai::app::kGammaWater;
    m::Project pr;
    m::Material s; s.name = "Silt"; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = 0.0; s.gamma_unsat = 20.0; s.gamma_sat = 20.0; s.e_init = 0.5;
    s.kx = s.ky = 1.0; s.gw_ga = 1.0; s.gw_gn = 2.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, L, L};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Head, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed};
    P.edge_head = {L, 0.0, 0.0, 0.0};
    pr.polygons.push_back(P);
    pr.x_min = 0; pr.x_max = W; pr.y_min = 0; pr.y_max = L;
    pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {L, L};

    m::Phase drain; drain.name = "Drawdown"; drain.type = m::PhaseType::TransientFlow;
    drain.duration = 400.0; drain.time_steps = 40;
    drain.water_override = true; drain.wx = {0.0, W}; drain.wy = {hd, hd};
    pr.phases.push_back(drain);
    m::Phase load; load.name = "Effective stress"; load.type = m::PhaseType::Plastic;
    load.pore_from_flow = true;
    load.water_override = true; load.wx = {0.0, W}; load.wy = {hd, hd};
    pr.phases.push_back(load);

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "column meshed");
    if (!M.ok) return 1;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[1].ok && res[2].ok, "initial, drawdown and Plastic phases ran");
    if (res.size() != 3 || !res[2].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str());
        return 1;
    }
    const auto& T = res[1];
    const auto& Pl = res[2];
    double worst = 0.0, scale = 1e-30;
    const bool same_size = Pl.pore.size() == T.pore.size();
    if (same_size)
        for (size_t n = 0; n < T.pore.size(); ++n) {
            worst = std::fmax(worst, std::fabs(Pl.pore[n] - T.pore[n]));
            scale = std::fmax(scale, std::fabs(T.pore[n]));
        }
    std::printf("  pore pressures handed over: worst difference %.3e kPa (largest %.3f kPa)\n", worst, scale);
    check(same_size && worst <= 1e-9 * scale, "the Plastic phase's pore pressures are those the drawdown ended with");
    bool said = false;
    for (const auto& d : Pl.diagnostics) said = said || d.code == "K2D-A038";
    check(said, "and it says where they came from (K2D-A038)");
    check(Pl.flow_head.size() == Pl.pore.size(), "it records the head it used, for a phase after it");

    double top = 0.0; int ntop = 0;
    for (int n = 0; n < Pl.mesh.node_count; ++n)
        if (std::fabs(Pl.mesh.y[n] - L) < 1e-9) { top += -Pl.disp[2 * n + 1]; ++ntop; }
    top = ntop ? top / ntop : 0.0;
    const double want = gw / E * ((L - hd) * hd + 0.5 * (L - hd) * (L - hd));
    std::printf("  top settles %.6f m, closed form %.6f m (%+.3f%%)\n", top, want, 100.0 * (top - want) / want);
    check(std::fabs(top - want) < 0.01 * want, "the top settles by the effective stress the drawdown added (1%)");

    // A Safety phase after it reads the same field: the drawdown, not a settled water table.
    {
        m::Project ps = pr;
        ps.materials[0].model = m::SoilModel::MohrCoulomb;   // a strength to reduce
        ps.materials[0].c = 20.0; ps.materials[0].phi = 30.0;
        m::Phase sf; sf.name = "Safety"; sf.type = m::PhaseType::Safety;
        sf.water_override = true; sf.wx = {0.0, W}; sf.wy = {hd, hd};
        sf.load_steps = 5;   // the field is set before the search; a few steps are enough to read it
        ps.phases.push_back(sf);
        const auto rs = katai::app::solve_phases(ps, M.mesh, InitialPhase::K0Procedure);
        if (!rs.empty()) std::printf("  (the Safety phase: %s)\n", rs.back().message.c_str());
        double gap = 1e300;
        if (rs.size() == 4 && rs[3].pore.size() == T.pore.size()) {
            gap = 0.0;
            for (size_t n = 0; n < T.pore.size(); ++n) gap = std::fmax(gap, std::fabs(rs[3].pore[n] - T.pore[n]));
        }
        std::printf("  the Safety phase after it: pore pressures %.3e kPa from the drawdown's\n", gap);
        check(gap <= 1e-9 * scale, "a Safety phase after it reads the drawdown's pore pressures too");
    }
    // But a STEADY flow is not handed on: two staged phases whose water is lowered from one to the
    // next each solve their own (an excavation lowering the water stage by stage). Taking the
    // first one's head left the second at the old level.
    {
        m::Project st = pr;
        st.phases.clear();
        for (const double level : {8.0, 6.0}) {
            m::Phase p; p.name = "Water " + std::to_string((int)level); p.type = m::PhaseType::Plastic;
            p.pore_from_flow = true;
            p.water_override = true; p.wx = {0.0, W}; p.wy = {level, level};
            st.phases.push_back(p);
        }
        const auto r2 = katai::app::solve_phases(st, M.mesh, InitialPhase::K0Procedure);
        check(r2.size() == 3 && r2[2].ok, "two staged phases lowering the water ran");
        if (r2.size() == 3 && r2[2].ok) {
            double base = 0.0; int nb = 0;
            for (int n = 0; n < r2[2].mesh.node_count; ++n)
                if (std::fabs(r2[2].mesh.y[n]) < 1e-9) { base += r2[2].pore[n]; ++nb; }
            base = nb ? base / nb : 0.0;
            std::printf("  second stage: base pore pressure %.4f kPa (its own water: %.4f; the first stage's: %.4f)\n",
                        base, gw * 6.0, gw * 8.0);
            check(std::fabs(base - gw * 6.0) < 1e-3 * gw * 6.0, "the second stage solves its own water, not the first's");
        }
    }

    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all drawdown hand-over checks passed\n");
    return 0;
}
