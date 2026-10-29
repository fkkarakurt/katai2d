// A distributed load drawn along a plate that has an interface is carried by the plate.
//
// A plate with an interface is a wall: the mesh is split along its line, and the plate meets the
// ground on that side only through the interface. A point load drawn on such a plate was already
// taken by the plate; a distributed load was assembled on the mesh nodes of its line -- the ground
// on one side -- so a raft with an interface beneath it, loaded along its line, carried nothing:
// no moment, no contact pressure, the load went into the ground directly (until 2026-10-26).
//
// THE ORACLE IS THE RAFT'S EQUILIBRIUM. A weightless raft on the ground carries the load q drawn
// along it to the ground through its interface, and nothing else touches it: the contact pressure
// integrated over the raft equals the load, q B. With the interface on the ground side the
// contact is the interface's normal stress; with it on the other side, or none, the raft is bonded
// to the ground and its bending is the measure -- the same raft, the same moment within the
// interface's compliance.
//
// BOTH ENDS. A wall hangs from its deeper end, the node it shares with the ground; a raft has
// none, and the toe rule took its left end and welded it to the ground there. Measured before
// 2026-10-26: the contact carried -777.7 of the -800 kN/m (-2.8%), the left end's share passing
// through the weld. A level wall with ground on one side only is split at both ends now.
//
// verify: KV-STR-026
//   oracle:   closed_form
//   source:   vertical equilibrium of a weightless raft: the resultant of the contact pressure equals the applied load (statics); KATAI 2D input contract (docs/k2d-format.md, structs[i].iface_pos / iface_neg, loads[i])
//   locator:  a linear elastic block 20 x 10 m (E = 50000 kPa, nu = 0.3, weightless); a raft (EA = 1e7 kN/m, EI = 1e4 kNm2/m) on its surface from x = 6 to 14 m; a distributed load q = -100 kPa drawn along the raft; the interface beneath (the positive side of a line drawn left to right), above, or none
//   quantity: the normal stress of the interface beneath the raft integrated along it [kN/m]; the raft's largest moment [kNm/m] with the interface beneath, against the bonded raft's
//   expected: the contact force q B = -800 kN/m; the raft bends, as the bonded one does
//   band:     0.1% on the contact force (the joint's own Newton-Cotes rule over each element's three stations; measured 0.000%; a trapezoid over the shared points gives +1.7%, the contact peaking at the raft's ends); the moment within 25% of the bonded raft's (measured -2.0%: the interface's compliance redistributes the contact), and above half of it -- the raft that carried nothing had 0
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

struct Out { bool ok = false, simpson_ok = true; double m_max = 0.0, resultant = 0.0, trapezoid = 0.0; int stations = 0; };

// side: 0 no interface, 1 beneath (positive), 2 above (negative)
Out raft(int side) {
    m::Project pr;
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
    s.E = 50000.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon g;
    g.name = "Ground"; g.material = 0; g.x = {0, 20, 20, 0}; g.y = {0, 0, 10, 10};
    g.edge_bc.assign(4, 0);
    pr.polygons.push_back(g);
    m::PlateMaterial pl;
    pl.name = "Raft"; pl.EA = 1.0e7; pl.EI = 1.0e4;
    pr.plates.push_back(pl);
    m::StructElement e;
    e.kind = m::StructKind::Plate; e.name = "Raft";
    e.x1 = 6; e.y1 = 10; e.x2 = 14; e.y2 = 10;
    e.material = 0;
    e.iface_pos = side == 1; e.iface_neg = side == 2;
    pr.structs.push_back(e);
    m::Load L;
    L.kind = m::LoadKind::Distributed; L.name = "q";
    L.x1 = 6; L.y1 = 10; L.x2 = 14; L.y2 = 10;
    L.qx1 = L.qx2 = 0.0; L.qy1 = L.qy2 = -100.0;
    pr.loads.push_back(L);
    pr.initial.load_active = {1};
    Out o;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return o;
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res[0].ok) return o;
    o.ok = true;
    for (const auto& f : res[0].struct_forces) o.m_max = std::max(o.m_max, f.max_M);
    // The contact force: sigma_n integrated along the raft (trapezoid over the stations in order
    // of x; a point two elements share is reported by both and counted once).
    std::vector<std::pair<double, double>> pts;
    for (const auto& f : res[0].interface_forces)
        for (const auto& st : f.stations) { pts.push_back({st.x, st.sigma_n}); ++o.stations; }
    std::sort(pts.begin(), pts.end());
    std::vector<std::pair<double, double>> uniq;
    for (const auto& p : pts)
        if (uniq.empty() || p.first - uniq.back().first > 1e-9) uniq.push_back(p);
    for (size_t i = 1; i < uniq.size(); ++i)
        o.trapezoid += 0.5 * (uniq[i].second + uniq[i - 1].second) * (uniq[i].first - uniq[i - 1].first);
    // The joint's own rule: each element's three stations (its nodes, in order) by Newton-Cotes
    // -- Simpson. This is the force the joint passes, which the raft's equilibrium fixes; the
    // trapezoid over the shared points is not, where the contact peaks at the raft's ends.
    for (const auto& f : res[0].interface_forces) {
        std::vector<std::pair<double, double>> el;
        for (const auto& st : f.stations) el.push_back({st.s, st.sigma_n});
        if (el.size() % 3 != 0) { o.simpson_ok = false; continue; }
        for (size_t i = 0; i + 2 < el.size(); i += 3)
            o.resultant += std::fabs(el[i + 2].first - el[i].first) / 6.0 *
                           (el[i].second + 4.0 * el[i + 1].second + el[i + 2].second);
    }
    return o;
}

}  // namespace

int main() {
    std::printf("== a load drawn on a plate with an interface (KV-STR-026) ==\n");
    const Out bonded = raft(0), beneath = raft(1), above = raft(2);
    check(bonded.ok && beneath.ok && above.ok, "the three rafts solve");
    std::printf("      bonded: M %.3f kNm/m; interface beneath: M %.3f, contact force %.3f kN/m over %d stations "
                "(q B = -800; trapezoid over the points %.3f); interface above: M %.3f\n",
                bonded.m_max, beneath.m_max, beneath.resultant, beneath.stations, beneath.trapezoid, above.m_max);
    check(beneath.stations > 0 && beneath.simpson_ok && std::fabs(beneath.resultant / -800.0 - 1.0) < 1e-3,
          "the interface beneath the raft carries the load drawn on it: its contact force is q B, "
          "both ends included");
    check(beneath.m_max > 0.5 * bonded.m_max && std::fabs(beneath.m_max / bonded.m_max - 1.0) < 0.25,
          "...and the raft bends under it, as the bonded raft does");
    check(std::fabs(above.m_max / bonded.m_max - 1.0) < 1e-6,
          "an interface on the other side leaves the bonded raft as it was");
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
