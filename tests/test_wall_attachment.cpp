// What is drawn on a wall with interfaces acts on the wall, and the ground a phase removed holds
// nothing.
//
// A wall with interfaces is a plate on DOFs of its own, joined to the soil on each side through a
// Coulomb joint; the mesh node at a point of the wall line is the soil beside it. Four defects hid in
// that construction, and an anchored sheet pile measured against an external worked example showed
// all of them at once -- anchors 1.7 / 1.3 kN where 35 / 88 were expected, and the wall's largest
// moment at its free toe:
//
//   (a) the toe's rotation was FIXED, so the free end of an embedded wall was a clamp and carried a
//       moment (-34.8 kNm/m, the wall's largest);
//   (b) an anchor, a point load, a geogrid or a plate drawn on the wall found its node by position
//       and was attached to the SOIL beside the wall, which it tensioned against itself;
//   (c) where a phase excavated one side, that side's joints stayed on, tied to nodes that the
//       phase had fixed in space, and went on pressing the wall with the earth pressure of the
//       removed ground: the excavation never unloaded the wall (the lower anchor ended in
//       compression);
//   (d) a one-sided interface was built on both sides (iface_pos / iface_neg were only or-ed).
//
// Rebuilt with embedded beams as the tiebacks' grout bodies, the same wall showed two more, in
// the pile: an anchor ending on a free beam pulled the soil beside it rather than the beam, and the
// beam's base carried TENSION -- its formulation says it bears only -- so a pulled grout body held
// the ground by one point at its tip.
//
// Every oracle below is exact. Above the ground, a wall whose soil is excavated on both sides is a
// cantilever carrying only what is drawn on it, so statics gives its moment at every station
// whatever the ground below does: M(y) = F (y_top - y) under a point force F at its head, and
// M(y) = N (y_a - y) below an anchor carrying N. A wall whose interface sits only on the side the
// phase has excavated is bonded to the other side, which is the plain plate: the two must agree to
// round-off. A joint with no plastic slip carries tau = k_s * slip, so an entered shear stiffness is
// read back exactly. A pile on linear-elastic ground pulled up is the mirror image of the same pile
// pushed down, except for its base, the one one-sided law in the model.
//
// verify: KV-STR-014
//   oracle:   closed_form
//   source:   statics of a free cantilever (the part of the wall above an excavation on both sides carries only the load drawn on it: M(y) = F (y_top - y), Q = F; below an anchor of force N, M(y) = N_x (y_a - y)); the identity of a one-sided interface on an excavated side with a plate bonded to the remaining side; the elastic branch of the Coulomb interface, tau = k_s (u_s - u_s^p) with u_s^p = 0 before any slip (docs/k2d-format.md, structs[i].iface_ks); the KATAI 2D input contract for iface_pos / iface_neg (positive = the right of the direction from (x1, y1) to (x2, y2)); the embedded-beam base bears in compression only (docs/references/embedded-beam-formulation.md sec 4), so on a linear ground a pile pulled up mirrors the same pile pushed down in everything but its base
//   locator:  linear-elastic block 20 m x 10 m (E = 30000 kPa, nu = 0.3, gamma 18 kN/m3, c 50 kPa, phi 30 deg, K0 procedure), base fixed, sides on rollers; a plate wall x = 10 m from y = 10 m down to a toe at y = 2 m (EA = 1.2e7 kN/m, EI = 1.6e5 kNm2/m) with interfaces; the soil above y = 4 m removed on both sides, then (i) a point force of 10 kN at the wall head, or (ii) a fixed-end anchor from (10, 9) to (-5, 9) with a 50 kN lock-off force; a second block with the soil removed only on the left (x < 10) down to y = 2 m, a one-sided wall against a plain plate; a free embedded beam (E = 2e7 kPa, D = 0.2 m) at the end of a node-to-node anchor with a 50 kN lock-off force; a hinged vertical pile 6 m long (E = 3e7 kPa, D = 0.6 m, spacing 2 m) pushed down and pulled up by 100 kN at its head
//   quantity: the wall's bending moment [kNm/m] at every station above the excavated level, its moment at the free toe, the one-sided wall's displacements against the plain plate's, the joint's tau / slip ratio [kN/m3], the grout body's axial force against the anchor's [kN], and the pile's axial force at its tip pulled against pushed [kN]
//   expected: M equal to the statics above the excavated level; |M_toe| a small fraction of the peak; the one-sided wall equal to the plain plate; tau / slip = k_s as entered; the grout body carrying the anchor's force; the pulled pile's tip holding nothing
//   band:     1e-6 relative on the cantilever moments (the plate interpolates a linear moment exactly; the run converges to 1e-10); 1e-9 on the identities; 1e-9 on tau / slip. The toe, the grout body and the pile tip are discretised, and their bands are fractions stated beside the assertions, each far from the defect: the free toe measured 1.56% of the peak moment (the clamp 94%), the grout body 99.6% of the anchor's force (10.6% when the anchor pulled the soil), the pulled tip 0.81% of the pushed one (exactly 100% when the base held in tension)

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;
using katai::app::SolveResult;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

m::Material soil() {
    m::Material s;
    s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.e_init = 0.6; s.c = 50.0; s.phi = 30.0;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.rinter_rigid = false; s.Rinter = 1.0;
    s.k0_auto = false; s.k0 = 0.5;
    return s;
}

m::SoilPolygon rect(double x0, double y0, double x1, double y1, std::vector<int> bc) {
    m::SoilPolygon P; P.material = 0;
    P.x = {x0, x1, x1, x0}; P.y = {y0, y0, y1, y1};
    P.edge_bc = std::move(bc);
    return P;
}

constexpr int kFix = (int)m::BCType::FullyFixed, kH = (int)m::BCType::HorizontallyFixed,
              kFree = (int)m::BCType::Free;

// The block: 0 base (y 0..y_cut), 1 left top (x 0..10), 2 right top (x 10..20). The wall is drawn
// from its head (10, 10) down to its toe (10, 2), so its POSITIVE side is the left (x < 10).
m::Project block(double y_cut) {
    m::Project pr;
    pr.name = "wall attachment";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.mesh.elem_size = 1.0; pr.mesh.order = 6; pr.mesh.auto_refine = false;
    pr.materials.push_back(soil());
    pr.polygons.push_back(rect(0, 0, 20, y_cut, {kFix, kH, kFree, kH}));
    pr.polygons.push_back(rect(0, y_cut, 10, 10, {kFree, kFree, kFree, kH}));
    pr.polygons.push_back(rect(10, y_cut, 20, 10, {kFree, kH, kFree, kFree}));
    m::PlateMaterial pm; pm.name = "Wall";
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = 0.0;
    pr.plates.push_back(pm);
    m::StructElement wall; wall.kind = m::StructKind::Plate; wall.name = "Wall"; wall.material = 0;
    wall.x1 = 10; wall.y1 = 10; wall.x2 = 10; wall.y2 = 2;
    wall.iface_pos = wall.iface_neg = true;
    pr.structs.push_back(wall);
    return pr;
}

m::Phase phase(const char* name, std::vector<char> poly, std::vector<char> st, std::vector<char> ld) {
    m::Phase ph; ph.name = name; ph.type = m::PhaseType::Plastic;
    ph.poly_active = std::move(poly); ph.struct_active = std::move(st); ph.load_active = std::move(ld);
    return ph;
}

std::vector<SolveResult> solve(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { check(false, "meshed: " + M.message); return {}; }
    auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    for (const auto& r : res)
        if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    return res;
}

bool all_ok(const std::vector<SolveResult>& R, size_t n) {
    return R.size() == n && std::all_of(R.begin(), R.end(), [](const SolveResult& r) { return r.ok; });
}

const katai::core::StructForce* force_of(const SolveResult& r, const char* name) {
    for (const auto& f : r.struct_forces)
        if (f.name == name) return &f;
    return nullptr;
}

double peak_M(const katai::core::StructForce& f) {
    double p = 0.0;
    for (const auto& s : f.stations) p = std::max(p, std::fabs(s.M));
    return p;
}

// Largest |M(y)| - M_statics(y) over the stations of the elements lying wholly in (y_lo, y_hi],
// relative to the statics' peak. Element by element because the plate's moment is piecewise
// linear: an element that straddles the ground level carries interface tractions over part of it,
// and one with a point force at its MIDDLE node cannot represent the kink there. Both are the
// discretisation, and the statics holds exactly on every other element. The diagram lists each
// element's stations bottom to top, so an element boundary is a height listed twice.
template <class Fn>
double cantilever_gap(const katai::core::StructForce& f, double y_lo, double y_hi, Fn statics,
                      int& n) {
    std::vector<double> ends;   // element boundaries: heights that two elements share
    for (size_t k = 1; k < f.stations.size(); ++k)
        if (std::fabs(f.stations[k].y - f.stations[k - 1].y) < 1e-9) ends.push_back(f.stations[k].y);
    double hi = y_lo;           // the highest element boundary at or below y_hi
    for (double e : ends) if (e <= y_hi + 1e-9) hi = std::max(hi, e);
    if (y_hi >= 10.0 - 1e-9) hi = y_hi;   // up to the head: the last element ends there
    double gap = 0.0, top = 0.0;
    n = 0;
    bool at_hi = false;
    for (const auto& s : f.stations) {
        if (s.y > hi + 1e-9 || (at_hi && std::fabs(s.y - hi) < 1e-9)) break;   // the next element
        if (std::fabs(s.y - hi) < 1e-9) at_hi = true;
        if (s.y < y_lo + 1e-6) continue;
        gap = std::max(gap, std::fabs(std::fabs(s.M) - statics(s.y)));
        top = std::max(top, statics(s.y));
        ++n;
    }
    return top > 0.0 ? gap / top : 1e300;
}

// (i) A point force at the head of a wall whose ground is gone above y = 4: M = F (10 - y).
// (a) is read at the toe of the same run.
void point_force_on_a_cantilever() {
    std::printf("\n(i) a point force at the head of the wall, the ground removed above y = 4\n");
    m::Project pr = block(4.0);
    constexpr double F = 10.0;
    m::Load pt; pt.kind = m::LoadKind::Point; pt.name = "head";
    pt.x1 = pt.x2 = 10.0; pt.y1 = pt.y2 = 10.0; pt.qx1 = pt.qx2 = -F; pt.qy1 = pt.qy2 = 0.0;
    pr.loads.push_back(pt);
    pr.initial.struct_active = {0};
    pr.initial.load_active = {0};
    pr.phases = {phase("wall", {1, 1, 1}, {1}, {0}), phase("excavate", {1, 0, 0}, {1}, {0}),
                 phase("load", {1, 0, 0}, {1}, {1})};
    const auto R = solve(pr);
    check(all_ok(R, 4), "the chain solves: initial, wall, excavation, load");
    if (!all_ok(R, 4)) return;
    const auto* w = force_of(R[3], "Wall");
    check(w && !w->stations.empty(), "the wall reports its forces");
    if (!w || w->stations.empty()) return;
    int n = 0;
    const double gap = cantilever_gap(*w, 4.0, 10.0, [&](double y) { return F * (10.0 - y); }, n);
    std::printf("      %d stations above the ground: max |M| %.6f kNm/m (statics %.6f), gap %.2e\n", n,
                peak_M(*w), F * 6.0, gap);
    check(n >= 6, "stations above the ground (teeth)");
    // The defect (b) put the force on a soil node the phase had fixed, and (c) held the wall on
    // the removed side: either leaves this moment far from F (10 - y).
    check(gap <= 1e-6, "(b)(c) the free part of the wall carries exactly the statics, F (10 - y)");

    // (a) The toe is a free end. Its M = 0 is a natural boundary condition, met by the
    // discretisation rather than imposed, so the band is a fraction of the peak: measured 1.56%,
    // where the clamped toe of the defect carried 94% of it (57.05 of 60.71 kNm/m, this run).
    const auto toe = std::min_element(w->stations.begin(), w->stations.end(),
                                      [](const auto& a, const auto& b) { return a.y < b.y; });
    const double M_toe = std::fabs(toe->M), M_pk = peak_M(*w);
    std::printf("      toe (y = %.3f): M %.6f kNm/m, %.3f%% of the peak\n", toe->y, toe->M,
                100.0 * M_toe / M_pk);
    check(std::fabs(toe->y - 2.0) < 1e-6, "the lowest station is the toe");
    check(M_toe <= 0.02 * M_pk, "(a) the free toe carries no moment");
}

// (ii) An anchor from the wall: below it, M = N_x (9 - y) with N the anchor's own force.
void anchor_on_a_cantilever() {
    std::printf("\n(ii) a prestressed anchor on the wall, the ground removed above y = 4\n");
    m::Project pr = block(4.0);
    m::AnchorMaterial am; am.name = "Tie"; am.EA = 2.0e5; am.Lspacing = 1.0; am.prestress = 50.0;
    pr.anchors.push_back(am);
    m::StructElement tie; tie.kind = m::StructKind::Anchor; tie.name = "Tie"; tie.material = 0;
    tie.x1 = 10; tie.y1 = 9; tie.x2 = -5; tie.y2 = 9;   // horizontal, fixed far end off the soil
    pr.structs.push_back(tie);
    pr.initial.struct_active = {0, 0};
    pr.phases = {phase("wall", {1, 1, 1}, {1, 0}, {}), phase("excavate", {1, 0, 0}, {1, 0}, {}),
                 phase("anchor", {1, 0, 0}, {1, 1}, {})};
    const auto R = solve(pr);
    check(all_ok(R, 4), "the chain solves: initial, wall, excavation, anchor");
    if (!all_ok(R, 4)) return;
    const auto* w = force_of(R[3], "Wall");
    const auto* a = force_of(R[3], "Tie");
    check(w && a && !w->stations.empty() && !a->stations.empty(), "wall and anchor report forces");
    if (!w || !a || w->stations.empty() || a->stations.empty()) return;
    // The anchor acts at the node it attaches to -- the wall node nearest its drawn head -- and
    // pulls along the line from there to its fixed point; statics is written from that node.
    const double N = a->stations[0].N, xa = a->stations[0].x, ya = a->stations[0].y;
    const double Nx = N * std::fabs(-5.0 - xa) / std::hypot(-5.0 - xa, 9.0 - ya);
    int n = 0;
    const double gap =
        cantilever_gap(*w, 4.0, ya, [&](double y) { return Nx * (ya - y); }, n);
    std::printf("      anchor N %.6f kN at (%.3f, %.3f); wall max |M| %.6f kNm/m (statics %.6f), "
                "gap %.2e\n", N, xa, ya, peak_M(*w), Nx * (ya - 4.0), gap);
    // The lock-off force relaxes as the cantilever deflects towards the anchor (EA/L = 13 MN/m
    // against a 5 m cantilever on elastic ground); what matters is that the WALL resists it. An
    // anchor attached to the soil beside the wall held its force against a node the phase had
    // fixed, and the wall carried nothing: that is what (b) did.
    check(std::fabs(xa - 10.0) < 1e-9, "the anchor attaches on the wall line");
    check(N > 1.0, "the anchor is in tension (teeth)");
    check(n >= 4, "stations between the ground and the anchor (teeth)");
    check(gap <= 1e-6, "(b) below the anchor the wall carries exactly N_x (y_a - y)");
}

// (iii) One-sided interfaces, with the ground removed on the left down to the toe. The wall is
// drawn head to toe, so its positive side is the left. An interface ONLY on the left sits against
// nothing: the wall is bonded to the right-hand soil, which is the plain plate. An interface on
// both sides is the right-hand one alone, because the left one has no ground (c). And an entered
// shear stiffness is what the joint carries (tau = k_s slip before any slip). The head is propped
// by a fixed-end anchor: a wall with a free toe, nothing in front and a joint that cannot pull
// behind is a mechanism -- the retained ground pushes it over -- and the solver says so.
void one_sided_interfaces() {
    std::printf("\n(iii) one-sided interfaces, the ground removed on the left down to y = 2\n");
    const auto run = [](bool pos, bool neg, double ks) {
        m::Project pr = block(2.0);
        pr.structs[0].iface_pos = pos; pr.structs[0].iface_neg = neg;
        pr.structs[0].iface_ks = ks;
        m::Load pt; pt.kind = m::LoadKind::Point; pt.name = "head";
        pt.x1 = pt.x2 = 10.0; pt.y1 = pt.y2 = 10.0; pt.qx1 = pt.qx2 = 20.0; pt.qy1 = pt.qy2 = 0.0;
        pr.loads.push_back(pt);
        m::AnchorMaterial am; am.name = "Prop"; am.EA = 2.0e5; am.Lspacing = 1.0;
        pr.anchors.push_back(am);
        m::StructElement prop; prop.kind = m::StructKind::Anchor; prop.name = "Prop"; prop.material = 0;
        prop.x1 = 10; prop.y1 = 10; prop.x2 = -5; prop.y2 = 10;
        pr.structs.push_back(prop);
        pr.initial.struct_active = {0, 0};
        pr.initial.load_active = {0};
        pr.phases = {phase("wall", {1, 1, 1}, {1, 1}, {0}), phase("excavate", {1, 0, 1}, {1, 1}, {0}),
                     phase("load", {1, 0, 1}, {1, 1}, {1})};
        return solve(pr);
    };
    std::printf("      plain:\n");
    const auto plain = run(false, false, 0.0);
    std::printf("      positive only:\n");
    const auto left = run(true, false, 0.0);    // positive = left, the excavated side
    std::printf("      negative only:\n");
    const auto right = run(false, true, 0.0);   // negative = right, the retained side
    std::printf("      both:\n");
    const auto both = run(true, true, 0.0);
    check(all_ok(plain, 4) && all_ok(left, 4) && all_ok(right, 4) && all_ok(both, 4),
          "plain, positive-only, negative-only and both solve");
    if (!all_ok(plain, 4) || !all_ok(left, 4) || !all_ok(right, 4) || !all_ok(both, 4)) return;
    // Compared on the WALL: its moment and its own displacement, station by station. (The field's
    // max|u| is not the same measure in every variant: where the bonded side is the left, the wall
    // moves on the left twins' DOFs and they report its movement; with two joints they are
    // orphaned and fixed.)
    const auto gap = [](const SolveResult& a, const SolveResult& b) {
        const auto* wa = force_of(a, "Wall");
        const auto* wb = force_of(b, "Wall");
        if (!wa || !wb || wa->stations.size() != wb->stations.size() || wa->stations.empty()) return 1e300;
        // A plain plate lists its stations in the drawn direction (head to toe), a wall with
        // interfaces from its toe up: bring both to bottom-up before pairing them.
        auto sa = wa->stations, sb = wb->stations;
        if (sa.front().y > sa.back().y) std::reverse(sa.begin(), sa.end());
        if (sb.front().y > sb.back().y) std::reverse(sb.begin(), sb.end());
        double dM = 0.0, dU = 0.0, M = 0.0, U = 0.0;
        for (size_t k = 0; k < sa.size(); ++k) {
            const auto &p = sa[k], &q = sb[k];
            if (std::fabs(p.y - q.y) > 1e-9) return 1e300;
            dM = std::max(dM, std::fabs(std::fabs(p.M) - std::fabs(q.M)));   // M's sign follows the element direction
            dU = std::max(dU, std::hypot(p.ux - q.ux, p.uy - q.uy));
            M = std::max({M, std::fabs(p.M), std::fabs(q.M)});
            U = std::max({U, std::hypot(p.ux, p.uy), std::hypot(q.ux, q.uy)});
        }
        return std::max(dM / M, dU / U);
    };
    std::printf("      wall max|M| plain %.9f | positive %.9f | negative %.9f | both %.9f\n",
                force_of(plain[3], "Wall") ? peak_M(*force_of(plain[3], "Wall")) : -1.0,
                force_of(left[3], "Wall") ? peak_M(*force_of(left[3], "Wall")) : -1.0,
                force_of(right[3], "Wall") ? peak_M(*force_of(right[3], "Wall")) : -1.0,
                force_of(both[3], "Wall") ? peak_M(*force_of(both[3], "Wall")) : -1.0);
    std::printf("      gaps: positive-plain %.2e | negative-plain %.2e | both-negative %.2e\n",
                gap(left[3], plain[3]), gap(right[3], plain[3]), gap(both[3], right[3]));
    check(gap(left[3], plain[3]) <= 1e-9,
          "(d) an interface only on the excavated (positive) side is the plain bonded plate");
    check(gap(right[3], plain[3]) > 1e-3,
          "(d) an interface on the retained (negative) side is not (teeth: the flag is read)");
    check(gap(both[3], right[3]) <= 1e-9,
          "(c) with the left ground removed, both sides is the right-hand joint alone");

    // An entered shear stiffness, read back from the joint's own stations.
    constexpr double ks = 1234.5;
    const auto own = run(false, true, ks);
    check(all_ok(own, 4), "the same wall with k_s entered solves");
    if (!all_ok(own, 4)) return;
    double worst = 0.0; int n = 0;
    for (const auto& it : own[3].interface_forces)
        for (const auto& st : it.stations) {
            if (st.slipping || std::fabs(st.slip) < 1e-7) continue;
            worst = std::max(worst, std::fabs(st.tau / st.slip - ks) / ks);
            ++n;
        }
    std::printf("      k_s entered %.1f: %d elastic stations, worst |tau/slip - k_s| / k_s %.2e\n", ks, n,
                worst);
    check(n >= 4, "stations with slip to read the stiffness from (teeth)");
    check(worst <= 1e-9, "the joint's shear stiffness is the one entered");
}

// (iv) An anchor ending on a FREE embedded beam pulls the beam: the grout body of a ground
// anchor. The beam's largest axial force is the anchor's less the skin traction over the head's
// share of the first element: 99.6% of it measured, against 10.6% when the anchor pulled the soil
// beside the beam (the beam then took only what the ground passed to it through its skin). The
// beam is horizontal, drawn from its head: its base is the FAR end, which the pull lifts off --
// a base bears and does not hold (84.9% when this beam's base sat at its head and held).
void anchor_on_a_free_beam() {
    std::printf("\n(iv) an anchor ending on a free embedded beam\n");
    m::Project pr;
    pr.name = "grout body";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.mesh.elem_size = 1.0; pr.mesh.order = 6; pr.mesh.auto_refine = false;
    pr.materials.push_back(soil());
    pr.polygons.push_back(rect(0, 0, 20, 10, {kFix, kH, kFree, kH}));
    m::AnchorMaterial am; am.name = "Tie"; am.EA = 2.0e5; am.Lspacing = 1.0; am.prestress = 50.0;
    pr.anchors.push_back(am);
    m::StructElement tie; tie.kind = m::StructKind::Anchor; tie.name = "Tie"; tie.material = 0;
    tie.x1 = 5; tie.y1 = 5; tie.x2 = 10; tie.y2 = 5;
    pr.structs.push_back(tie);
    m::EmbeddedBeamMaterial em; em.name = "Grout"; em.E = 2.0e7; em.gamma = 0.0; em.diameter = 0.2;
    em.Lspacing = 1.0;
    pr.embedded.push_back(em);
    m::StructElement g; g.kind = m::StructKind::EmbeddedBeam; g.name = "Grout"; g.material = 0;
    g.x1 = 10; g.y1 = 5; g.x2 = 16; g.y2 = 5; g.conn = 1;   // free: the grout body
    pr.structs.push_back(g);
    pr.initial.struct_active = {0, 0};
    pr.phases = {phase("tension", {1}, {1, 1}, {})};
    const auto R = solve(pr);
    check(all_ok(R, 2), "the anchor and its grout body solve");
    if (!all_ok(R, 2)) return;
    const auto* a = force_of(R[1], "Tie");
    const auto* b = force_of(R[1], "Grout");
    check(a && b && !a->stations.empty() && !b->stations.empty(), "anchor and beam report forces");
    if (!a || !b || a->stations.empty() || b->stations.empty()) return;
    const double N = a->stations[0].N;
    double Nb = 0.0;
    for (const auto& s : b->stations) Nb = std::max(Nb, std::fabs(s.N));
    std::printf("      anchor N %.4f kN, beam max |N| %.4f kN (%.2f%%)\n", N, Nb, 100.0 * Nb / N);
    check(N > 1.0, "the anchor is in tension (teeth)");
    check(Nb >= 0.8 * N && Nb <= 1.0 * N + 1e-9, "the free beam carries the anchor's force into the ground");
}

// (v) A pile base bears and does not hold. On a linear-elastic ground every spring of the pile is
// linear, so a pile pushed down and the same pile pulled up by the same force are mirror images --
// unless the base lifts off in tension, which is the only one-sided law in the model. The beam's
// axial force at its tip is the base's reaction less the skin over the tip's share of the first
// element: pushed, the base carries it; pulled, only that skin share is left.
void pile_base_in_tension() {
    std::printf("\n(v) a pile pushed down and pulled up by the same force\n");
    const auto run = [](double qy) {
        m::Project pr;
        pr.name = "pile base";
        pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
        pr.has_water = false;
        pr.mesh.elem_size = 1.0; pr.mesh.order = 6; pr.mesh.auto_refine = false;
        pr.materials.push_back(soil());
        pr.polygons.push_back(rect(0, 0, 20, 10, {kFix, kH, kFree, kH}));
        m::EmbeddedBeamMaterial em; em.name = "Pile"; em.E = 3.0e7; em.gamma = 0.0;
        em.diameter = 0.6; em.Lspacing = 2.0;
        pr.embedded.push_back(em);
        m::StructElement p; p.kind = m::StructKind::EmbeddedBeam; p.name = "Pile"; p.material = 0;
        p.x1 = 10; p.y1 = 10; p.x2 = 10; p.y2 = 4; p.conn = 0;   // hinged head at the surface
        pr.structs.push_back(p);
        m::Load pt; pt.kind = m::LoadKind::Point; pt.name = "head";
        pt.x1 = pt.x2 = 10.0; pt.y1 = pt.y2 = 10.0; pt.qx1 = pt.qx2 = 0.0; pt.qy1 = pt.qy2 = qy;
        pr.loads.push_back(pt);
        pr.initial.struct_active = {1};
        pr.initial.load_active = {0};
        pr.phases = {phase("load", {1}, {1}, {1})};
        return solve(pr);
    };
    const auto push = run(-100.0), pull = run(100.0);
    check(all_ok(push, 2) && all_ok(pull, 2), "the pushed and the pulled pile solve");
    if (!all_ok(push, 2) || !all_ok(pull, 2)) return;
    const auto tip_N = [](const SolveResult& r) {
        const auto* f = force_of(r, "Pile");
        if (!f || f->stations.empty()) return 1e300;
        return std::min_element(f->stations.begin(), f->stations.end(),
                                [](const auto& a, const auto& b) { return a.y < b.y; })->N;
    };
    const double Np = tip_N(push[1]), Nu = tip_N(pull[1]);
    std::printf("      tip N: pushed %.4f kN, pulled %.4f kN (ratio %.4f)\n", Np, Nu, std::fabs(Nu / Np));
    // A base that held would make the two mirror images: ratio 1 to round-off.
    check(std::fabs(Np) > 1.0, "the pushed pile's base bears (teeth)");
    check(std::fabs(Nu) < 0.25 * std::fabs(Np), "the pulled pile's base lifts off: it holds nothing");
}

}  // namespace

int main() {
    point_force_on_a_cantilever();
    anchor_on_a_cantilever();
    one_sided_interfaces();
    anchor_on_a_free_beam();
    pile_base_in_tension();
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures,
                g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
