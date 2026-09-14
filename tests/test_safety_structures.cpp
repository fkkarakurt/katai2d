// A Safety run with the structures in it, checked against the one problem where the answer is
// known exactly: a stiff block on an interface.
//
// Until 2026-09 the strength-reduction search was handed no structural element at all, and every
// active one was refused (K2D-G016). The search now solves them with the soil, and three things
// have to be true of how it does so, each of which has a closed form on this block:
//
//   (a) an INTERFACE's strength is a soil strength, so it is reduced with the soil's -- c_i / SRF
//       and tan(phi_i) / SRF. The block slides when the reduced joint capacity meets the horizontal
//       load, so FoS = (B c_w + W tan(phi_w)) / H. An interface left unreduced never slides; one
//       whose wished-in-place normal stress survived into the unstressed search carries twice the
//       normal force.
//   (b) a PLATE carries its weight into the search. A plate of weight w on the block's top adds
//       w B to the normal force: FoS = (B c_w + (W + w B) tan(phi_w)) / H. A search without the
//       weight returns the answer of (a). The same w as a distributed load on the same line must
//       give the same factor, which separates "the weight is in the search" from "something that
//       happens to weigh the same is".
//   (c) an ANCHOR's capacity is the structure's own, NOT a soil strength, so it is not reduced. An
//       elastoplastic anchor pulling back against H yields before the block can slide, so
//       FoS = (B c_w + W tan(phi_w)) / (H - F). Had its capacity been reduced as well, the answer
//       would be (B c_w + W tan(phi_w) + F) / H; without the anchor, that of (a). The same F as a
//       constant point force at the same node must give the same factor.
//
// verify: KV-STR-010
//   oracle:   closed_form
//   source:   Griffiths, D.V. & Lane, P.A. (1999). Slope stability analysis by finite elements. Geotechnique 49(3), 387-403 -- the strength reduction factor of the finite-element strength-reduction method, c_f = c / SRF and tan(phi_f) = tan(phi) / SRF; combined with Coulomb friction with adhesion on a planar joint, the same statics as KV-STR-002 (the horizontal force that makes a rigid block slide is the joint's adhesion over the contact width plus the normal force times the tangent of the joint friction angle); KATAI 2D input contract (docs/k2d-format.md: iface_material / Rinter for the interface strength, plates w, anchors Fmax_tens / Lspacing / prestress)
//   locator:  the checked-in KV-STR-002 block (B = 4 m wide, 1 m high, gamma = 25 kN/m3 so W = 100 kN/m; joint c_w = 2.5 kN/m2, phi_w = 26.6 deg, R_inter = 1) with its imposed slip replaced by a horizontal distributed load on its left face from y = 0.25 to 1 m (H = 40 kN/m unless stated) and run as a Safety procedure. Stated in full: (a) FoS = (B c_w + W tan(phi_w)) / H; (b) FoS = (B c_w + (W + w B) tan(phi_w)) / H with w = 10 kN/m/m; (c) FoS = (B c_w + W tan(phi_w)) / (H - F) with F = 10 kN/m, the anchor horizontal from the block's top-left corner (0, 1) to a fixed point at (-2, 1); an anchor that attaches off its drawn line pulls along the line from the node it attached to, with F cos(theta) against H and F sin(theta) on the normal force, and the closed form is evaluated on that direction
//   quantity: factor of safety of the Safety procedure [-]
//   expected: (a) 1.501907 at H = 40, 2.002542 at H = 30, 1.251907 with c_w = 0; (b) 2.002669; (c) 2.002542, 1.716465 at an out-of-plane spacing of 2 m
//   band:     one bisection interval of the search, (3.0 - 0.4) / 2^12 = 6.35e-4 absolute, as asserted below -- the search reports the midpoint of its last bracket, and nothing else separates it from the closed form on this problem: measured -2.71e-4 / -0.77e-4 / +2.66e-4 for (a), -2.04e-4 for (b) and -0.77e-4 / -2.78e-4 for (c). Four further checks: the plate's weight and the same distributed load give the same factor bit for bit, and so do the anchor and the same point force; an elastic anchor (no capacity) holds the block through the whole search, which reports a lower bound; and the anchor drawn at mid-height attached to the node at y = 0.53125 m, pulled 0.9 deg downward, and read +0.123% against the horizontal closed form -- six bisection intervals -- while the closed form written on the direction it actually pulled along matches within one. The last one is the check that this band can see a structure's geometry, not only its presence

#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kB = 4.0, kW = 100.0, kCw = 2.5, kPhiw = 26.6;
// The search brackets [0.4, 3.0] and halves it twelve times (phase_solver/safety.hpp).
constexpr double kInterval = (3.0 - 0.4) / 4096.0;

double tanphi() { return std::tan(kPhiw * kPi / 180.0); }

bool raised(const katai::app::SolveResult& R, const char* code) {
    for (const auto& d : R.diagnostics)
        if (d.code == code) return true;
    return false;
}

// The checked-in KV-STR-002 block, turned from a push into a Safety problem: the imposed slip goes,
// a horizontal load on the left face comes in, and the initial procedure is the Safety search.
m::Project block(double H = 40.0) {
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-str-002-sliding-block.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        check(false, "load " + path + ": " + err);
        return pr;
    }
    pr.initial_procedure = m::InitialProcedure::Safety;
    pr.phases.clear();
    pr.disps.clear();
    pr.initial.disp_active.clear();
    m::Load h;
    h.kind = m::LoadKind::Distributed;
    h.name = "H";
    h.x1 = 0.0; h.y1 = 0.25; h.x2 = 0.0; h.y2 = 1.0;
    h.qx1 = h.qx2 = H / 0.75;
    h.qy1 = h.qy2 = 0.0;
    pr.loads = {h};
    pr.initial.load_active = {1};
    pr.initial.struct_active = {1};
    return pr;
}

void add_plate(m::Project& pr, double w) {
    m::PlateMaterial pm;
    pm.name = "Slab"; pm.EA = 1.0e7; pm.EI = 1.0e4; pm.nu = 0.0; pm.w = w;
    pr.plates = {pm};
    m::StructElement s;
    s.kind = m::StructKind::Plate; s.name = "Slab"; s.material = 0;
    s.x1 = 0.0; s.y1 = 1.0; s.x2 = kB; s.y2 = 1.0;
    pr.structs.push_back(s);
    pr.initial.struct_active.push_back(1);
}

void add_anchor(m::Project& pr, double y, bool elastoplastic, double F, double spacing = 1.0,
                double prestress = 0.0) {
    m::AnchorMaterial am;
    am.name = "Tie"; am.EA = 1.0e5; am.Lspacing = spacing; am.prestress = prestress;
    am.elastoplastic = elastoplastic;
    am.Fmax_tens = F;
    pr.anchors = {am};
    m::StructElement s;
    s.kind = m::StructKind::Anchor; s.name = "Tie"; s.material = 0;
    s.x1 = 0.0; s.y1 = y; s.x2 = -2.0; s.y2 = y;     // fixed end outside the block
    pr.structs.push_back(s);
    pr.initial.struct_active.push_back(1);
}

void add_load(m::Project& pr, m::LoadKind kind, double x1, double y1, double x2, double y2,
              double qx, double qy) {
    m::Load L;
    L.kind = kind; L.name = "Equivalent";
    L.x1 = x1; L.y1 = y1; L.x2 = x2; L.y2 = y2;
    L.qx1 = L.qx2 = qx;
    L.qy1 = L.qy2 = qy;
    pr.loads.push_back(L);
    pr.initial.load_active.push_back(1);
}

katai::app::SolveResult safety(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) {
        check(false, "meshed: " + M.message);
        return {};
    }
    katai::app::PhaseIO io;
    io.config = &pr.initial;   // the file's own activation flags, as a run from the file reads them
    return katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety, nullptr, io);
}

// Run, print, and hold the factor to one bisection interval of its closed form.
double against(const char* what, const m::Project& pr, double exact) {
    const auto R = safety(pr);
    if (!R.ok) std::printf("      (%s)\n", R.message.c_str());
    std::printf("      %-44s FoS %.10f  closed form %.6f  (%+.2e)\n", what, R.fos, exact,
                R.fos - exact);
    check(R.ok && !R.fos_lower_bound && std::fabs(R.fos - exact) <= kInterval,
          std::string(what) + ": within one bisection interval of the closed form");
    check(raised(R, "K2D-A018"), std::string(what) + ": the run says the structures took part");
    return R.fos;
}

void interface_is_reduced_with_the_soil() {
    std::printf("\n(a) the interface's strength is reduced with the soil's\n");
    const double cap = kB * kCw + kW * tanphi();
    against("H = 40 kN/m", block(40.0), cap / 40.0);
    against("H = 30 kN/m", block(30.0), cap / 30.0);
    // The friction term alone: the adhesion is not what carries the answer.
    m::Project no_c = block(40.0);
    no_c.materials[1].c = 0.0;
    against("c_w = 0, friction alone", no_c, kW * tanphi() / 40.0);
}

void plate_weight_enters_the_search() {
    std::printf("\n(b) a plate carries its weight into the search\n");
    constexpr double w = 10.0;
    const double exact = (kB * kCw + (kW + w * kB) * tanphi()) / 40.0;
    m::Project slab = block();
    add_plate(slab, w);
    const double f_plate = against("plate w = 10 kN/m/m on the top", slab, exact);
    // Identity: the plate weightless, and its weight as a distributed load on the same line. The
    // stiffness is the same in both, so only where the load vector came from differs.
    m::Project loaded = block();
    add_plate(loaded, 0.0);
    add_load(loaded, m::LoadKind::Distributed, 0.0, 1.0, kB, 1.0, 0.0, -w);
    const double f_load = against("weightless plate + 10 kN/m/m load", loaded, exact);
    std::printf("      plate weight vs the same load: %s\n",
                f_plate == f_load ? "bit for bit" : "different");
    check(std::fabs(f_plate - f_load) <= kInterval,
          "the plate's weight and the same distributed load give the same factor");
    // What the band would have seen had the weight been left out: the answer of (a).
    const double without = (kB * kCw + kW * tanphi()) / 40.0;
    check(std::fabs(f_plate - without) > 100.0 * kInterval,
          "a search without the plate's weight is far outside the band (it is not a close call)");
}

void anchor_capacity_is_not_reduced() {
    std::printf("\n(c) an anchor's capacity is the structure's own, not reduced\n");
    constexpr double F = 10.0;
    const double cap = kB * kCw + kW * tanphi();
    m::Project tied = block();
    add_anchor(tied, 1.0, true, F);
    const double f_anchor = against("elastoplastic anchor, F = 10 kN", tied, cap / (40.0 - F));
    // Identity: the anchor's capacity as a constant point force at the node it attaches to.
    m::Project pushed = block();
    add_load(pushed, m::LoadKind::Point, 0.0, 1.0, 0.0, 1.0, -F, 0.0);
    const auto RP = safety(pushed);
    std::printf("      constant point force -10 kN at (0, 1): FoS %.10f (%s the anchor)\n", RP.fos,
                RP.fos == f_anchor ? "bit for bit" : "differs from");
    check(RP.ok && std::fabs(RP.fos - f_anchor) <= kInterval,
          "a yielded anchor and a constant force equal to its capacity give the same factor");
    // The two answers the band has to be able to reject.
    const double reduced = (cap + F) / 40.0, absent = cap / 40.0;
    std::printf("      had the capacity been reduced: %.6f; without the anchor: %.6f\n", reduced,
                absent);
    check(std::fabs(f_anchor - reduced) > 100.0 * kInterval &&
              std::fabs(f_anchor - absent) > 100.0 * kInterval,
          "a reduced capacity, or no anchor, would be far outside the band");
    // The capacity is per metre of wall: two metres between anchors halve it.
    m::Project spaced = block();
    add_anchor(spaced, 1.0, true, F, 2.0);
    against("F = 10 kN at 2 m spacing", spaced, cap / (40.0 - 0.5 * F));

    // An elastic anchor has no capacity, so the block cannot slide at any strength the search
    // reaches: the factor is a lower bound, reported as one.
    m::Project stiff = block();
    add_anchor(stiff, 1.0, false, 0.0);
    const auto RS = safety(stiff);
    std::printf("      elastic anchor: FoS %.6f, lower bound %d\n", RS.fos, (int)RS.fos_lower_bound);
    check(RS.ok && RS.fos_lower_bound, "an elastic anchor holds the block through the whole search");

    // THE BAND CAN SEE GEOMETRY. Drawn at mid-height the anchor attaches to the nearest node, which
    // on this mesh is not on its line, so it pulls slightly downward: that adds F sin(theta) to the
    // normal force and takes F (1 - cos(theta)) off its resistance. The horizontal closed form is
    // then several intervals wrong, and the one written on the actual direction is not.
    m::Project tilted = block();
    add_anchor(tilted, 0.5, true, F);
    const auto RT = safety(tilted);
    check(RT.ok, "the anchor drawn at mid-height runs");
    if (!RT.ok) return;
    int node = -1;
    double best = 1e300;
    for (int n = 0; n < RT.mesh.node_count; ++n) {
        const double d = std::hypot(RT.mesh.x[n] - 0.0, RT.mesh.y[n] - 0.5);
        if (d < best) { best = d; node = n; }
    }
    const double dx = -2.0 - RT.mesh.x[node], dy = 0.5 - RT.mesh.y[node];
    const double L = std::hypot(dx, dy);
    const double along_H = F * (-dx / L), on_normal = F * (-dy / L);
    const double tilted_exact = (kB * kCw + (kW + on_normal) * tanphi()) / (40.0 - along_H);
    const double level_exact = cap / (40.0 - F);
    std::printf("      attached at (%.5f, %.5f), %.2f deg below its line: FoS %.10f; "
                "horizontal closed form %.6f (%+.3f%%), on the actual direction %.6f (%+.3f%%)\n",
                RT.mesh.x[node], RT.mesh.y[node], std::atan2(-dy, -dx) * 180.0 / kPi, RT.fos,
                level_exact, 100.0 * (RT.fos - level_exact) / level_exact, tilted_exact,
                100.0 * (RT.fos - tilted_exact) / tilted_exact);
    check(std::fabs(RT.fos - tilted_exact) <= kInterval,
          "the closed form on the direction the anchor actually pulls along is met");
    check(std::fabs(RT.fos - level_exact) > 3.0 * kInterval,
          "and the horizontal one is not -- the band sees a 0.9 degree tilt");
}

void prestressed_anchor_is_refused() {
    std::printf("\n(d) a prestressed anchor cannot enter a search from the unstressed state\n");
    m::Project pre = block();
    add_anchor(pre, 1.0, true, 10.0, 1.0, 5.0);
    const auto R = safety(pre);
    check(!R.ok && R.fos < 0.0 && raised(R, "K2D-G016"),
          "the engine refuses it (K2D-G016), with no factor of safety");
    bool schema = false;
    for (const auto& is : katai::io::validate_project(pre).issues)
        schema |= is.severity == katai::io::Severity::Error && is.path == "initial.struct";
    check(schema, "and the .k2d contract refuses it at initial.struct");
    // The remedy the message names runs.
    pre.initial.struct_active = {1, 0};
    const auto RO = safety(pre);
    check(RO.ok && !raised(RO, "K2D-G016"), "deactivated in the Safety run, it runs");
}

}  // namespace

int main() {
    std::printf("KV-STR-010: a Safety run with the structures in it, on the sliding block\n");
    interface_is_reduced_with_the_soil();
    plate_weight_enters_the_search();
    anchor_capacity_is_not_reduced();
    prestressed_anchor_is_refused();
    if (g_failures == 0) {
        std::printf("\nOK: interface reduced, plate weight and anchor capacity carried unreduced\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
