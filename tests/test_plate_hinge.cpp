// A plate's hinged end: no moment through a hinge, and a pin at a support.
//
// Plates and walls that meet at a point are joined rigidly there -- one rotation for all of them --
// which is what a cast-in-place joint is. A segment joint of a lining, a strut seated on a wall, a
// slab resting on a ledge are not: they pass force and no moment. A hinged end keeps its own
// rotation, so the moment at that end is zero by its own equilibrium, whatever the ground and the
// loads -- the closed form this test holds the engine to. At a supported boundary, where a plate's
// rotation is held (the standard fixities), a hinged end is a pin: the moment there is zero too,
// where the rigid end carries the largest moment of the plate.
//
// verify: KV-STR-025
//   oracle:   closed_form
//   source:   the equilibrium of a hinge: a member end that is free to rotate carries no bending moment (any structural analysis text, e.g. Hibbeler, R.C. Structural Analysis, ch. 2, support and connection idealisations -- the pin transmits force, not moment)
//   locator:  (a) a slab on elastic ground (E = 10000 kPa, nu = 0.3, weightless; 20 x 10 m), drawn as two plates meeting at x = 10 m under a 100 kN/m point load there; (b) a 6 m pile wall standing on the fixed base of the same ground, a 50 kN/m horizontal load at its head; (c) the joint of (a) between two walls with interfaces on both sides, 3 m deep, under the same load. Each rigid, then with the end at the joint (or at the base) hinged
//   quantity: the bending moment a force station reports at the hinged end [kNm/m], against the moment the rigid joint carries there; its fall with the element size
//   expected: M = 0 at a hinged end; the rigid joint carries a moment, the same on both sides
//   band:     the station moment is recovered from the element's field, not from the end's own equation, so a hinge reads zero to the element's accuracy: below 5% of the rigid joint's at 1 m elements (measured 3.6%) and falling more than 2.5x per halving of the element size (measured 3.0x and 3.2x: 2.62, 0.88, 0.28 kNm/m) -- a moment passed through would not fall; the pinned base below 1% of the clamped one (measured 0.12%); the walls' hinge below 2% (measured 1.0%); the rigid joints' two sides alike to 0.1% (slab, measured 4e-6) and 2% (walls, measured 0.6%)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
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
}

m::Project ground(double h = 1.0) {
    m::Project pr;
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = h;
    pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
    s.E = 10000.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P;
    P.name = "Ground"; P.material = 0; P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc.assign(4, 0);   // none drawn: the standard fixities
    pr.polygons.push_back(P);
    m::PlateMaterial pl;
    pl.name = "Slab"; pl.EA = 5.0e6; pl.EI = 5.0e4; pl.nu = 0.0;
    pr.plates.push_back(pl);
    return pr;
}

m::StructElement plate(const char* name, double x1, double y1, double x2, double y2) {
    m::StructElement e;
    e.kind = m::StructKind::Plate;
    e.name = name;
    e.x1 = x1; e.y1 = y1; e.x2 = x2; e.y2 = y2;
    e.material = 0;
    return e;
}

void point_load(m::Project& pr, double x, double y, double qx, double qy) {
    m::Load L;
    L.kind = m::LoadKind::Point;
    L.name = "P";
    L.x1 = x; L.y1 = y; L.qx1 = qx; L.qy1 = qy;
    pr.loads.push_back(L);
    pr.initial.load_active.assign(pr.loads.size(), 1);
}

katai::app::SolveResult solve(const m::Project& pr, const std::string& what) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) {
        check(false, std::string(what) + ": meshes (" + M.message + ")");
        return {};
    }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    const bool ok = !res.empty() && res[0].ok;
    check(ok, std::string(what) + ": solves");
    return ok ? res[0] : katai::app::SolveResult{};
}

// The force station of structure `name` at (x, y): its moment and the rotation the plate has there.
struct End { bool found = false; double M = 0.0; };
End end_of(const katai::app::SolveResult& r, const std::string& name, double x, double y) {
    End e;
    for (const auto& f : r.struct_forces) {
        if (f.name != name) continue;
        for (const auto& st : f.stations)
            if (std::fabs(st.x - x) < 1e-6 && std::fabs(st.y - y) < 1e-6) { e.found = true; e.M = st.M; }
    }
    return e;
}

}  // namespace

int main() {
    std::printf("== a plate's hinged end (KV-STR-025) ==\n");

    // (a) A slab on the ground in two pieces, loaded at the joint.
    // The moment a station reports is recovered from the element's field, not from the end's own
    // equation, so at a hinge it is zero only to the element's accuracy: it must vanish as the
    // mesh is refined -- which is what tells a recovery residual from a moment passed through.
    std::vector<double> residual;
    for (const double h : {1.0, 0.5, 0.25}) {
        m::Project pr = ground(h);
        pr.structs.push_back(plate("Left", 2, 10, 10, 10));
        pr.structs.push_back(plate("Right", 10, 10, 18, 10));
        pr.structs[0].hinge2 = true;
        point_load(pr, 10, 10, 0, -100);
        const auto r = solve(pr, "slab, hinged joint, element size " + std::to_string(h).substr(0, 4));
        residual.push_back(std::fabs(end_of(r, "Left", 10, 10).M));
        std::printf("      element size %.2f m: M at the hinge %.6g kNm/m\n", h, residual.back());
    }
    double m_rigid = 0.0;
    for (const bool hinged : {false, true}) {
        m::Project pr = ground();
        pr.structs.push_back(plate("Left", 2, 10, 10, 10));
        pr.structs.push_back(plate("Right", 10, 10, 18, 10));
        pr.structs[0].hinge2 = hinged;
        point_load(pr, 10, 10, 0, -100);
        const auto r = solve(pr, hinged ? "slab, hinged joint" : "slab, rigid joint");
        if (r.struct_forces.empty()) continue;
        const End l = end_of(r, "Left", 10, 10), rr = end_of(r, "Right", 10, 10);
        std::printf("      %s: M at the joint, left %.6g, right %.6g kNm/m\n", hinged ? "hinged" : "rigid ", l.M, rr.M);
        if (!hinged) {
            m_rigid = std::fabs(l.M);
            check(l.found && rr.found && m_rigid > 1.0 && std::fabs(l.M - rr.M) < 1e-3 * m_rigid,
                  "a rigid joint carries its moment through, the same on both sides");
        } else {
            check(l.found && rr.found && std::fabs(l.M) < 0.05 * m_rigid && std::fabs(rr.M) < 0.05 * m_rigid,
                  "a hinged joint carries no moment, on either side (to the element's accuracy)");
        }
    }

    check(residual.size() == 3 && residual[0] > 2.5 * residual[1] && residual[1] > 2.5 * residual[2],
          "...and what remains of it vanishes with the element size (more than 2.5x per halving)");

    // (b) A pile wall on the fixed base: clamped, then pinned.
    double m_clamp = 0.0;
    for (const bool hinged : {false, true}) {
        m::Project pr = ground();
        pr.structs.push_back(plate("Wall", 10, 0, 10, 6));
        pr.structs[0].hinge1 = hinged;
        point_load(pr, 10, 6, 50, 0);
        const auto r = solve(pr, hinged ? "wall, pinned base" : "wall, clamped base");
        if (r.struct_forces.empty()) continue;
        const End b = end_of(r, "Wall", 10, 0);
        std::printf("      %s: M at the base %.6g kNm/m\n", hinged ? "pinned " : "clamped", b.M);
        if (!hinged) {
            m_clamp = std::fabs(b.M);
            check(b.found && m_clamp > 1.0, "a plate on a supported base is clamped there: it carries a moment");
        } else {
            check(b.found && std::fabs(b.M) < 0.01 * m_clamp, "a hinged end on the base is a pin: no moment");
        }
    }

    // (c) The joint of (a) between two walls with interfaces.
    double m_wall = 0.0;
    for (const bool hinged : {false, true}) {
        m::Project pr = ground();
        pr.structs.push_back(plate("Upper", 10, 9, 10, 6));
        pr.structs.push_back(plate("Lower", 10, 6, 10, 3));
        for (auto& s : pr.structs) s.iface_pos = s.iface_neg = true;
        pr.structs[0].hinge2 = hinged;
        point_load(pr, 10, 9, 40, 0);
        const auto r = solve(pr, hinged ? "walls, hinged joint" : "walls, rigid joint");
        if (r.struct_forces.empty()) continue;
        const End u = end_of(r, "Upper", 10, 6), l = end_of(r, "Lower", 10, 6);
        std::printf("      %s: M at the joint, upper %.6g, lower %.6g kNm/m\n", hinged ? "hinged" : "rigid ", u.M, l.M);
        if (!hinged) {
            m_wall = std::fabs(u.M);
            check(u.found && l.found && m_wall > 1.0 && std::fabs(std::fabs(u.M) - std::fabs(l.M)) < 0.02 * m_wall,
                  "two walls with interfaces joined rigidly carry the moment through");
        } else {
            check(u.found && l.found && std::fabs(u.M) < 0.02 * m_wall && std::fabs(l.M) < 0.02 * m_wall,
                  "...and hinged, none");
        }
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
