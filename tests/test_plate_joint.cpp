// Plates that meet at a node are joined rigidly: a plate drawn in two pieces bends as one.
//
// Found rebuilding a shield-tunnel tutorial, whose lining is an arc and so has to be drawn as a
// chain of straight plates: each drawn plate numbered its own rotations, two plates meeting at a
// node kept two independent rotations there, and every joint was a hinge -- the ring carried its
// hoop force with no bending, and nothing said so. On the case below, a plate on elastic ground
// loaded at mid-span, the moment under the load was +56.65 kNm/m in one piece and -0.23 kNm/m
// with the plate drawn as two pieces meeting under the load, the sign reversing a metre away.
//
// verify: KV-STR-015
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, structs of type plate): a plate is a beam on the mesh nodes it is drawn along, and plates that share a node are rigidly connected there (translations and rotation), so a plate drawn as two consecutive lines is the same structure as the plate drawn as one
//   locator:  a weightless elastic block 12 m x 6 m (E = 5000 kPa, nu = 0.3), base fixed, sides on rollers; a plate EA = 1e7 kN/m, EI = 1e4 kNm2/m along the surface from x = 3 to 9 m, drawn once as one line and once as two lines meeting at x = 6 m, with a point load of 100 kN/m down at (6, 0); both meshes carry a vertex at x = 6 (stated in full)
//   quantity: the plate's bending moment under the load and its largest moment, one piece against two [kNm/m]
//   expected: the two pieces give the moments of the one piece; a hinge at the joint would give a moment near zero under the load. The same with interfaces on both sides of the plate (walls on their own degrees of freedom, joined at their common end), and the force report of each wall inside its drawn ends
//   band:     1e-6 relative without interfaces, as asserted below -- the two models assemble the same mesh and, with the joint rotation shared, the same system (measured 0); 1e-3 with interfaces, where the joint's interface nodes differ between the two models (measured 3.8e-4)

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
}

m::Project model(bool two_pieces, bool interfaces = false) {
    m::Project pr;
    pr.name = "plate joint";
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "Elastic ground";
    s.model = m::SoilModel::LinearElastic;
    s.E = 5000.0; s.nu = 0.3;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P;
    P.name = "Ground"; P.material = 0;
    P.x = {0, 12, 12, 0};
    P.y = {-6, -6, 0, 0};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    m::PlateMaterial pm;
    pm.name = "Slab"; pm.EA = 1.0e7; pm.EI = 1.0e4;
    pr.plates.push_back(pm);
    const auto add = [&pr, interfaces](const char* name, double x1, double x2) {
        m::StructElement e;
        e.kind = m::StructKind::Plate; e.name = name; e.material = 0;
        e.x1 = x1; e.y1 = 0.0; e.x2 = x2; e.y2 = 0.0;
        e.iface_pos = e.iface_neg = interfaces;
        pr.structs.push_back(e);
    };
    if (two_pieces) { add("A", 3.0, 6.0); add("B", 6.0, 9.0); }
    else add("A", 3.0, 9.0);
    m::Load L;
    L.kind = m::LoadKind::Point; L.name = "P";
    L.x1 = 6.0; L.y1 = 0.0; L.qx1 = 0.0; L.qy1 = -100.0;
    pr.loads.push_back(L);
    // Both meshes carry the same vertices: a zero line load from x = 5 to 6 m is a mesh constraint
    // in both, so the one-piece plate has a node at x = 6 too and the only difference between the
    // two models is how the plate is drawn.
    m::Load Z;
    Z.kind = m::LoadKind::Distributed; Z.name = "Vertex at the joint";
    Z.x1 = 5.0; Z.y1 = 0.0; Z.x2 = 6.0; Z.y2 = 0.0;
    pr.loads.push_back(Z);
    return pr;
}

struct Moments { bool ok = false; double at_load = 0.0, peak = 0.0, x_min = 1e300, x_max = -1e300; };

Moments solve(const m::Project& pr) {
    Moments r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return r; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res.back().ok) {
        std::printf("      (solve: %s)\n", res.empty() ? "no result" : res.back().message.c_str());
        return r;
    }
    double best = 1e300;
    for (const auto& sf : res.back().struct_forces) {
        if (sf.kind != 0) continue;
        for (const auto& st : sf.stations) {
            r.peak = std::fmax(r.peak, std::fabs(st.M));
            r.x_min = std::fmin(r.x_min, st.x);
            r.x_max = std::fmax(r.x_max, st.x);
            const double d = std::fabs(st.x - 6.0);
            if (d < best - 1e-12) { best = d; r.at_load = st.M; }
        }
    }
    r.ok = true;
    return r;
}

}  // namespace

int main() {
    std::printf("KV-STR-015: plates meeting at a node are joined rigidly\n");
    const Moments one = solve(model(false)), two = solve(model(true));
    check(one.ok && two.ok, "both models solve");
    if (!one.ok || !two.ok) return 1;
    std::printf("   one piece : M under the load %.6f, max|M| %.6f kNm/m\n", one.at_load, one.peak);
    std::printf("   two pieces: M under the load %.6f, max|M| %.6f kNm/m\n", two.at_load, two.peak);
    const double e1 = std::fabs(two.at_load - one.at_load) / std::fabs(one.at_load);
    const double e2 = std::fabs(two.peak - one.peak) / one.peak;
    std::printf("   relative differences %.2e (under the load), %.2e (largest)\n", e1, e2);
    check(std::fabs(one.at_load) > 10.0, "the load bends the plate (the check below has teeth)");
    check(e1 <= 1e-6 && e2 <= 1e-6,
          "a plate drawn in two pieces carries the moments of the plate drawn in one");

    // The same with interfaces on both sides: the walls then move on their own degrees of freedom,
    // and their ends are joined as one point. The seam used to run a metre past a wall's end, so
    // the force report of each piece also ran into the next (and past x = 9 m).
    const Moments one_i = solve(model(false, true)), two_i = solve(model(true, true));
    check(one_i.ok && two_i.ok, "both models with interfaces solve");
    if (one_i.ok && two_i.ok) {
        const double ei = std::fabs(two_i.at_load - one_i.at_load) / std::fabs(one_i.at_load);
        std::printf("   with interfaces: M under the load %.6f (one piece) vs %.6f (two), relative "
                    "%.2e; stations x = %.4f .. %.4f and %.4f .. %.4f\n", one_i.at_load,
                    two_i.at_load, ei, one_i.x_min, one_i.x_max, two_i.x_min, two_i.x_max);
        check(ei <= 1e-3,
              "with interfaces too, two pieces carry the moment of one (1e-3: the joint's "
              "interface nodes differ)");
        check(one_i.x_min >= 3.0 - 1e-9 && one_i.x_max <= 9.0 + 1e-9 && two_i.x_min >= 3.0 - 1e-9 &&
                  two_i.x_max <= 9.0 + 1e-9,
              "and the walls are built -- and report forces -- only where they are drawn");
    }
    if (g_failures == 0) {
        std::printf("\nOK: the joint is rigid\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
