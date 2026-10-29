// A plate in an axisymmetric model, through the whole program: a tank wall under internal
// pressure.
//
// The element (test_axisym_plate, KV-STR-029) is a shell of revolution; this asks whether the
// program delivers it -- the plate drawn in an axisymmetric project, its hoop terms switched on by
// the driver, the pressure of a distributed load on its line integrated per radian, and the hoop
// force in the structural output. A vertical plate at radius R = 2 m, 4 m high, standing on a
// smooth base inside ground a million times softer than itself, loaded by p = 10 kPa outward along
// its whole height: away from its ends the wall is a free cylinder, N_t = p R (statics) and
// u_r = p R^2 / (EA (1 - nu^2)). The same wall in plane strain is a beam with no hoop at all: it
// carries the pressure in bending, and moves far more -- the check that the hoop is really there.
//
// verify: KV-STR-030
//   oracle:   closed_form
//   source:   Timoshenko, S.P. & Woinowsky-Krieger, S. (1959). Theory of Plates and Shells, 2nd ed. McGraw-Hill, New York, ch. 15 -- the cylinder under internal pressure: N_t = p R, w = p R^2 / (E t) (here EA (1 - nu^2) with the plate's plane strain EA)
//   locator:  axisymmetric project, ground 0 <= r <= 6 m, 0 <= z <= 4 m, E = 1 kPa, nu = 0.3, weightless, smooth base, the axis held radially; a plate EA = 1e6 kN/m, EI = 1e3 kNm/m, nu = 0.3 along r = 2 m; a distributed load of 10 kPa in +r on the plate's line; tri6 and tri15 at 0.25 m
//   quantity: the hoop force at mid-height [kN/m] and the radial displacement there [m]
//   expected: N_t = p R = 20 kN/m; u_r = p R^2 / (EA (1 - nu^2)) = 4.3956e-5 m; in plane strain, a displacement at least ten times that
//   band:     1e-3 relative (the ground's stiffness, a millionth of the wall's, and the ends' bending, which decays over sqrt(R t) ~ 0.2 m)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

struct Out { bool ok = false; double Nt = 0.0, ur = 0.0; std::string msg; };

Out run(bool axisym, int order) {
    constexpr double R = 2.0, H = 4.0, p = 10.0;
    m::Project pr;
    pr.axisymmetric = axisym;
    pr.x_min = 0; pr.x_max = 6; pr.y_min = 0; pr.y_max = H;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    m::Material s;
    s.name = "Soft"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::NonPorous;
    s.E = 1.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 0.0;
    s.k0_auto = false; s.k0 = 0.5;
    pr.materials.push_back(s);
    m::SoilPolygon g;
    g.name = "Ground"; g.material = 0; g.x = {0, 6, 6, 0}; g.y = {0, 0, H, H};
    g.edge_bc = {(int)m::BCType::VerticallyFixed, 0, 0, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(g);
    m::PlateMaterial pm;
    pm.name = "Wall"; pm.EA = 1.0e6; pm.EI = 1.0e3; pm.nu = 0.3; pm.w = 0.0;
    pr.plates.push_back(pm);
    m::StructElement w;
    w.kind = m::StructKind::Plate; w.name = "Tank wall"; w.material = 0;
    w.x1 = R; w.y1 = 0.0; w.x2 = R; w.y2 = H;
    pr.structs.push_back(w);
    m::Load q;
    q.kind = m::LoadKind::Distributed; q.name = "Pressure";
    q.x1 = R; q.y1 = 0.0; q.x2 = R; q.y2 = H; q.qx1 = q.qx2 = p; q.qy1 = q.qy2 = 0.0;
    pr.loads = {q};
    pr.initial.struct_active = {1};
    pr.initial.load_active = {0};
    m::Phase ph = pr.initial;
    ph.name = "Fill"; ph.type = m::PhaseType::Plastic; ph.load_active = {1}; ph.struct_active = {1};
    pr.phases.push_back(ph);
    Out o;
    const auto M = katai::app::mesh_from_project(pr, 0.25, order);
    if (!M.ok) { o.msg = M.message; return o; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) { o.msg = res.empty() ? "no result" : res.back().message; return o; }
    o.ok = true;
    // the station nearest mid-height
    double best = 1e300;
    for (const auto& sf : res[1].struct_forces)
        for (const auto& st : sf.stations)
            if (std::fabs(st.y - 0.5 * H) < best) { best = std::fabs(st.y - 0.5 * H); o.Nt = st.Nt; o.ur = st.ux; }
    return o;
}

}  // namespace

int main() {
    std::printf("== A plate in an axisymmetric model: a tank wall (KV-STR-030) ==\n");
    constexpr double R = 2.0, p = 10.0, EA = 1.0e6, nu = 0.3;
    const double want_u = p * R * R / (EA * (1.0 - nu * nu));
    for (int order : {6, 15}) {
        const Out a = run(true, order);
        if (!a.ok) std::printf("      (%s)\n", a.msg.substr(0, 300).c_str());
        std::printf("      %d-noded, axisymmetric: N_t %.6f kN/m (p R %.3f), u_r %.6e m (%.6e)\n", order, a.Nt, p * R, a.ur, want_u);
        char lbl[200];
        std::snprintf(lbl, sizeof lbl, "%d-noded: the tank wall carries its pressure in the hoop, N_t = p R, and expands p R^2 / (EA (1 - nu^2))", order);
        check(a.ok && rel(a.Nt, p * R) < 1e-3 && rel(a.ur, want_u) < 1e-3, lbl);
    }
    {
        const Out b = run(false, 6);
        std::printf("      plane strain, the same wall: u %.6e m (%.0f times the shell's)\n", b.ur, b.ur / want_u);
        check(b.ok && b.ur > 10.0 * want_u && b.Nt == 0.0, "in plane strain the same wall has no hoop: it bends, and moves at least ten times as far");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
