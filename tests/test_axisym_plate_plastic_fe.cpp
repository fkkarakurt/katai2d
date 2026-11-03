// A shell of revolution with a state, through the whole program: a tank wall that yields in its
// hoop, and a tank wall whose section a phase changes.
//
// The element (test_axisym_plate_plastic, KV-STR-034) caps the hoop force of a shell of revolution
// at Np apart from the meridional diamond; this asks whether the program delivers it -- an
// elastoplastic plate accepted in an axisymmetric project, its state carried per radian, and the
// hoop force in the structural output read from that state. A vertical plate at radius R = 2 m,
// 4 m high, standing on a smooth base inside ground a million times softer than itself:
//   (a) expanded radially by a prescribed displacement of three times its yield displacement: the
//       wall is a free cylinder at its hoop capacity, N_t = Np along its height and N_s = 0; with no
//       capacity the same expansion gives three times as much;
//   (b) under an internal pressure p (KV-STR-030), then a phase that gives the wall a section twice
//       as stiff and changes nothing else: the wall keeps its forces, N_t = p R (statics), and the
//       phase moves nothing. Until 2026-10-30 a plate whose section a phase changed went on to be
//       solved without its hoop and without its radius in an axisymmetric model: a plane strain
//       plate inside a shell.
//   (c) the same expansion of a geogrid of revolution with a tensile capacity Np: its hoop force
//       stops at Np (the hoop has a plastic strain of its own); with no capacity it carries 3 Np.
//       Until 2026-10-30 a geogrid with a capacity was refused in axisymmetry.
//
// verify: KV-STR-035
//   oracle:   closed_form
//   source:   Timoshenko, S.P. & Woinowsky-Krieger, S. (1959). Theory of Plates and Shells, 2nd ed. McGraw-Hill, New York, ch. 15 -- the free cylinder: N_t = p R from statics, N_t = EA (1 - nu^2) u_r / R with N_s = 0
//   locator:  axisymmetric project, ground 0 <= r <= 6 m, 0 <= z <= 4 m, E = 1 kPa (1e-3 kPa under the geogrid), nu = 0.3, weightless, smooth base, the axis held radially; a plate EA = 1e6 kN/m, EI = 1e3 kNm/m, nu = 0.3 along r = 2 m; (a) Np = 50 kN/m, Mp = 10 kNm/m, u_r = 3 Np R / (EA (1 - nu^2)) prescribed on the plate's line; (b) a distributed load of 10 kPa in +r, then the section EA = 2e6, EI = 2e3; tri6 and tri15 at 0.25 m
//   quantity: the hoop and meridional force at mid-height [kN/m]; (b) the largest displacement of the second phase [m]; (c) the geogrid's hoop force at mid-height [kN/m]
//   expected: (a) N_t = Np = 50 kN/m, N_s = 0 (3 Np with no capacity); (b) N_t = p R = 20 kN/m after the change, and no displacement in that phase; (c) N_t = Np = 30 kN/m (3 Np with no capacity), a geogrid EA = 1000 kN/m
//   band:     1e-3 relative (the ground's stiffness, a millionth of the wall's); (b) the second phase's displacement below 1e-3 of the first's
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

constexpr double R = 2.0, H = 4.0, EA = 1.0e6, nu = 0.3, Np = 50.0, p = 10.0;

m::Project tank() {
    m::Project pr;
    pr.axisymmetric = true;
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
    pm.name = "Wall"; pm.EA = EA; pm.EI = 1.0e3; pm.nu = nu; pm.w = 0.0;
    pr.plates.push_back(pm);
    m::StructElement w;
    w.kind = m::StructKind::Plate; w.name = "Tank wall"; w.material = 0;
    w.x1 = R; w.y1 = 0.0; w.x2 = R; w.y2 = H;
    pr.structs.push_back(w);
    pr.initial.struct_active = {1};
    return pr;
}

// the hoop and meridional force at the station nearest mid-height of a phase's result
void mid(const katai::app::SolveResult& r, double& Nt, double& Ns) {
    double best = 1e300;
    for (const auto& sf : r.struct_forces)
        for (const auto& st : sf.stations)
            if (std::fabs(st.y - 0.5 * H) < best) { best = std::fabs(st.y - 0.5 * H); Nt = st.Nt; Ns = st.N; }
}

void test_yield(int order) {
    const double uy = Np * R / (EA * (1.0 - nu * nu));
    for (bool capped : {true, false}) {
        m::Project pr = tank();
        if (capped) { pr.plates[0].elastoplastic = true; pr.plates[0].Np = Np; pr.plates[0].Mp = 10.0; }
        m::PrescribedDisp d;
        d.name = "Expansion"; d.x1 = R; d.y1 = 0.0; d.x2 = R; d.y2 = H;
        d.set_ux = true; d.ux = 3.0 * uy; d.set_uy = false;
        pr.disps = {d};
        pr.initial.disp_active = {0};
        m::Phase ph = pr.initial;
        ph.name = "Expand"; ph.type = m::PhaseType::Plastic; ph.disp_active = {1};
        pr.phases.push_back(ph);
        const auto M = katai::app::mesh_from_project(pr, 0.25, order);
        const auto res = M.ok ? katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure))
                              : std::vector<katai::app::SolveResult>{};
        const bool ok = res.size() == 2 && res[1].ok;
        if (!ok) std::printf("      (%s)\n", (M.ok ? (res.empty() ? std::string("no result") : res.back().message) : M.message).substr(0, 300).c_str());
        double Nt = 0.0, Ns = 0.0;
        if (ok) mid(res[1], Nt, Ns);
        std::printf("      %d-noded, %s: N_t %.6f kN/m, N_s %.2e kN/m\n", order, capped ? "Np = 50" : "no capacity", Nt, Ns);
        char lbl[200];
        if (capped) {
            std::snprintf(lbl, sizeof lbl, "(a) %d-noded: expanded three times past its yield, the wall holds N_t = Np with N_s = 0", order);
            check(ok && rel(Nt, Np) < 1e-3 && std::fabs(Ns) < 1e-3 * Np, lbl);
        } else {
            std::snprintf(lbl, sizeof lbl, "(a) %d-noded: with no capacity the same expansion carries 3 Np", order);
            check(ok && rel(Nt, 3.0 * Np) < 1e-3, lbl);
        }
    }
}

void test_section_change(int order) {
    m::Project pr = tank();
    m::PlateMaterial thick = pr.plates[0];
    thick.name = "Thicker wall"; thick.EA = 2.0 * EA; thick.EI = 2.0e3;
    pr.plates.push_back(thick);
    m::Load q;
    q.kind = m::LoadKind::Distributed; q.name = "Pressure";
    q.x1 = R; q.y1 = 0.0; q.x2 = R; q.y2 = H; q.qx1 = q.qx2 = p; q.qy1 = q.qy2 = 0.0;
    pr.loads = {q};
    pr.initial.load_active = {0};
    m::Phase fill = pr.initial;
    fill.name = "Fill"; fill.type = m::PhaseType::Plastic; fill.load_active = {1};
    pr.phases.push_back(fill);
    m::Phase thicker = fill;
    thicker.name = "Thicker wall"; thicker.struct_material = {1};
    pr.phases.push_back(thicker);
    const auto M = katai::app::mesh_from_project(pr, 0.25, order);
    const auto res = M.ok ? katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure))
                          : std::vector<katai::app::SolveResult>{};
    const bool ok = res.size() == 3 && res[1].ok && res[2].ok;
    if (!ok) std::printf("      (%s)\n", (M.ok ? (res.empty() ? std::string("no result") : res.back().message) : M.message).substr(0, 300).c_str());
    double Nt1 = 0.0, Ns1 = 0.0, Nt2 = 0.0, Ns2 = 0.0;
    if (ok) { mid(res[1], Nt1, Ns1); mid(res[2], Nt2, Ns2); }
    const double d1 = ok ? res[1].max_disp : 0.0, d2 = ok ? res[2].max_disp : 0.0;
    std::printf("      %d-noded: N_t %.6f -> %.6f kN/m (p R %.1f); max |u| of the phases %.4e, %.4e m\n",
                order, Nt1, Nt2, p * R, d1, d2);
    char lbl[200];
    std::snprintf(lbl, sizeof lbl, "(b) %d-noded: a wall given a stiffer section keeps its hoop force, N_t = p R, and does not move", order);
    check(ok && rel(Nt1, p * R) < 1e-3 && rel(Nt2, p * R) < 1e-3 && d2 < 1e-3 * d1, lbl);
}

void test_geogrid(int order) {
    constexpr double EAg = 1000.0, Npg = 30.0;
    const double uy = Npg * R / EAg;   // no Poisson coupling in a grid: N_t = EA u_r / R
    for (bool capped : {true, false}) {
        m::Project pr = tank();
        pr.materials[0].E = 1.0e-3;    // a millionth of the grid's stiffness, as for the wall
        pr.structs.clear();
        pr.plates.clear();
        m::GeogridMaterial gm;
        gm.name = "Grid"; gm.EA = EAg;
        if (capped) { gm.elastoplastic = true; gm.Np = Npg; }
        pr.geogrids.push_back(gm);
        m::StructElement g;
        g.kind = m::StructKind::Geogrid; g.name = "Grid cylinder"; g.material = 0;
        g.x1 = R; g.y1 = 0.0; g.x2 = R; g.y2 = H;
        pr.structs.push_back(g);
        m::PrescribedDisp d;
        d.name = "Expansion"; d.x1 = R; d.y1 = 0.0; d.x2 = R; d.y2 = H;
        d.set_ux = true; d.ux = 3.0 * uy; d.set_uy = false;
        pr.disps = {d};
        pr.initial.disp_active = {0};
        m::Phase ph = pr.initial;
        ph.name = "Expand"; ph.type = m::PhaseType::Plastic; ph.disp_active = {1};
        pr.phases.push_back(ph);
        const auto M = katai::app::mesh_from_project(pr, 0.25, order);
        const auto res = M.ok ? katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure))
                              : std::vector<katai::app::SolveResult>{};
        const bool ok = res.size() == 2 && res[1].ok;
        if (!ok) std::printf("      (%s)\n", (M.ok ? (res.empty() ? std::string("no result") : res.back().message) : M.message).substr(0, 300).c_str());
        double Nt = 0.0, Ns = 0.0;
        if (ok) mid(res[1], Nt, Ns);
        std::printf("      %d-noded geogrid, %s: N_t %.6f kN/m, N %.2e kN/m\n", order, capped ? "Np = 30" : "no capacity", Nt, Ns);
        char lbl[200];
        if (capped)
            std::snprintf(lbl, sizeof lbl, "(c) %d-noded: a geogrid of revolution expanded past its capacity holds N_t = Np", order);
        else
            std::snprintf(lbl, sizeof lbl, "(c) %d-noded: with no capacity the same grid carries 3 Np", order);
        check(ok && rel(Nt, capped ? Npg : 3.0 * Npg) < 1e-3, lbl);
    }
}

}  // namespace

int main() {
    std::printf("== A shell of revolution with a state, through the program (KV-STR-035) ==\n");
    for (int order : {6, 15}) {
        test_yield(order);
        test_section_change(order);
        test_geogrid(order);
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
