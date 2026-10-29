// A geogrid and an anchor in an axisymmetric model, through the whole program.
//
// A GEOGRID OF REVOLUTION is a membrane that stretches round its circumference as well as along its
// meridian: the hoop strain u_r / r carries a hoop force N_t = EA u_r / r beside the meridional
// N = EA du_s/ds, both tension only, the grid's two directions with the same stiffness and no
// Poisson coupling. Two closed forms, in ground a hundred thousand times softer than the grid:
//   (a) a vertical grid at radius R under an outward pressure p: a free cylinder, N_t = p R and
//       u_r = p R^2 / EA;
//   (b) a horizontal grid from the axis to r = b, its outer edge free, under an outward load q per
//       unit area of the grid: the membrane equation d(r N_r)/dr - N_t + q r = 0 with N_r = EA u',
//       N_t = EA u / r, bounded at the axis and N_r(b) = 0, is solved by
//       u = q r (2 b - r) / (3 EA), N_r = 2 q (b - r) / 3, N_t = q (2 b - r) / 3 -- a quadratic,
//       which the quadratic element holds exactly. In plane strain the same grid is a bar,
//       N = q (b - x), and has no hoop.
// AN ANCHOR IN AXISYMMETRY is a row of bars round a circle at the spacing Ls measured at its soil end,
// r / Ls of them per radian: a fixed-end anchor of length L under a ring load F per radian at that
// end stretches F Ls L / (EA r), and each anchor carries F Ls / r; its capacity and its lock-off
// force are per anchor, as in plane strain.
//
// verify: KV-STR-031
//   oracle:   closed_form
//   source:   Timoshenko, S.P. & Woinowsky-Krieger, S. (1959). Theory of Plates and Shells, 2nd ed. McGraw-Hill, New York -- ch. 15, the cylinder under internal pressure (N_t = p R), and the membrane equilibrium of a surface of revolution, d(r N_r)/dr - N_t + q r = 0, solved here with N_r = EA u', N_t = EA u / r for a grid with no Poisson coupling
//   locator:  axisymmetric project, ground 0 <= r <= 6 m, 0 <= z <= 4 m, E = 1 kPa, nu = 0.3, weightless, smooth base, the axis held radially; a geogrid EA = 1e5 kN/m (a) along r = 2 m, 0 <= z <= 4 m under 10 kPa in +r, (b) along z = 2 m, 0 <= r <= 4 m under 10 kPa in +r; tri6 and tri15 at 0.25 m; (b) also in plane strain
//   quantity: the hoop and the meridional force [kN/m] at the grid's stations; the radial displacement [m]
//   expected: (a) N_t = p R = 20 kN/m, u_r = p R^2 / EA = 4e-4 m at mid-height; (b) N_r = 2 q (b - r) / 3 and N_t = q (2 b - r) / 3 at every station, u(r) = q r (2 b - r) / (3 EA); in plane strain N = q (b - x) and no hoop
//   band:     1e-3 relative on (a) (the ground, 1e-5 of the grid's stiffness, and the grid's ends); 1e-4 of q b on (b) (an exact state of the element, less the ground)
//
// verify: KV-STR-032
//   oracle:   closed_form
//   source:   the statics and the stiffness of a row of bars: n bars of axial stiffness EA / L share a load F, each stretching F L / (n EA) and carrying F / n (any structural analysis text, e.g. Hibbeler, R.C. Structural Analysis); round a circle of radius r at the spacing Ls there are n = r / Ls per radian
//   locator:  the ground of KV-STR-031; a fixed-end anchor EA = 1e5 kN from (2, 2) to a fixed point (2, 6), Ls = 1 m, under a point load of -10 kN per radian at (2, 2); then with a tensile capacity of 3 kN per anchor; then with a lock-off force of 6 kN per anchor and no load; each also in plane strain (Ls = 1 m, the load per metre)
//   quantity: the anchor force per anchor [kN] and the settlement of the loaded node [m]
//   expected: axisymmetric: N = F Ls / r = 5 kN, settlement F Ls L / (EA r) = 2e-4 m; capped at 3 kN and yielded; 6 kN at lock-off. Plane strain: N = F Ls = 10 kN, settlement 4e-4 m; 3 kN yielded; 6 kN
//   band:     1e-4 relative (the ground carries 1e-5 of the load beside the anchor); the capacity and the lock-off force exactly
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

constexpr double H = 4.0;

m::Project ground(bool axisym) {
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
    // Four regions meeting at (2, 2), so that the anchor's soil end and the ring load are on a node.
    const double xs[3] = {0.0, 2.0, 6.0}, ys[3] = {0.0, 2.0, H};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            m::SoilPolygon g;
            g.name = "Ground"; g.material = 0;
            g.x = {xs[i], xs[i + 1], xs[i + 1], xs[i]}; g.y = {ys[j], ys[j], ys[j + 1], ys[j + 1]};
            g.edge_bc = {j == 0 ? (int)m::BCType::VerticallyFixed : 0, 0, 0,
                         i == 0 ? (int)m::BCType::HorizontallyFixed : 0};
            pr.polygons.push_back(g);
        }
    return pr;
}

// The initial phase with nothing on, then one Plastic phase with the structure and the load.
std::vector<katai::app::SolveResult> solve(m::Project& pr, int order, std::string& msg) {
    pr.initial.struct_active.assign(pr.structs.size(), 0);
    pr.initial.load_active.assign(pr.loads.size(), 0);
    m::Phase ph = pr.initial;
    ph.name = "Load"; ph.type = m::PhaseType::Plastic;
    ph.struct_active.assign(pr.structs.size(), 1);
    ph.load_active.assign(pr.loads.size(), 1);
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr, 0.25, order);
    if (!M.ok) { msg = M.message; return {}; }
    auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) { msg = res.empty() ? "no result" : res.back().message; return {}; }
    return res;
}

void add_geogrid(m::Project& pr, double x1, double y1, double x2, double y2, double EA, double q) {
    m::GeogridMaterial gm;
    gm.name = "Grid"; gm.EA = EA;
    pr.geogrids.push_back(gm);
    m::StructElement e;
    e.kind = m::StructKind::Geogrid; e.name = "Grid"; e.material = 0;
    e.x1 = x1; e.y1 = y1; e.x2 = x2; e.y2 = y2;
    pr.structs.push_back(e);
    m::Load l;
    l.kind = m::LoadKind::Distributed; l.name = "Pressure";
    l.x1 = x1; l.y1 = y1; l.x2 = x2; l.y2 = y2; l.qx1 = l.qx2 = q; l.qy1 = l.qy2 = 0.0;
    pr.loads.push_back(l);
}

struct AnchorOut { bool ok = false, yielded = false; double N = 0.0, uy = 0.0; std::string msg; };

AnchorOut run_anchor(bool axisym, int order, double F, double Fmax, double prestress) {
    m::Project pr = ground(axisym);
    m::AnchorMaterial am;
    am.name = "Bar"; am.EA = 1.0e5; am.Lspacing = 1.0; am.prestress = prestress;
    if (Fmax > 0.0) { am.elastoplastic = true; am.Fmax_tens = Fmax; }
    pr.anchors.push_back(am);
    m::StructElement e;
    e.kind = m::StructKind::Anchor; e.name = "Anchor"; e.material = 0;
    e.x1 = 2.0; e.y1 = 2.0; e.x2 = 2.0; e.y2 = 6.0;
    pr.structs.push_back(e);
    if (F != 0.0) {
        m::Load l;
        l.kind = m::LoadKind::Point; l.name = "Ring load";
        l.x1 = 2.0; l.y1 = 2.0; l.qx1 = 0.0; l.qy1 = F;
        pr.loads.push_back(l);
    }
    AnchorOut o;
    const auto res = solve(pr, order, o.msg);
    if (res.empty()) return o;
    for (const auto& sf : res[1].struct_forces)
        if (sf.kind == 1 && !sf.stations.empty()) {
            o.ok = true; o.N = sf.stations[0].N; o.uy = sf.stations[0].uy; o.yielded = sf.yielded;
        }
    return o;
}

}  // namespace

int main() {
    std::printf("== A geogrid and an anchor in an axisymmetric model (KV-STR-031, KV-STR-032) ==\n");
    constexpr double EA = 1.0e5, p = 10.0, R = 2.0, b = 4.0;

    for (int order : {6, 15}) {
        // (a) a cylinder of grid under internal pressure
        {
            m::Project pr = ground(true);
            add_geogrid(pr, R, 0.0, R, H, EA, p);
            std::string msg;
            const auto res = solve(pr, order, msg);
            double Nt = 0.0, ur = 0.0, best = 1e300;
            for (const auto& sf : res.empty() ? std::vector<katai::core::StructForce>{} : res[1].struct_forces)
                for (const auto& st : sf.stations)
                    if (std::fabs(st.y - 0.5 * H) < best) { best = std::fabs(st.y - 0.5 * H); Nt = st.Nt; ur = st.ux; }
            if (res.empty()) std::printf("      (%s)\n", msg.substr(0, 300).c_str());
            std::printf("      %d-noded, cylinder: N_t %.6f kN/m (p R %.1f), u_r %.6e m (%.6e)\n", order, Nt, p * R, ur, p * R * R / EA);
            char lbl[200];
            std::snprintf(lbl, sizeof lbl, "%d-noded: a cylinder of geogrid carries its pressure in the hoop, N_t = p R, u_r = p R^2 / EA", order);
            check(!res.empty() && rel(Nt, p * R) < 1e-3 && rel(ur, p * R * R / EA) < 1e-3, lbl);
        }
        // (b) a disk of grid from the axis, its edge free, under an outward load
        for (bool axisym : {true, false}) {
            m::Project pr = ground(axisym);
            add_geogrid(pr, 0.0, 2.0, b, 2.0, EA, p);
            std::string msg;
            const auto res = solve(pr, order, msg);
            double eN = 0.0, eNt = 0.0, eu = 0.0;
            int n = 0;
            for (const auto& sf : res.empty() ? std::vector<katai::core::StructForce>{} : res[1].struct_forces)
                for (const auto& st : sf.stations) {
                    const double r = st.x;
                    const double N = axisym ? 2.0 * p * (b - r) / 3.0 : p * (b - r);
                    const double Nt = axisym ? p * (2.0 * b - r) / 3.0 : 0.0;
                    const double u = axisym ? p * r * (2.0 * b - r) / (3.0 * EA) : p * r * (2.0 * b - r) / (2.0 * EA);
                    eN = std::max(eN, std::fabs(st.N - N));
                    eNt = std::max(eNt, std::fabs(st.Nt - Nt));
                    eu = std::max(eu, std::fabs(st.ux - u) / (p * b * b / EA));
                    ++n;
                }
            if (res.empty()) std::printf("      (%s)\n", msg.substr(0, 300).c_str());
            std::printf("      %d-noded, disk, %s: %d stations, max |N - N*| %.3e, |N_t - N_t*| %.3e kN/m, |u - u*| %.3e of q b^2 / EA\n",
                        order, axisym ? "axisymmetric" : "plane strain", n, eN, eNt, eu);
            char lbl[260];
            std::snprintf(lbl, sizeof lbl, axisym
                ? "%d-noded: a disk of geogrid, N_r = 2 q (b - r) / 3 and N_t = q (2 b - r) / 3 at every station, u = q r (2 b - r) / (3 EA)"
                : "%d-noded: the same grid in plane strain is a bar, N = q (b - x), and has no hoop", order);
            check(!res.empty() && n >= 8 && eN < 1e-4 * p * b && eNt < 1e-4 * p * b && eu < 1e-4, lbl);
        }
    }

    // The anchor: a row of bars round a circle
    {
        constexpr double F = -10.0, r = 2.0, Ls = 1.0, L = 4.0, EAa = 1.0e5;
        for (int order : {6, 15}) {
            const AnchorOut a = run_anchor(true, order, F, 0.0, 0.0);
            const AnchorOut ps = run_anchor(false, order, F, 0.0, 0.0);
            if (!a.ok) std::printf("      (%s)\n", a.msg.substr(0, 300).c_str());
            if (!ps.ok) std::printf("      (%s)\n", ps.msg.substr(0, 300).c_str());
            const double Na = -F * Ls / r, ua = F * Ls * L / (EAa * r);
            const double Np = -F * Ls, up = F * Ls * L / EAa;
            std::printf("      %d-noded: axisymmetric N %.6f kN (%.3f), u %.6e m (%.6e); plane strain N %.6f (%.3f), u %.6e (%.6e)\n",
                        order, a.N, Na, a.uy, ua, ps.N, Np, ps.uy, up);
            char lbl[260];
            std::snprintf(lbl, sizeof lbl, "%d-noded: a row of anchors round a circle, r / Ls per radian: each carries F Ls / r and the ring settles F Ls L / (EA r)", order);
            check(a.ok && rel(a.N, Na) < 1e-4 && rel(a.uy, ua) < 1e-4, lbl);
            std::snprintf(lbl, sizeof lbl, "%d-noded: the same anchor in plane strain carries F Ls and settles F Ls L / EA", order);
            check(ps.ok && rel(ps.N, Np) < 1e-4 && rel(ps.uy, up) < 1e-4, lbl);
        }
        const AnchorOut c = run_anchor(true, 6, F, 3.0, 0.0);
        const AnchorOut cp = run_anchor(false, 6, F, 3.0, 0.0);
        std::printf("      capacity 3 kN per anchor: axisymmetric N %.6f (yielded %d), plane strain N %.6f (yielded %d)\n",
                    c.N, (int)c.yielded, cp.N, (int)cp.yielded);
        check(c.ok && c.yielded && std::fabs(c.N - 3.0) < 1e-9 && cp.ok && cp.yielded && std::fabs(cp.N - 3.0) < 1e-9,
              "the capacity is per anchor in both: each anchor stops at 3 kN and is reported yielded");
        const AnchorOut l = run_anchor(true, 6, 0.0, 0.0, 6.0);
        const AnchorOut lp = run_anchor(false, 6, 0.0, 0.0, 6.0);
        std::printf("      lock-off 6 kN per anchor: axisymmetric N %.6f, plane strain N %.6f\n", l.N, lp.N);
        check(l.ok && std::fabs(l.N - 6.0) < 1e-9 && lp.ok && std::fabs(lp.N - 6.0) < 1e-9,
              "the lock-off force is per anchor in both: 6 kN at the end of the phase that installs it");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
