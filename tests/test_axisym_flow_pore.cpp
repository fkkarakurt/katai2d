// The pore pressures of a phase taken from a groundwater flow, and region water, in an
// axisymmetric model. Both were refused there (K2D-S002) because the loads they make were
// integrated as plane strain; they are now per radian, and the pore load carries the hoop
// component the phreatic one has always carried.
//
// Two identities and one boundary. (a) A level water table with no flow: the steady flow a phase
// solves is the table itself, so the phase with its pore pressures from the flow must be the phase
// with them from the phreatic line -- the same loads, the same displacement. (b) Region water at
// the table's own level: the same again. (c) The axis is a line of symmetry: no water crosses it,
// so it is closed, not held at the table's head as a model side is. With the table at 6 m on the
// axis and 8 m at the outer radius the only head boundary is the outer one, and the steady head is
// 8 m everywhere; held at the axis it would be 6 m there.
//
// verify: KV-FLW-013
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, phases[i].flowpore and pwater / pwlevel): a phase's pore pressures from a steady groundwater flow on its active ground, below the water table every boundary of that ground held at the table's head and the base closed, or from region water levels; in an axisymmetric model the axis x = 0 is a line of symmetry across which no water flows. The independent path is the same phase with its pore pressures from the phreatic line, which the static solver has integrated per radian since axisymmetry was introduced
//   locator:  a cylinder of radius 10 m and height 10 m (linear elastic, gamma 18 / 20 kN/m3), the water table level at 8 m, a 30 kPa load on the inner 2 m of its top, K0 initial phase; one Plastic phase applying the load with its pore pressures (a) from the phreatic line, (a') from a steady flow, (b) from region water at 8 m. (c): the table at 6 m on the axis and 8 m at the outer radius, pore pressures from a steady flow
//   quantity: the displacement field of the phase [m]; the steady head at the axis [m]
//   expected: (a') and (b) give the field of (a); in (c) the head is the outer table's, 8 m, at every node
//   band:     1e-9 m on the fields (the loads are the same functions integrated by the same rule); 1e-9 m on the head
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

constexpr double kR = 10.0, kH = 10.0, kWl = 8.0;

// mode: 0 phreatic, 1 flow, 2 region water at the table's level
m::Project cylinder(int mode, double wl_axis = kWl) {
    m::Project pr;
    pr.name = "axisymmetric pore pressures";
    pr.axisymmetric = true;
    pr.x_min = 0; pr.x_max = kR; pr.y_min = 0; pr.y_max = kH;
    pr.has_water = true; pr.wx = {0.0, kR}; pr.wy = {wl_axis, kWl};
    pr.mesh.elem_size = 1.0;
    m::Material s; s.name = "Soil"; s.model = m::SoilModel::LinearElastic; s.E = 2e4; s.nu = 0.3;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.kx = s.ky = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Ground"; P.material = 0;
    P.x = {0, kR, kR, 0}; P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    m::Load q; q.kind = m::LoadKind::Distributed; q.name = "Load";
    q.x1 = 0.0; q.y1 = kH; q.x2 = 2.0; q.y2 = kH; q.qy1 = q.qy2 = -30.0;
    pr.loads.push_back(q);
    pr.initial.load_active = {0};
    m::Phase p; p.name = "Load"; p.type = m::PhaseType::Plastic; p.load_active = {1};
    if (mode == 1) p.pore_from_flow = true;
    if (mode == 2) { p.poly_water = {m::Phase::kRegionWaterLevel}; p.poly_water_level = {kWl}; }
    pr.phases.push_back(p);
    return pr;
}

std::vector<katai::app::SolveResult> solve(const m::Project& pr, const katai::mesh::Mesh& mesh) {
    return katai::app::solve_phases(pr, mesh, katai::app::InitialPhase::K0Procedure);
}
double field_difference(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    if (a.size() != b.size() || a.size() == 0) return 1e300;
    return (a - b).cwiseAbs().maxCoeff();
}
}  // namespace

int main() {
    std::printf("== pore pressures from a flow and region water, axisymmetric (KV-FLW-013) ==\n");
    const m::Project base = cylinder(0);
    const auto M = katai::app::mesh_from_project(base);
    check(M.ok, "the cylinder meshes");
    if (!M.ok) return 1;
    const auto ra = solve(base, M.mesh), rf = solve(cylinder(1), M.mesh), rw = solve(cylinder(2), M.mesh);
    const bool ok = ra.size() == 2 && rf.size() == 2 && rw.size() == 2 && ra[1].ok && rf[1].ok && rw[1].ok;
    for (const auto* r : {&ra, &rf, &rw})
        if (!r->empty() && !r->back().ok) std::printf("      (%s)\n", r->back().message.c_str());
    check(ok, "the phreatic, flow and region-water phases solve in axisymmetry");
    if (ok) {
        const double dmax = ra[1].disp.cwiseAbs().maxCoeff();
        const double df = field_difference(ra[1].disp, rf[1].disp);
        const double dw = field_difference(ra[1].disp, rw[1].disp);
        std::printf("  max|u| %.6e m; flow against phreatic %.3e m, region water against phreatic %.3e m\n",
                    dmax, df, dw);
        check(dmax > 1e-4, "the load moves the ground (the comparison is not vacuous)");
        check(df < 1e-9, "(a') pore pressures from a flow on a level table = from the phreatic line");
        check(dw < 1e-9, "(b) region water at the table's level = the phreatic line");
    }

    // (c) The axis is closed: with the table at 6 m on it and 8 m outside, the only head the flow is
    // held at is the outer one.
    const m::Project tilt = cylinder(1, 6.0);
    const auto rt = solve(tilt, M.mesh);
    check(rt.size() == 2 && rt[1].ok, "the phase with a tilted table solves");
    if (rt.size() == 2 && rt[1].ok) {
        const auto& h = rt[1].flow_head;
        double lo = 1e300, hi = -1e300;
        for (double v : h) { lo = std::min(lo, v); hi = std::max(hi, v); }
        std::printf("  steady head: %.9f .. %.9f m (the outer table 8 m; held at the axis it would be 6 m there)\n",
                    lo, hi);
        check(!h.empty() && std::fabs(lo - kWl) < 1e-9 && std::fabs(hi - kWl) < 1e-9,
              "(c) the axis is a line of symmetry: the head is the outer table's everywhere");
    }
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
