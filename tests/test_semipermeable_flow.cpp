// A SEMI-PERMEABLE SCREEN IN EVERY GROUNDWATER FLOW, against the same screen drawn as soil.
//
// A wall or interface with a hydraulic resistance R = d/k passes q = dh / R per unit area across
// its line. The steady flow solved it only while the ground stayed saturated; with a free surface,
// in a transient flow phase and in the flow of a staged phase it was refused (K2D-S004), so a
// leaky cut-off could not be drawn under a dam, in a drawdown or beside an excavation.
//
// Oracle (independent path): the same screen drawn as what it stands for -- a strip of soil of
// width d whose permeability gives the strip the resistance of the screen plus that of the soil it
// replaces, k_s = d / (R + d/k) -- solved by the ordinary flow elements. Both models carry the strip
// as a region (soil of permeability k in one, k_s in the other), so they share the geometry and
// differ only in where the resistance sits: on a line, or spread over the strip. In the
// unsaturated zone the screen's conductance follows the relative permeability of the soil it is
// drawn in, as the strip's does.
//
// verify: KV-FLW-011
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, structs[i].flow_barrier / hyd_res): a semi-permeable barrier passes q = dh / R per unit area, R = d/k [day]; Darcy's law in series across a strip of width d: the strip of permeability d / (R + d/k) has the resistance of the screen plus that of the soil it replaces
//   locator:  a 20 x 6 m block of k = 1 m/day between a reservoir at 5 m (left) and one at 1 m (right), the Head edges drawn over the whole faces, a 0.4 m strip at mid-length, R = 20 days, 6-noded elements of 0.1 m ((A), (B)) and 0.03 m ((C)): (A) steady flow with a free surface; (B) a transient flow phase lowering the right reservoir from 5 to 1 m, 2 days; (C) a Plastic phase taking its pore pressures from the steady flow of its own water (5 m behind the screen, 1 m in front of it)
//   quantity: (A) discharge [m3/day per m] and the head on either side of the screen [m]; (B) the head at four points at the end of the phase [m]; (C) the pore pressure at four points [kPa]
//   expected: the screen model equals the strip model
//   band:     2% of the discharge and 2% of the head drop (4 m) on the heads; 2% of the largest pore pressure difference on (C) -- the strip spreads over 0.4 m what the screen holds on a line, and the mesh resolves the two differently
#include <katai/jobs/driver.hpp>
#include <katai/jobs/flow_driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kL = 20.0, kH = 6.0, kH1 = 5.0, kH2 = 1.0, kK = 1.0, kR = 20.0;
constexpr double kX0 = 9.8, kX1 = 10.2, kD = kX1 - kX0;
// Element sizes. The screen and the strip agree to 0.12% (discharge) and 0.6% of the drop
// (transient heads) at 0.1 m; the staged phase's pore pressures, which compare a seam held at the
// pit's water through its resistance with the strip, agree to 1.2% of their scale at 0.03 m, 1.6%
// at 0.06 and 1.75% at 0.1 -- too close to the 2% the check allows, so that part keeps 0.03.
// (At 0.03 throughout the test took 750 s, the whole suite's wall clock.)
constexpr double kMeshOpen = 0.1, kMeshStaged = 0.03;

// `screen`: the strip is soil and a semi-permeable interface runs down its middle; otherwise the
// strip carries the screen's resistance itself.
m::Project block(bool screen) {
    m::Project pr;
    m::Material soil; soil.name = "Soil"; soil.model = m::SoilModel::LinearElastic;
    soil.E = 2.0e4; soil.nu = 0.3; soil.gamma_unsat = 18.0; soil.gamma_sat = 20.0;
    soil.kx = kK; soil.ky = kK;
    m::Material strip = soil; strip.name = "Strip";
    if (!screen) strip.kx = strip.ky = kD / (kR + kD / kK);
    pr.materials = {soil, strip};
    const int C = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
    const auto poly = [&](const char* name, int mat, double x0, double x1, int left, int right,
                          double hl, double hr) {
        m::SoilPolygon P; P.name = name; P.material = mat;
        P.x = {x0, x1, x1, x0}; P.y = {0, 0, kH, kH};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
        P.edge_flow = {C, right, C, left};
        P.edge_head = {0.0, hr, 0.0, hl};
        return P;
    };
    pr.polygons = {poly("Upstream", 0, 0, kX0, Hd, C, kH1, 0.0), poly("Strip", 1, kX0, kX1, C, C, 0.0, 0.0),
                   poly("Downstream", 0, kX1, kL, C, Hd, 0.0, kH2)};
    if (screen) {
        m::StructElement w; w.kind = m::StructKind::Interface; w.name = "Screen";
        w.x1 = 0.5 * (kX0 + kX1); w.y1 = -1.0; w.x2 = w.x1; w.y2 = kH + 1.0;
        w.material = -1; w.flow_barrier = 2; w.hydraulic_resistance = kR;
        pr.structs.push_back(w);
    }
    pr.has_water = false;
    return pr;
}

double head_near(const katai::mesh::Mesh& mesh, const std::vector<double>& head, double x, double y) {
    int best = -1; double bd = 1e300;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double d = std::hypot(mesh.x[n] - x, mesh.y[n] - y);
        if (d < bd) { bd = d; best = n; }
    }
    return best >= 0 ? head[(size_t)best] : 0.0;
}

void steady() {
    std::printf("-- (A) steady flow with a free surface --\n");
    double q[2] = {0, 0}, hu[2] = {0, 0}, hd[2] = {0, 0};
    for (int s = 0; s < 2; ++s) {
        const m::Project pr = block(s == 0);
        const auto M = katai::app::mesh_from_project(pr, kMeshOpen, 6);
        if (!M.ok) { check(false, "meshed"); return; }
        const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
        check(F.ok, std::string(s == 0 ? "screen" : "strip") + " model solved" + (F.ok ? "" : " (" + F.message.substr(0, 160) + ")"));
        if (!F.ok) return;
        std::vector<double> h(F.head.data(), F.head.data() + M.mesh.node_count);
        q[s] = F.discharge;
        hu[s] = head_near(M.mesh, h, kX0 - 0.5, 0.5);
        hd[s] = head_near(M.mesh, h, kX1 + 0.5, 0.5);
    }
    std::printf("   discharge %.5f (screen) vs %.5f (strip); head 0.5 m upstream %.4f vs %.4f, "
                "downstream %.4f vs %.4f\n", q[0], q[1], hu[0], hu[1], hd[0], hd[1]);
    check(std::fabs(q[0] - q[1]) < 0.02 * q[1], "the discharge through the screen is the strip's within 2%");
    check(std::fabs(hu[0] - hu[1]) < 0.02 * (kH1 - kH2) && std::fabs(hd[0] - hd[1]) < 0.02 * (kH1 - kH2),
          "the heads on either side are the strip's within 2% of the drop");
    check(hu[0] - hd[0] > 0.5, "the screen holds a head difference (teeth)");
}

void transient() {
    std::printf("-- (B) a transient flow phase --\n");
    std::vector<double> got[2];
    const double pts[4][2] = {{5.0, 0.5}, {9.0, 0.5}, {11.0, 0.5}, {15.0, 0.5}};
    for (int s = 0; s < 2; ++s) {
        m::Project pr = block(s == 0);
        pr.has_water = true; pr.wx = {0, kL}; pr.wy = {kH1, kH1};
        m::Phase p; p.type = m::PhaseType::TransientFlow; p.name = "Drawdown";
        p.duration = 2.0; p.time_steps = 20;
        pr.phases = {p};
        const auto M = katai::app::mesh_from_project(pr, kMeshOpen, 6);
        if (!M.ok) { check(false, "meshed"); return; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 2 && res[1].ok;
        check(ok, std::string(s == 0 ? "screen" : "strip") + " model solved" +
                      (ok ? "" : " (" + (res.empty() ? std::string("no result") : res.back().message.substr(0, 160)) + ")"));
        if (!ok) return;
        for (const auto& p2 : pts) got[s].push_back(head_near(res[1].mesh, res[1].flow_head, p2[0], p2[1]));
    }
    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        std::printf("   (%.0f, %.1f): head %.4f (screen) vs %.4f (strip)\n", pts[i][0], pts[i][1], got[0][i], got[1][i]);
        worst = std::fmax(worst, std::fabs(got[0][i] - got[1][i]));
    }
    check(worst < 0.02 * (kH1 - kH2), "the transient heads are the strip's within 2% of the drop");
    check(got[0][1] - got[0][2] > 0.3, "the screen holds a head difference in time (teeth)");
}

void staged() {
    std::printf("-- (C) the flow of a staged phase --\n");
    std::vector<double> got[2];
    const double pts[4][2] = {{5.0, 0.5}, {9.0, 0.5}, {11.0, 0.5}, {15.0, 0.5}};
    double scale = 0.0;
    for (int s = 0; s < 2; ++s) {
        m::Project pr = block(s == 0);
        for (auto& P : pr.polygons) { P.edge_flow.clear(); P.edge_head.clear(); }
        pr.has_water = true; pr.wx = {0, kL}; pr.wy = {kH1, kH1};
        m::Phase p; p.type = m::PhaseType::Plastic; p.name = "Lowered in front";
        p.pore_from_flow = true;
        p.water_override = true; p.wx = {0, kX0, kX1, kL}; p.wy = {kH1, kH1, kH2, kH2};
        pr.phases = {p};
        const auto M = katai::app::mesh_from_project(pr, kMeshStaged, 6);
        if (!M.ok) { check(false, "meshed"); return; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 2 && res[1].ok;
        check(ok, std::string(s == 0 ? "screen" : "strip") + " model solved" +
                      (ok ? "" : " (" + (res.empty() ? std::string("no result") : res.back().message.substr(0, 160)) + ")"));
        if (!ok) return;
        for (const auto& p2 : pts) {
            got[s].push_back(head_near(res[1].mesh, res[1].pore, p2[0], p2[1]));
            scale = std::fmax(scale, got[s].back());
        }
    }
    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        std::printf("   (%.0f, %.1f): pore %.3f (screen) vs %.3f (strip) kPa\n", pts[i][0], pts[i][1], got[0][i], got[1][i]);
        worst = std::fmax(worst, std::fabs(got[0][i] - got[1][i]));
    }
    check(worst < 0.02 * scale, "the pore pressures of the phase are the strip's within 2%");
    check(got[0][1] - got[0][2] > 3.0, "the screen holds a pore pressure difference (teeth)");
}
}  // namespace

int main() {
    std::printf("A semi-permeable screen in every groundwater flow, against the screen drawn as soil\n\n");
    steady();
    transient();
    staged();
    if (g_failures == 0) {
        std::printf("\nOK: the screen is the strip in steady, transient and staged flow\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", g_failures);
    return 1;
}
