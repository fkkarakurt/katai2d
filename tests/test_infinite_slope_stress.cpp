// test_infinite_slope_stress -- the initial stresses of ground on a slope, against the infinite slope.
//
// Every slope model starts here: the stresses gravity gives to ground whose surface is inclined.
// On level ground the K0 procedure states them; on a slope it cannot (it has no equilibrium to
// offer), and the initial phase is gravity loading -- either chosen, or taken when the K0
// procedure is asked for and fails on the slope (K2D-A022). Until 2026-10-22 that path had no
// independent check: the K0 cases are level ground, and the slope cases read only a factor of
// safety.
//
// Oracle (closed form, equilibrium alone): in an infinite slope of inclination beta the traction on
// a plane parallel to the surface at vertical depth z is
//     sigma_n = gamma z cos^2(beta)        tau = gamma z sin(beta) cos(beta),
// whatever the elastic constants, because the column between two vertical sections carries its
// own weight and nothing else crosses the sections (they carry the same stress on both sides).
// The model is a layer of uniform vertical thickness on a fixed base, twenty times as long as it is
// deep; the traction is read in its middle third, away from the ends.
//
// verify: KV-SLP-004
//   oracle:   closed_form
//   source:   the infinite slope (Taylor, D.W. (1948). Fundamentals of Soil Mechanics. Wiley, New York, ch. 16): on planes parallel to the surface sigma_n = gamma z cos^2 beta and tau = gamma z sin beta cos beta, from the equilibrium of a vertical slice whose two sides carry equal and opposite forces
//   locator:  the traction of the recovered nodal stress on the slope-parallel plane, evaluated at every node of the middle third of the model, against the two expressions above at that node's vertical depth below the surface
//   quantity: normal and shear traction on slope-parallel planes [kPa], by gravity loading and by the K0 procedure (whose equilibrium step on sloping ground, or its fall-back to gravity loading, K2D-A022, must reach the same state)
//   expected: both tractions within 1% of gamma z cos^2 beta and gamma z sin beta cos beta (relative to gamma z cos beta) at every node of the middle third deeper than 1 m, for slopes of 10 and 25 degrees, on 6- and 15-noded elements
//   band:     1% -- the ends of a finite model and the nodal recovery of a Gauss-point field; the answer has no material dependence to approximate. Measured: at most 0.012% normal and 0.028% shear, on every combination

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

struct Worst { bool ok = false; double n = 0.0, s = 0.0; int nodes = 0; std::string msg; bool fallback = false; };

static Worst run(double beta_deg, int order, InitialPhase proc) {
    constexpr double PI = 3.14159265358979323846;
    const double beta = beta_deg * PI / 180.0, t = std::tan(beta);
    constexpr double H = 5.0, L = 100.0, gamma = 18.0;
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic; s.E = 3.0e4; s.nu = 0.3;
    s.gamma_unsat = gamma; s.gamma_sat = gamma; s.e_init = 0.5;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    // base from (0,0) to (L, L tan b); surface H above it, parallel
    P.x = {0.0, L, L, 0.0};
    P.y = {0.0, L * t, L * t + H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    pr.initial_procedure = proc == InitialPhase::GravityLoading ? m::InitialProcedure::GravityLoading
                                                                : m::InitialProcedure::K0Procedure;
    Worst w;
    const auto M = katai::app::mesh_from_project(pr, 1.0, order);
    if (!M.ok) { w.msg = M.message; return w; }
    const auto res = katai::app::solve_phases(pr, M.mesh, proc);
    if (res.empty() || !res[0].ok) { w.msg = res.empty() ? "no result" : res[0].message; return w; }
    const auto& R = res[0];
    for (const auto& d : R.diagnostics) w.fallback = w.fallback || d.code == "K2D-A022";
    w.ok = true;
    // slope-parallel plane: unit normal n = (-sin b, cos b), along-slope s = (cos b, sin b)
    const double nx = -std::sin(beta), ny = std::cos(beta), sx = std::cos(beta), sy = std::sin(beta);
    for (int k = 0; k < R.mesh.node_count; ++k) {
        const double x = R.mesh.x[k], y = R.mesh.y[k];
        if (x < L / 3.0 || x > 2.0 * L / 3.0) continue;
        const double z = (x * t + H) - y;   // vertical depth below the surface
        if (z < 1.0) continue;
        const auto& sg = R.stress.stress[(size_t)k];   // [sxx, syy, sxy], tension positive
        const double tx = sg(0) * nx + sg(2) * ny, ty = sg(2) * nx + sg(1) * ny;
        const double tn = -(tx * nx + ty * ny);          // compression positive
        const double ts = tx * sx + ty * sy;             // along the slope, uphill positive
        const double ref = gamma * z * std::cos(beta);
        const double want_n = gamma * z * std::cos(beta) * std::cos(beta);
        const double want_s = gamma * z * std::sin(beta) * std::cos(beta);
        w.n = std::max(w.n, std::fabs(tn - want_n) / ref);
        w.s = std::max(w.s, std::fabs(std::fabs(ts) - want_s) / ref);
        ++w.nodes;
    }
    return w;
}

int main() {
    std::printf("Initial stresses on a slope against the infinite slope\n\n");
    for (double beta : {10.0, 25.0})
        for (int order : {6, 15})
            for (InitialPhase proc : {InitialPhase::GravityLoading, InitialPhase::K0Procedure}) {
                const Worst w = run(beta, order, proc);
                const char* how = proc == InitialPhase::GravityLoading ? "gravity loading" : "K0 procedure";
                char label[200];
                std::snprintf(label, sizeof(label), "beta %.0f deg, %d-noded, %s", beta, order, how);
                if (!w.ok) { check(false, std::string(label) + ": solved (" + w.msg + ")"); continue; }
                std::printf("  %s: worst normal %.3f%%, shear %.3f%% over %d nodes%s\n", label, 100.0 * w.n,
                            100.0 * w.s, w.nodes, w.fallback ? " (fell back to gravity loading, K2D-A022)" : "");
                check(w.nodes > 50 && w.n < 0.01 && w.s < 0.01,
                      std::string(label) + ": the slope-parallel tractions are the infinite slope's within 1%");
            }
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
