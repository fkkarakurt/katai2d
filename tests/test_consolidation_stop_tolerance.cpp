// A consolidation phase that ends on a state -- "until the excess pore pressure is below 1 kPa" --
// answers with a TIME, and that time must be the one the model has, not the one its equilibrium
// tolerance lets it reach early.
//
// The phase's Newton iteration stops when the residual is a set fraction of the forces in play.
// In the tail, where the excess pore pressure creeps under its target, what is left to dissipate is
// a small part of those forces, so an error tolerated against them is a large part of what each
// time step does. At a third of the phase tolerance (the rule for every consolidation phase) the
// tail dissipated early: this embankment stopped four time steps before the time every tighter
// run agrees on. A phase that ends on a state now takes a tenth (driver: where its stop is read).
//
// verify: KV-NUM-015
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, tol, cstop / cminp): a Consolidation phase with a state stop ends at the first time step whose largest excess pore pressure is at or below cminp, and its answer is that time; the tolerated error is relative to the forces in play. The independent path is the same phase solved at a tolerance more than three decades tighter, where the stop time has stopped moving
//   locator:  a strip embankment on soft clay, plane strain, half model: a 30 m wide, 6 m deep Soft Soil clay (lambda* 0.10, kappa* 0.02, nu_ur 0.15, c 1 kPa, phi 25, POP 5 kPa, Undrained (A), k 0.05 m/day, c_k 0.5), drained at its base (head -1 m) and its surface, under two 1.5 m lifts of a Mohr-Coulomb fill (E 15 MPa, phi 30) placed in consolidation phases of 1 day each, then consolidated until the largest excess pore pressure is 1 kPa; tri6, element size 2 m. Solved at the default tolerance and with every phase at 3e-6
//   quantity: the time [day] at which the state-stop phase ends, and the largest excess pore pressure there [kPa]
//   expected: the default tolerance stops at the time the tight run stops at. Before the rule the default stopped at 11.254 day against 11.350 day (four steps of the march early, -0.85%), with 0.976 kPa left where the converged run has 0.993 kPa
//   band:     0.1% on the stop time (measured 0.000%: the same step of the march), 0.01 kPa on the pressure left
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// tol = 0: every phase at its default tolerance.
m::Project embankment(double tol) {
    m::Project pr;
    m::Material soft; soft.name = "Soft"; soft.model = m::SoilModel::SoftSoil;
    soft.drainage = m::Drainage::Undrained;
    soft.lam_star = 0.10; soft.kap_star = 0.02; soft.nu_ur = 0.15;
    soft.c = 1.0; soft.phi = 25.0; soft.psi = 0.0;
    soft.gamma_unsat = 14.0; soft.gamma_sat = 16.0; soft.e_init = 1.5;
    soft.kx = soft.ky = 0.05; soft.ck = 0.5;
    soft.oc_mode = 2; soft.POP = 5.0;
    m::Material fill; fill.name = "Fill"; fill.model = m::SoilModel::MohrCoulomb;
    fill.E = 15e3; fill.nu = 0.3; fill.c = 1.0; fill.phi = 30.0; fill.psi = 0.0;
    fill.gamma_unsat = 16.0; fill.gamma_sat = 19.0; fill.e_init = 0.5; fill.kx = fill.ky = 3.5;
    pr.materials = {soft, fill};

    const int Fr = (int)m::BCType::Free, Hz = (int)m::BCType::HorizontallyFixed,
              Fx = (int)m::BCType::FullyFixed;
    const int Cl = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head,
              Sp = (int)m::FlowBCType::Seepage;
    m::SoilPolygon S; S.name = "Soft"; S.material = 0;
    S.x = {0, 30, 30, 10, 0}; S.y = {-6, -6, 0, 0, 0};
    S.edge_bc = {Fx, Hz, Fr, Fr, Hz};
    S.edge_flow = {Hd, Sp, Sp, Sp, Cl};
    S.edge_head = {-1.0, 0, 0, 0, 0};
    m::SoilPolygon F1; F1.name = "Fill 1"; F1.material = 1;
    F1.x = {0, 10, 7, 0}; F1.y = {0, 0, 1.5, 1.5};
    F1.edge_bc = {Fr, Fr, Fr, Hz}; F1.edge_flow = {Cl, Sp, Sp, Cl}; F1.edge_head = {0, 0, 0, 0};
    m::SoilPolygon F2 = F1; F2.name = "Fill 2";
    F2.x = {0, 7, 4, 0}; F2.y = {1.5, 1.5, 3.0, 3.0};
    pr.polygons = {S, F1, F2};
    pr.has_water = true; pr.wx = {0, 30}; pr.wy = {-1, -1};
    pr.mesh.elem_size = 2.0; pr.mesh.order = 6;
    pr.mesh.corner_refine = false;   // the mesh the locator names (the .k2d v28 one at the fill's slopes)
    pr.initial.poly_active = {1, 0, 0};

    m::Phase p1; p1.name = "Fill 1 (1 d)"; p1.type = m::PhaseType::Consolidation;
    p1.duration = 1.0; p1.time_steps = 10; p1.poly_active = {1, 1, 0}; p1.tolerance = tol;
    m::Phase p2 = p1; p2.name = "Fill 2 (1 d)"; p2.poly_active = {1, 1, 1};
    m::Phase p3 = p2; p3.name = "To 1 kPa";
    p3.consol_stop = m::ConsolStop::MinExcessPore; p3.consol_min_pore = 1.0;
    pr.phases = {p1, p2, p3};
    return pr;
}

struct Stop { bool ok = false; double t = 0.0, p_left = 0.0; int steps = 0; double seconds = 0.0; };

Stop solve(double tol) {
    const m::Project pr = embankment(tol);
    const auto M = katai::app::mesh_from_project(pr);
    Stop s;
    if (!M.ok) { check(false, "the embankment meshes"); return s; }
    const auto t0 = std::chrono::steady_clock::now();
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    s.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (res.size() != 4 || !res[1].ok || !res[2].ok || !res[3].ok) {
        if (!res.empty()) std::printf("   (%s)\n", res.back().message.c_str());
        return s;
    }
    const auto& R = res[3];
    s.ok = R.consol_stop_met && !R.consol_time.empty();
    if (!s.ok) return s;
    s.t = R.consol_time.back();
    s.steps = (int)R.consol_time.size() - 1;
    for (double p : R.excess_pore) s.p_left = std::max(s.p_left, std::fabs(p));
    return s;
}
}  // namespace

int main() {
    std::printf("A consolidation phase that ends on a state: its time at the default tolerance\n\n");
    const Stop def = solve(0.0);
    const Stop tight = solve(3e-6);
    check(def.ok, "the default run reaches its target");
    check(tight.ok, "the tight run reaches its target");
    if (!def.ok || !tight.ok) return 1;
    std::printf("   default  : t = %.5f day, %d steps, %.4f kPa left (%.1f s)\n", def.t, def.steps,
                def.p_left, def.seconds);
    std::printf("   tol 3e-6 : t = %.5f day, %d steps, %.4f kPa left (%.1f s)\n", tight.t, tight.steps,
                tight.p_left, tight.seconds);
    const double dt = def.t / tight.t - 1.0;
    std::printf("   stop time %+.3f%%, pressure left %+.4f kPa\n", 100.0 * dt, def.p_left - tight.p_left);
    check(def.p_left <= 1.0 && tight.p_left <= 1.0, "both stop at or below the 1 kPa they were asked for");
    check(std::fabs(dt) < 1e-3, "the default tolerance stops when the converged run does (0.1%)");
    check(std::fabs(def.p_left - tight.p_left) < 0.01,
          "with the pressure the converged run leaves (0.01 kPa)");
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
