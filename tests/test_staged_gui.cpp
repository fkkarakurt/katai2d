// Multi-phase staged construction through the GUI compute path (solve_phases / PhaseIO).
// Every phase ramps the configuration imbalance d = f(active) - f_int(committed) -- the PLAXIS
// SumMstage -- so fill placement, excavation unloading and newly activated loads all follow from
// one rule. Closed forms (1D laterally-confined column, E_oed = E(1-nu)/((1+nu)(1-2nu))):
//   FILL:      settlement at the old surface  u = gamma_f h_f H / E_oed, fill top adds h_f^2/(2E_oed);
//   EXCAVATE:  heave at the new surface       u = +gamma h_exc H / E_oed (elastic unloading);
//   NEW LOAD:  next-phase increment           u = q H_total / E_oed.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kE = 1.0e4, kNu = 0.3;
const double kEoed = kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));   // 13461.5

// Two stacked layers: lower 0..6 (gamma 18), upper 6..10 (gamma 17, the fill).
m::Project two_layers() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);
    m::Material f = s; f.gamma_unsat = 17.0;
    pr.materials.push_back(f);
    m::SoilPolygon L; L.material = 0;
    L.x = {0, 20, 20, 0}; L.y = {0, 0, 6, 6};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    m::SoilPolygon U; U.material = 1;
    U.x = {0, 20, 20, 0}; U.y = {6, 6, 10, 10};
    U.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(L);
    pr.polygons.push_back(U);
    pr.has_water = false;
    return pr;
}

// Nearest-node vertical displacement / vertical stress.
double uy_at(const katai::app::SolveResult& R, double x, double y) {
    int best = -1; double bd = 1e300;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        const double d = std::hypot(R.mesh.x[n] - x, R.mesh.y[n] - y);
        if (d < bd) { bd = d; best = n; }
    }
    return R.disp[best * 2 + 1];
}
double sv_base(const katai::app::SolveResult& R) {
    double sv = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] < 1e-6) sv = std::fmin(sv, R.stress.stress[n](1));
    return sv;
}

void test_fill_then_load() {
    std::printf("-- fill placement, then a surcharge in a third phase --\n");
    m::Project pr = two_layers();
    pr.initial.poly_active = {1, 0};            // initial: lower layer only (K0, u = 0)
    m::Phase fill; fill.name = "Fill"; fill.poly_active = {1, 1};
    pr.phases.push_back(fill);
    // Surcharge q = 50 kPa over the full new surface, activated AFTER the fill phase.
    m::Load q; q.kind = m::LoadKind::Distributed; q.name = "q";
    q.x1 = 0; q.y1 = 10; q.x2 = 20; q.y2 = 10; q.qy1 = -50; q.qy2 = -50; q.qx1 = q.qx2 = 0;
    pr.loads.push_back(q);
    pr.initial.load_active = {0};               // not in the initial phase
    pr.phases[0].load_active = {0};             // not in the fill phase either
    m::Phase ld; ld.name = "Surcharge"; ld.poly_active = {1, 1}; ld.load_active = {1};
    pr.phases.push_back(ld);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "two-layer model meshed");
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3, "initial + 2 phases ran");
    if (res.size() != 3) { for (const auto& r : res) std::printf("  (%s)\n", r.message.c_str()); return; }
    for (const auto& r : res) check(r.ok, "phase converged");

    // Initial: undisturbed lower layer -> zero displacement (geostatic K0 on the ACTIVE config).
    std::printf("   initial max|u| = %.3e\n", res[0].max_disp);
    check(res[0].max_disp < 1e-9, "initial phase (lower only): u = 0 exactly");

    // Fill phase: settlement under the new fill weight (gamma_f h_f = 68 kPa).
    const double u6 = -17.0 * 4.0 * 6.0 / kEoed;                       // old surface
    const double u10 = u6 - 17.0 * 4.0 * 4.0 / (2.0 * kEoed);          // fill top (self-compression)
    const double f6 = uy_at(res[1], 10.0, 6.0), f10 = uy_at(res[1], 10.0, 10.0);
    std::printf("   fill: u(y=6)=%.5f (exact %.5f)  u(top)=%.5f (exact %.5f)\n", f6, u6, f10, u10);
    check(std::fabs(f6 - u6) < 0.02 * std::fabs(u6), "old-surface settlement = gamma_f h_f H/E_oed (2%)");
    check(std::fabs(f10 - u10) < 0.02 * std::fabs(u10), "fill-top settlement incl. self-compression (2%)");
    const double sv1 = sv_base(res[1]), sv1_ex = -(18.0 * 6.0 + 17.0 * 4.0);
    std::printf("   fill: base sigma_v=%.2f (exact %.2f)\n", sv1, sv1_ex);
    check(std::fabs(sv1 - sv1_ex) < 0.02 * std::fabs(sv1_ex), "base stress carries the fill weight (2%)");

    // Surcharge phase: INCREMENTAL settlement of the full 10 m column under q = 50.
    const double uq = -50.0 * 10.0 / kEoed;
    const double g10 = uy_at(res[2], 10.0, 10.0);
    std::printf("   surcharge: u(top)=%.5f (exact %.5f)\n", g10, uq);
    check(std::fabs(g10 - uq) < 0.02 * std::fabs(uq), "phase increment = q H/E_oed (2%)");
    const double sv2 = sv_base(res[2]), sv2_ex = sv1_ex - 50.0;
    check(std::fabs(sv2 - sv2_ex) < 0.02 * std::fabs(sv2_ex), "stress chain: base adds the surcharge (2%)");
}

void test_excavation() {
    std::printf("-- excavation (deactivate the upper layer) --\n");
    m::Project pr = two_layers();
    pr.materials[1].gamma_unsat = 17.0;
    m::Phase exc; exc.name = "Excavate"; exc.poly_active = {1, 0};
    pr.phases.push_back(exc);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2, "initial + excavation ran");
    if (res.size() != 2) { for (const auto& r : res) std::printf("  (%s)\n", r.message.c_str()); return; }
    check(res[0].ok && res[1].ok, "both phases converged");
    check(res[0].max_disp < 1e-9, "initial (full geometry, level): u = 0 exactly");

    // Heave of the pit floor: elastic unloading of gamma_up h = 68 kPa over the remaining 6 m.
    const double uh = +17.0 * 4.0 * 6.0 / kEoed;
    const double e6 = uy_at(res[1], 10.0, 6.0);
    std::printf("   excavation: heave u(y=6)=%.5f (exact %.5f)\n", e6, uh);
    check(std::fabs(e6 - uh) < 0.02 * uh, "pit-floor heave = unloading/E_oed (2%)");
    const double sv = sv_base(res[1]), sv_ex = -18.0 * 6.0;
    std::printf("   excavation: base sigma_v=%.2f (exact %.2f)\n", sv, sv_ex);
    check(std::fabs(sv - sv_ex) < 0.02 * std::fabs(sv_ex), "base stress sheds the excavated weight (2%)");
    check(!res[1].active.empty(), "phase reports its element activity (GUI dims the excavation)");
}

}  // namespace

int main() {
    std::printf("Staged construction phases through the GUI compute path\n\n");
    test_fill_then_load();
    test_excavation();
    if (g_failures == 0) {
        std::printf("\nOK: phases = SumMstage chaining (fill settlement, excavation heave, load increments)\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
