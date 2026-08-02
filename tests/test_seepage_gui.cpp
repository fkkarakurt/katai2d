// Groundwater flow through the GUI compute path (build_flow + build_problem coupling).
// The seepage CORE is validated quantitatively in test_seepage (Darcy, MMS, Thiem, Charny,
// Terzaghi gradient coupling); here the INTEGRATED path -- project edge flow-BCs -> mesh ->
// solve_groundwater_flow -> (optionally) solve_gravity_le(flow_head) -- must reproduce the
// same closed forms:
//   (1) confined 1D Darcy between two reservoirs: h linear, Q = k i A (exact);
//   (2) unconfined dam with tailwater + seepage face: Charny q = k(h1^2-h2^2)/(2L);
//   (3) upward-seepage column coupled into effective stress: sigma'_v = -(gamma'-i gamma_w)(H-y).
#include <katai/jobs/flow_driver.hpp>
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

constexpr int kClosed = (int)m::FlowBCType::Closed;
constexpr int kHead = (int)m::FlowBCType::Head;
constexpr int kSeep = (int)m::FlowBCType::Seepage;

// (1) Confined flow between two reservoirs across a 20 x 10 block (both heads above the top,
// so the domain stays fully saturated): h(x) = 15 - 3x/20 is linear -> the FE solution is exact;
// discharge Q = k * i * A = 1 * (3/20) * 10 = 1.5 m3/day/m.
void test_confined_darcy() {
    std::printf("-- (1) confined Darcy reservoir-to-reservoir --\n");
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.kx = 1.0; s.ky = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {kClosed, kHead, kClosed, kHead};   // right = 12 m, left = 15 m
    P.edge_head = {0.0, 12.0, 0.0, 15.0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "block meshed");
    const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(F.ok, "confined flow solved");
    if (!F.ok) { std::printf("  (%s)\n", F.message.c_str()); return; }

    double max_h_err = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        const double h_exact = 15.0 - 3.0 * M.mesh.x[n] / 20.0;
        max_h_err = std::fmax(max_h_err, std::fabs(F.head[n] - h_exact));
    }
    const double q_exact = 1.0 * (3.0 / 20.0) * 10.0;
    std::printf("   max|h-h_exact|=%.3e  Q=%.6f (exact %.6f)  balance=%.2e  iters=%d\n",
                max_h_err, F.discharge, q_exact, F.balance_err, F.iterations);
    check(max_h_err < 1e-6, "linear head field reproduced exactly");
    check(std::fabs(F.discharge - q_exact) < 1e-5 * q_exact, "discharge Q = k i A exact");
    check(F.balance_err < 1e-8, "global mass balance ~ 0");
}

// (2) Unconfined dam between reservoirs h1=5 / tailwater h2=1 with a SEEPAGE FACE above the
// tailwater (the canonical free-surface benchmark). Charny's theorem: the discharge is exact,
// q = k (h1^2 - h2^2) / (2L), regardless of the free-surface shape. The polygon carries the
// BCs the way a user would draw them: extra vertices at the reservoir / tailwater levels split
// the faces, so each edge takes one condition.
void test_unconfined_charny() {
    std::printf("-- (2) unconfined dam + seepage face (Charny discharge) --\n");
    constexpr double L = 10.0, D = 6.0, h1 = 5.0, h2 = 1.0, k = 0.5;
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.kx = k; s.ky = k;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    // CCW from origin; vertices at (L, h2) and (0, h1) split the vertical faces.
    P.x = {0, L, L, L, 0, 0};
    P.y = {0, 0, h2, D, D, h1};
    P.edge_bc.assign(6, (int)m::BCType::Free);
    P.edge_bc[0] = (int)m::BCType::FullyFixed;
    //            bottom   right<h2   right>h2  top      left>h1   left<h1
    P.edge_flow = {kClosed, kHead,     kSeep,    kClosed, kClosed,  kHead};
    P.edge_head = {0.0,     h2,        0.0,      0.0,     0.0,      h1};
    pr.polygons.push_back(P);
    pr.has_water = false;

    const auto M = katai::app::mesh_from_project(pr, 0.06, 6);   // ~0.35 m elements
    check(M.ok, "dam meshed");
    const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(F.ok, "unconfined flow with seepage face solved");
    if (!F.ok) { std::printf("  (%s)\n", F.message.c_str()); return; }

    const double q_exact = k * (h1 * h1 - h2 * h2) / (2.0 * L);
    const double relerr = std::fabs(F.discharge - q_exact) / q_exact;
    std::printf("   q=%.5f  q_charny=%.5f  relerr=%.2e  iters=%d  balance=%.2e\n",
                F.discharge, q_exact, relerr, F.iterations, F.balance_err);
    check(relerr < 0.03, "discharge matches Charny q = k(h1^2-h2^2)/(2L) within 3%");

    // Phreatic surface (pore = 0) must start at the reservoir level and stay inside the dam.
    int n_top_sat = 0;
    for (int n = 0; n < M.mesh.node_count; ++n)
        if (M.mesh.y[n] > D - 1e-6 && F.pore[n] > 1.0) ++n_top_sat;
    check(n_top_sat == 0, "crest stays unsaturated (free surface below the top)");
}

// (3) Coupling into the mechanical solve: upward seepage in a saturated column reduces the
// effective stress, sigma'_v(y) = -(gamma' - i gamma_w)(H - y) (Terzaghi). Flow: bottom head
// 12 m, top head 10 m (= surface) -> i = 0.2 upward, h linear (exact). The K0 phase under flow
// coupling always runs its equilibrium nil-step, so the converged state is the true equilibrium
// with the seepage pore field; in a 1D column sigma'_v is statically determinate, so the base
// value must match the closed form.
void test_coupled_upward_seepage() {
    std::printf("-- (3) flow -> effective stress coupling (upward seepage column) --\n");
    constexpr double H = 10.0, W = 4.0, dh = 2.0, gamma_sat = 20.0;
    const double gw = katai::app::kGammaWater;
    const double i = dh / H;
    const double sv_base = -((gamma_sat - gw) - i * gw) * H;   // -(gamma' - i gamma_w) H

    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = gamma_sat;
    s.kx = 1.0; s.ky = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0};
    P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {kHead, kClosed, kHead, kClosed};   // bottom = 12 m, top = 10 m
    P.edge_head = {H + dh, 0.0, H, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;   // saturation comes from the flow head field, not a polyline

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "column meshed");
    const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(F.ok, "column flow solved");
    if (!F.ok) { std::printf("  (%s)\n", F.message.c_str()); return; }

    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure, &F.head);
    check(R.ok, "coupled K0 solve converged");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.nil_step, "flow coupling runs the equilibrium nil-step");

    double base_sv = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] < 1e-6) base_sv = std::fmin(base_sv, R.stress.stress[n](1));
    std::printf("   base sigma'_v=%.3f  exact=%.3f  (i=%.2f upward)\n", base_sv, sv_base, i);
    check(std::fabs(base_sv - sv_base) < 0.01 * std::fabs(sv_base),
          "effective stress reduced by upward seepage (Terzaghi, <1%)");

    // Pore field for display must come from the flow head, not a (missing) hydrostatic polyline.
    double base_pore = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] < 1e-6) base_pore = std::fmax(base_pore, R.pore[n]);
    const double pore_exact = gw * (H + dh);
    std::printf("   base pore=%.3f  exact=%.3f\n", base_pore, pore_exact);
    check(std::fabs(base_pore - pore_exact) < 0.01 * pore_exact, "pore display from the flow head");
}

}  // namespace

int main() {
    std::printf("Groundwater flow through the GUI compute path\n\n");
    test_confined_darcy();
    test_unconfined_charny();
    test_coupled_upward_seepage();
    if (g_failures == 0) {
        std::printf("\nOK: GUI flow path = Darcy exact, Charny discharge, Terzaghi coupling\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
