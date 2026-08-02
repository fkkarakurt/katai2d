// K0 procedure on NON-LEVEL ground through the GUI compute path (build_problem).
//
// PLAXIS Reference Manual: the K0 procedure is only correct when the ground surface, the layer
// boundaries and the water table are ALL horizontal; otherwise the column-overburden field leaves
// genuine unbalanced forces (principal directions must rotate, shear must develop near a slope
// face). PLAXIS resolves this with a "plastic nil-step" whose equilibrium is comparable to gravity
// loading. build_problem now builds that nil-step into the K0 phase: the imbalance
// d = f_body - f_int(seed) is ramped together with the external loads, so the converged state is
// the TRUE equilibrium. On level ground d is round-off and is dropped (u = 0 exactly, the
// undisturbed-case identity) -- including LAYERED and WATER-TABLE columns, which the
// strata-break-exact overburden integral keeps exact.
//
// History note: before the unified-B K0 (commit 3b99321) the GUI ramped full gravity from the K0
// seed, so a slope showed this redistribution as colored |u| bands ("slip-circle rings"); the
// unification silently hid it (residual(0) == 0 by construction). This test pins the restored,
// correct behavior: level => zero, slope => genuine equilibrated redistribution with shear.
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

m::Material mc_soil(double c, double phi) {
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3; s.gamma_unsat = 20.2; s.gamma_sat = 22.0;
    s.c = c; s.phi = phi; s.psi = 0.0;
    return s;
}

// 20x10 level block, one material.
m::Project level_block() {
    m::Project pr;
    pr.materials.push_back(mc_soil(20.0, 25.0));
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

// 20x10 level block, TWO horizontal layers with different gamma/E (interface at y=6).
m::Project level_two_layers() {
    m::Project pr;
    pr.materials.push_back(mc_soil(20.0, 25.0));            // lower, gamma 20.2
    m::Material up = mc_soil(15.0, 22.0); up.gamma_unsat = 17.0; up.E = 4.0e4;
    pr.materials.push_back(up);                              // upper, gamma 17
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
    return pr;
}

// Griffiths & Lane 1:2 slope geometry (as test_safety_gui), strength c/phi chosen per test.
m::Project slope(double c, double phi) {
    m::Project pr;
    pr.materials.push_back(mc_soil(c, phi));
    m::SoilPolygon P; P.material = 0;
    P.x = {20, 70, 70, 50, 30, 20};
    P.y = {20, 20, 35, 35, 25, 25};
    P.edge_bc = {(int)m::BCType::FullyFixed,        // base
                 (int)m::BCType::HorizontallyFixed, // right
                 (int)m::BCType::Free,              // top
                 (int)m::BCType::Free,              // slope face
                 (int)m::BCType::Free,              // foundation top
                 (int)m::BCType::HorizontallyFixed};// back
    pr.polygons.push_back(P);
    return pr;
}

void test_level_exact_zero() {
    std::printf("-- level ground: K0 stays exactly zero-displacement (no nil-step) --\n");
    {   // single layer, dry
        const auto pr = level_block();
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
        std::printf("   single layer: max|u|=%.3e  nil_step=%d\n", R.max_disp, (int)R.nil_step);
        check(R.ok && !R.nil_step && R.max_disp < 1e-9, "single layer level: u = 0, no nil-step");
    }
    {   // two layers (gamma jump at y=6) -- strata-break-exact integral must keep the identity
        const auto pr = level_two_layers();
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
        std::printf("   two layers:   max|u|=%.3e  nil_step=%d\n", R.max_disp, (int)R.nil_step);
        check(R.ok && !R.nil_step && R.max_disp < 1e-9, "two-layer level: u = 0, no nil-step");
    }
    {   // water table at y=8 (horizontal) -- buoyancy switch must not fake an imbalance
        auto pr = level_block();
        pr.has_water = true; pr.wx = {0.0, 20.0}; pr.wy = {8.0, 8.0};
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
        std::printf("   water y=8:    max|u|=%.3e  nil_step=%d\n", R.max_disp, (int)R.nil_step);
        check(R.ok && !R.nil_step && R.max_disp < 1e-9, "level + horizontal water table: u = 0, no nil-step");
    }
    {   // SLOPED water table on level ground -> non-level geometry -> nil-step must engage
        auto pr = level_block();
        pr.has_water = true; pr.wx = {0.0, 20.0}; pr.wy = {9.0, 5.0};
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
        std::printf("   sloped water: max|u|=%.3e  nil_step=%d\n", R.max_disp, (int)R.nil_step);
        check(R.ok && R.nil_step, "sloped water table triggers the equilibrium nil-step");
    }
}

void test_slope_nil_step() {
    std::printf("-- slope: K0 resolves the genuine imbalance (PLAXIS plastic nil-step) --\n");
    const auto pr = slope(20.0, 25.0);   // stable slope (well above FoS 1)
    const auto M = katai::app::mesh_from_project(pr, 2.5, 6);
    check(M.ok, "slope meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "K0 on the slope converged (nil-step equilibrium reached)");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.nil_step, "the non-level surface was detected (nil-step active)");
    std::printf("   slope K0: max|u|=%.4e  nil_step=%d\n", R.max_disp, (int)R.nil_step);
    check(R.max_disp > 1e-5, "redistribution displacement is genuinely non-zero");
    // The raw K0 field has sigma_xy == 0 everywhere; equilibrium on a slope REQUIRES shear.
    double max_sxy = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        max_sxy = std::fmax(max_sxy, std::fabs(R.stress.stress[n](2)));
    std::printf("   slope K0: max|sigma_xy| = %.2f kPa\n", max_sxy);
    check(max_sxy > 5.0, "shear stress developed near the slope face (> 5 kPa)");

    // PLAXIS: the nil-step equilibrium is comparable to gravity loading. Anchor: sigma'_v at a
    // deep point under the level crest (x=60, y=22, overburden -gamma*(35-22) = -262.6) matches
    // between the two procedures and the overburden estimate.
    const auto G = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    check(G.ok, "gravity loading on the slope converged");
    if (G.ok) {
        const double sv_ref = -20.2 * (35.0 - 22.0);
        auto sv_at = [](const katai::app::SolveResult& S, double x, double y) {
            int best = -1; double bd = 1e300;
            for (int n = 0; n < S.mesh.node_count; ++n) {
                const double d = std::hypot(S.mesh.x[n] - x, S.mesh.y[n] - y);
                if (d < bd) { bd = d; best = n; }
            }
            return S.stress.stress[best](1);
        };
        const double sk = sv_at(R, 60.0, 22.0), sg = sv_at(G, 60.0, 22.0);
        std::printf("   crest sigma'_v: K0+nil=%.1f  gravity=%.1f  overburden=%.1f kPa\n", sk, sg, sv_ref);
        check(std::fabs(sk - sv_ref) < 0.10 * std::fabs(sv_ref), "K0+nil sigma'_v ~ overburden (10%)");
        check(std::fabs(sk - sg) < 0.10 * std::fabs(sv_ref), "K0+nil ~ gravity loading (PLAXIS claim)");
    }
}

void test_unstable_slope_reported() {
    std::printf("-- unstable slope: K0 nil-step reports collapse honestly --\n");
    // c = 1 kPa on the Griffiths & Lane slope -> FoS well below 1; gravity equilibrium does not
    // exist, so the nil-step must NOT silently return the (non-equilibrium) K0 field.
    const auto pr = slope(1.0, 19.6);
    const auto M = katai::app::mesh_from_project(pr, 2.5, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    std::printf("   c=1 slope: ok=%d load_factor=%.2f\n", (int)R.ok, R.load_factor);
    check(!R.ok && R.load_factor < 1.0, "collapsing slope is reported (no fake equilibrium)");
}

}  // namespace

int main() {
    std::printf("K0 procedure on non-level ground (nil-step) through the GUI path\n\n");
    test_level_exact_zero();
    test_slope_nil_step();
    test_unstable_slope_reported();
    if (g_failures == 0) {
        std::printf("\nOK: K0 = exact zero on level ground, true equilibrium (nil-step) on slopes\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
