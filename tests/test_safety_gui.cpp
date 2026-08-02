// Safety analysis (phi-c reduction / SRM) through the GUI compute path (build_problem, phase=Safety).
// PLAXIS "Safety" phase -> slope factor of safety + failure mechanism. The core SRM is validated
// directly (test_slope: FoS=1.01 vs ~0.99); here we verify the INTEGRATED GUI path reproduces it and
// returns a non-zero failure mechanism (the slip surface), answering the user's question: a slope DOES
// show collapse -- via the Safety analysis, not the K0 procedure (which is zero-displacement by design).
//
// Rocscience/Slide #1 (Griffiths & Lane): homogeneous 1:2 slope on a foundation, gamma=20.2, c=3 kPa,
// phi=19.6 deg, psi=0 -> FoS ~ 0.99 (Bishop 0.988, Spencer 0.987, Phase2 T6 0.997).
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

m::Project slope() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3; s.gamma_unsat = 20.2; s.c = 3.0; s.phi = 19.6; s.psi = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    // CCW: base, right, top, slope face, foundation top, back.
    P.x = {20, 70, 70, 50, 30, 20};
    P.y = {20, 20, 35, 35, 25, 25};
    P.edge_bc = {(int)m::BCType::FullyFixed,       // base
                 (int)m::BCType::HorizontallyFixed, // right
                 (int)m::BCType::Free,             // top
                 (int)m::BCType::Free,             // slope face
                 (int)m::BCType::Free,             // foundation top
                 (int)m::BCType::HorizontallyFixed};// back
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

void test_safety_le_rejected() {
    // Linear-elastic soil has no strength -> Safety must refuse honestly, not return garbage.
    m::Project pr = slope();
    pr.materials[0].model = m::SoilModel::LinearElastic;
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(!R.ok && R.fos < 0.0, "Safety on linear-elastic soil is rejected (no shear strength)");
}

void test_safety_slope_fos() {
    m::Project pr = slope();
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    check(M.ok, "slope meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(R.ok, "Safety analysis ran");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    const double ref = 0.99;
    const double err = std::fabs(R.fos - ref) / ref * 100.0;
    // Failure mechanism: the slip surface must show NON-zero displacement (unlike the K0 procedure).
    double max_u = R.max_disp;
    std::printf("  GUI Safety: FoS = %.3f  (ref ~0.99, T6 0.997)  err = %.1f%%   mechanism max|u| = %.3e\n",
                R.fos, err, max_u);
    check(err < 8.0, "GUI-path slope factor of safety within 8% of the benchmark");
    check(max_u > 1e-6, "Safety returns a non-zero failure mechanism (slip surface displacement)");
}

void test_safety_hs_gated() {
    // Hardening Soil Safety is GATED honestly: phi-c reduction from a stress-free state is path-
    // unstable with the HS cap (it hangs, or returns a FALSE collapse = a misleading factor of
    // safety). Until a path-stable HS Safety exists, the analysis must REFUSE with guidance (use
    // Mohr-Coulomb c'/phi') rather than return a fake number -- the worst outcome is a wrong FoS the
    // user trusts. (The strength-reduction fix itself -- factoring the HS sub-struct's c/phi -- is in
    // place for when HS Safety is re-enabled.)
    m::Project pr = slope();
    pr.materials[0].model = m::SoilModel::HardeningSoil;
    const auto M = katai::app::mesh_from_project(pr, 25.0, 6);
    check(M.ok, "HS slope meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    std::printf("  HS Safety: ok=%d fos=%.3f  msg=\"%.60s...\"\n", (int)R.ok, R.fos, R.message.c_str());
    check(!R.ok, "HS Safety is refused (gated) -- no misleading factor of safety");
    check(R.fos < 0.0, "HS Safety returns no factor-of-safety number");
}

void test_safety_confined_no_mechanism() {
    // A laterally-confined block under self-weight has NO shear-failure mechanism, so it stays stable
    // up to srf_max. The Safety analysis must report this honestly (fos_lower_bound) instead of
    // presenting the cap (3.0) as a definitive factor of safety.
    m::Project pr;
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3; s.gamma_unsat = 20.0; s.c = 5.0; s.phi = 30.0; s.psi = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 20, 20};            // CCW rectangle
    P.edge_bc = {(int)m::BCType::FullyFixed,        // base
                 (int)m::BCType::HorizontallyFixed, // right (roller)
                 (int)m::BCType::Free,             // top
                 (int)m::BCType::HorizontallyFixed};// left (roller)
    pr.polygons.push_back(P);
    pr.has_water = false;
    const auto M = katai::app::mesh_from_project(pr, 20.0, 6);
    check(M.ok, "confined block meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(R.ok, "Safety ran on the confined block");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    std::printf("  Confined block: FoS = %.3f  lower_bound = %d (expect: no mechanism)\n",
                R.fos, (int)R.fos_lower_bound);
    check(R.fos_lower_bound, "confined block reports NO failure mechanism (FoS is a lower bound, not the cap)");
}

}  // namespace

int main() {
    std::printf("Safety analysis (phi-c reduction) through the GUI compute path\n");
    test_safety_le_rejected();
    test_safety_slope_fos();
    test_safety_confined_no_mechanism();
    test_safety_hs_gated();
    if (g_failures == 0) {
        std::printf("\nOK: GUI Safety analysis -> factor of safety + failure mechanism\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
