// EC7 design approach through the GUI compute path (build_problem, phase design_approach) -- v0.3 B1.3.
// Verifies the design-code material factoring is correctly WIRED into the integrated solve:
//   (1) DesignApproach::None reproduces the plain Safety factor of safety (no-op wiring);
//   (2) EC7 DA3 returns the over-design factor ODF = FoS_characteristic / gamma_M (gamma_M = 1.25) --
//       the exact mesh-free composition identity proven in study_design_ec7, here reached through the
//       FULL GUI path (mesh_from_project -> solve_gravity_le with a phase carrying the design approach).
// PLAXIS applies its Design Approaches by the identical mechanism (reduce c/phi/psi by the factor).
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
namespace app = katai::app;

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
    P.x = {20, 70, 70, 50, 30, 20};
    P.y = {20, 20, 35, 35, 25, 25};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

// Run a Safety analysis through the GUI path with a given design approach; return the reported factor.
double safety_fos(const m::Project& pr, const katai::mesh::Mesh& mesh, m::DesignApproach da) {
    m::Phase ph; ph.type = m::PhaseType::Safety; ph.design_approach = da;
    app::PhaseIO io; io.config = &ph;
    const auto R = app::solve_gravity_le(pr, mesh, app::InitialPhase::Safety, nullptr, io);
    if (!R.ok) { std::printf("  solve failed: %s\n", R.message.c_str()); return -1.0; }
    return R.fos;
}

void test_ec7_da3_gui() {
    const m::Project pr = slope();
    const auto M = app::mesh_from_project(pr, 6.0, 6);   // the ODF/FoS ratio is mesh-independent
    check(M.ok, "slope meshed");
    if (!M.ok) return;

    const double fos_char = safety_fos(pr, M.mesh, m::DesignApproach::None);
    const double odf_da3 = safety_fos(pr, M.mesh, m::DesignApproach::EC7_DA3);
    check(fos_char > 0.0 && odf_da3 > 0.0, "both Safety analyses ran through the GUI path");
    if (fos_char <= 0.0 || odf_da3 <= 0.0) return;

    const double gamma_M = 1.25;  // EC7 M2: gamma_c' = gamma_phi' = 1.25 (drained)
    const double predicted = fos_char / gamma_M;
    const double err = std::fabs(odf_da3 - predicted) / predicted * 100.0;
    std::printf("  GUI  FoS(char) = %.4f   ODF(DA3) = %.4f   FoS/gamma_M = %.4f   err = %.2f%%\n",
                fos_char, odf_da3, predicted, err);
    // The identity is the wiring proof: the design approach reduced c'/phi' inside build_problem.
    check(err < 2.0, "GUI-path EC7 DA3 ODF = FoS_characteristic / gamma_M (factoring correctly wired)");
    // The correct design verdict for these low characteristic strengths (FoS ~ 1.0 < gamma_M).
    check(odf_da3 < 1.0, "GUI-path EC7 DA3 reports this slope as unsafe (ODF < 1)");
}

}  // namespace

int main() {
    std::printf("EC7 design approach through the GUI compute path (build_problem)\n");
    test_ec7_da3_gui();
    if (g_failures == 0) {
        std::printf("\nOK: GUI-path EC7 DA3 material-factoring verified\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
