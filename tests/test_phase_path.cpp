// The path of a phase: the curve a factor of safety is read off, and the load-displacement curve
// of a staged construction phase.
//
// A Plastic phase reports, at its start and at every accepted step, the share of the stage applied
// (Sum-Mstage) and the largest displacement of a node since the phase began; a Safety phase the
// strength reduction factor (Sum-Msf) and the same displacement, at every equilibrium its search
// reached. Checked on two models of the corpus:
//   (a) the staged excavation (Mohr-Coulomb ground): the Plastic path starts at (0, 0), its factor
//       rises to 1, and its last displacement is the largest nodal displacement of the phase's own
//       result;
//   (b) the Griffiths-Lane slope, gravity loading and then a Safety phase: the path starts at the
//       full strength (factor 1, no displacement), its factor never falls, it ends at the factor of
//       safety the phase reports, and the displacement it reaches there is many times that of its
//       first steps -- the plateau the factor is read off;
//   (c) both paths survive a results file (v12).
//
// verify: none -- the record of the solver's own steps, checked against the phase's own result; the
// factor of safety itself is verified elsewhere (KV-SLP-*)
#include <katai/io/project_io.hpp>
#include <katai/io/results_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
}  // namespace

int main() {
    std::printf("== The path of a phase ==\n");
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-exc-001-staged-excavation.k2d";
    if (!m::load_project(path, pr, &err)) { std::printf("FATAL: %s\n", err.c_str()); return 1; }
    for (auto& mat : pr.materials) { mat.model = m::SoilModel::MohrCoulomb; mat.c = 10.0; mat.phi = 30.0; }
    const auto M = katai::app::mesh_from_project(pr);
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) {
        std::printf("FATAL: %s\n", res.empty() ? "no result" : res.back().message.c_str());
        return 1;
    }
    m::Project slope;
    const std::string slope_path = std::string(KATAI_CORPUS_DIR) + "/kv-slp-001-griffiths-lane-slope.k2d";
    if (!m::load_project(slope_path, slope, &err)) { std::printf("FATAL: %s\n", err.c_str()); return 1; }
    slope.initial_procedure = m::InitialProcedure::GravityLoading;
    m::Phase safety;
    safety.name = "Safety";
    safety.type = m::PhaseType::Safety;
    slope.phases.push_back(safety);
    const auto MS = katai::app::mesh_from_project(slope);
    const auto sres =
        katai::app::solve_phases(slope, MS.mesh, katai::app::initial_phase_from(slope.initial_procedure));
    if (sres.size() != 2 || !sres[1].ok) {
        std::printf("FATAL: %s\n", sres.empty() ? "no result" : sres.back().message.c_str());
        return 1;
    }

    // (a) the Plastic phase
    const auto& P = res[1];
    double umax = 0.0;
    for (int n = 0; n < (int)(P.disp.size() / 2); ++n)
        umax = std::max(umax, std::hypot(P.disp(2 * n), P.disp(2 * n + 1)));
    bool rising = P.path_factor.size() >= 2;
    for (size_t k = 1; k < P.path_factor.size(); ++k) rising = rising && P.path_factor[k] > P.path_factor[k - 1];
    std::printf("   Plastic: %zu points, factor %.6f -> %.6f, displacement %.6e -> %.6e m (phase max %.6e)\n",
                P.path_factor.size(), P.path_factor.empty() ? 0.0 : P.path_factor.front(),
                P.path_factor.empty() ? 0.0 : P.path_factor.back(), P.path_disp.empty() ? 0.0 : P.path_disp.front(),
                P.path_disp.empty() ? 0.0 : P.path_disp.back(), umax);
    check(P.path_factor.size() == P.path_disp.size() && rising && P.path_factor.front() == 0.0 &&
              P.path_disp.front() == 0.0 && std::fabs(P.path_factor.back() - 1.0) < 1e-12,
          "(a) the Plastic path starts at (0, 0) and its factor rises step by step to 1");
    check(!P.path_disp.empty() && std::fabs(P.path_disp.back() - umax) <= 1e-12 * umax,
          "(a) its last displacement is the phase's largest nodal displacement");

    // (b) the Safety phase
    const auto& S = sres[1];
    bool never_falls = S.path_factor.size() >= 3;
    for (size_t k = 1; k < S.path_factor.size(); ++k) never_falls = never_falls && S.path_factor[k] >= S.path_factor[k - 1];
    const size_t n = S.path_factor.size();
    std::printf("   Safety: %zu points, factor %.6f -> %.6f (FoS %.6f), displacement %.3e (2nd) -> %.3e m\n", n,
                n ? S.path_factor.front() : 0.0, n ? S.path_factor.back() : 0.0, S.fos,
                n > 1 ? S.path_disp[1] : 0.0, n ? S.path_disp.back() : 0.0);
    check(n == S.path_disp.size() && never_falls && n > 0 && S.path_factor.front() == 1.0 && S.path_disp.front() == 0.0,
          "(b) the Safety path starts at full strength with no displacement, and its factor never falls");
    check(n > 0 && S.path_factor.back() <= S.fos * (1.0 + 1e-9) && S.path_factor.back() >= S.fos * 0.99,
          "(b) it ends at the factor of safety the phase reports");
    check(n > 2 && S.path_disp.back() > 10.0 * S.path_disp[1],
          "(b) and its displacement there is many times that of its first step: the plateau");

    // (c) through a results file
    const std::string file = "test_phase_path.res";
    const std::uint64_t hash = 0x1234;
    std::vector<katai::app::SolveResult> back, sback;
    check(katai::app::save_results(file, hash, res, &err) && katai::app::load_results(file, hash, back, &err) &&
              back.size() == 2 && back[1].path_factor == P.path_factor && back[1].path_disp == P.path_disp &&
              katai::app::save_results(file, hash, sres, &err) &&
              katai::app::load_results(file, hash, sback, &err) && sback.size() == 2 &&
              sback[1].path_factor == S.path_factor && sback[1].path_disp == S.path_disp,
          "(c) the paths come back from a results file exactly");
    std::remove(file.c_str());

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
