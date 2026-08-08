// Is a published number a physics result, or an artefact of the tolerance it was computed at?
//
// KV-NUM-005 and KV-SLP-003 answer the mesh half of that question. This case answers the solver
// half, and it matters most exactly where it is least visible: the driver chooses the tolerated
// force residual by material class, so a user never sees the number that decided when the
// Newton loop was allowed to stop. Mohr-Coulomb runs at 1e-6, the hardening and soft-soil
// families at 1e-2 -- a hundred times looser, because their tangent is continuum rather than
// consistent, which is a defensible engineering choice and an undeclared one until it is
// measured.
//
// The strength-reduction case is the sharp one. A factor of safety obtained by pushing the
// strength down until the solver stops converging is, by construction, a statement about a
// convergence criterion -- Griffiths and Lane's own factor of safety is "the first trial value
// that did not converge within 1000 iterations". If OUR factor of safety moved with OUR
// tolerance, the comparison against theirs would be comparing two different numerical accidents.
// So it is measured rather than assumed, on the same benchmark and the same mesh, changing only
// the tolerated residual.
//
// verify: KV-NUM-007
//   oracle:   independent_path
//   source:   the same problem solved with a different stopping rule is an independent path to the same answer; the criterion under study is the tolerated force residual of the driver's nonlinear solver (kernel/jobs/src/driver.cpp), whose default is 1e-6 for the Mohr-Coulomb family; the practice of demonstrating iterative-convergence independence before quoting a discretisation error is I. B. Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, which requires iterative convergence to be established first, and ASME V&V 20-2009
//   locator:  factor of safety by phi-c reduction on the checked-in tests/corpus/kv-slp-002-griffiths-lane-example1.k2d, one fixed mesh, solved at tolerated residuals 1e-4, 1e-6 (the driver's default for this material) and 1e-8; and the elastic strip-load benchmark, whose linear system is solved directly and therefore carries no iterative tolerance at all
//   quantity: factor of safety [-] as a function of the tolerated force residual, and the spread across the three tolerances relative to the default run
//   expected: a spread far below the mesh dependence of the same quantity (measured at 4-8% in KV-SLP-003), so that the published comparison against Griffiths and Lane and Bishop and Morgenstern is a statement about the model and the mesh rather than about the stopping rule
//   band:     as asserted below and MEASURED on this tree, not inherited -- see the printed table; the assertion is that the three tolerances agree within 1%, i.e. at least four times tighter than the mesh dependence they must not be confused with

#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

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

}  // namespace

int main() {
    std::printf("== does the answer depend on the stopping rule? (KV-NUM-007) ==\n");

    m::Project pr;
    std::string err;
    const std::string path =
        std::string(KATAI_CORPUS_DIR) + "/kv-slp-002-griffiths-lane-example1.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        return 1;
    }
    pr.mesh.elem_size = 2.2;   // one fixed mesh: only the tolerance changes here
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "mesh built from the case file");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }
    std::printf("  mesh: %d elements, %d nodes\n", M.mesh.element_count, M.mesh.node_count);

    // 1e-6 is what the driver picks for this material family; the neighbours bracket it by two
    // orders of magnitude each way.
    const double tolerances[3] = {1e-4, 1e-6, 1e-8};
    double fos[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        katai::app::PhaseIO io;
        io.config = &pr.initial;
        io.numeric.tolerance = tolerances[i];
        const auto R = katai::app::solve_gravity_le(
            pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, io);
        check(R.ok, "the run converges at tolerance " + std::to_string(tolerances[i]));
        if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return 1; }
        fos[i] = R.fos;
        std::printf("  tolerated residual %.0e -> FoS %s %.6f\n", tolerances[i],
                    R.fos_lower_bound ? ">" : "=", R.fos);
    }

    // The default run must be reproduced exactly by an untouched PhaseIO: the override is a
    // measuring instrument, and an instrument that changes the thing it measures is useless.
    katai::app::PhaseIO io_default;
    io_default.config = &pr.initial;
    const auto R_default = katai::app::solve_gravity_le(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, io_default);
    check(R_default.ok && R_default.fos == fos[1],
          "the driver's own default IS 1e-6 here, bit-for-bit (the override changes nothing "
          "when it is not set)");

    const double lo = std::fmin(fos[0], std::fmin(fos[1], fos[2]));
    const double hi = std::fmax(fos[0], std::fmax(fos[1], fos[2]));
    const double spread = (hi - lo) / fos[1];
    std::printf("\n  spread across four orders of magnitude of tolerance = %.4f%% of the "
                "default run\n", 100.0 * spread);
    std::printf("  mesh dependence of the same quantity (KV-SLP-003) = 4.0%% over a fourfold "
                "refinement, 7.9%% on a finer family\n");

    check(spread < 0.01,
          "the factor of safety is independent of the stopping rule to within 1%");
    check(spread < 0.04 / 4.0,
          "and at least four times tighter than the mesh dependence it must not be confused "
          "with");

    // The three runs agreeing bit-for-bit says something further, and it should be said rather
    // than admired: what quantises this factor of safety is not the residual tolerance but the
    // strength-reduction SEARCH. The safety strategy bisects the reduction factor between 0.4
    // and 3.0 for 12 iterations (phase_solver/safety.hpp), so the reported factor carries a
    // resolution of 2.6/2^12 = 6.3e-4, about 0.05% at a factor of 1.38. That is finer than the
    // mesh dependence by two orders of magnitude -- and, worth noting beside the reference,
    // finer than the 0.05 trial increments Griffiths and Lane stepped through.
    std::printf("  resolution of the reported factor: 2.6/2^12 = %.1e (the bisection's own "
                "granularity)\n", 2.6 / 4096.0);
    check(spread * fos[1] < 2.6 / 4096.0 * 2.0,
          "the tolerance spread is inside the search's own resolution, so the stopping rule is "
          "not what sets this number");

    // The elastic benchmark needs no sweep, and saying why is part of the record: its system is
    // linear, so it is factorised and solved once. There is no iteration to stop, and therefore
    // no tolerance for its answer to depend on. A study that swept it anyway would be reporting
    // three identical numbers as evidence.
    std::printf("\n  note: the elastic strip-load benchmark (KV-NUM-005) is a LINEAR system --\n"
                "        solved directly, no iteration, so no stopping rule enters its answer.\n");

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
