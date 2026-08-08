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
// 2026-08-09: THIS CASE ONCE PROVED NOTHING, AND THE RECORD SAYS SO. The strength-reduction
// strategy hard-coded its trial tolerance (phase_solver/safety.hpp: NewtonOptions{8, 120, 1e-3}),
// so the sweep below set a number the search never read. Three identical factors of safety came
// back and were reported as independence -- when what they showed was that nothing had changed.
// A sweep over an input that is silently dropped is indistinguishable, in its output, from a
// genuine insensitivity; that is what makes this failure worth recording rather than quietly
// fixing. The controls are now threaded into the search (and, since .k2d v7, into the file), and
// the re-measured answer is BETTER than the claim it replaces AND has a hard edge the original
// could never have found: below the search's own 1e-3 the factor is bit-identical over five
// orders of magnitude, while ABOVE it the factor climbs, monotonically and always UNSAFE-SIDED
// (+2.0% at 1e-2, +45.6% at 1e-1, and at 3e-1 the search returns its own cap because no trial
// ever fails to "converge"). A slope that is reported 45% safer than it is would be a fine
// example of the kind of number this project exists to refuse.
//
// verify: KV-NUM-007
//   oracle:   independent_path
//   source:   the same problem solved with a different stopping rule is an independent path to the same answer; the criterion under study is the tolerated force residual of the driver's nonlinear solver (kernel/jobs/src/driver.cpp), whose default is 1e-6 for the Mohr-Coulomb family; the practice of demonstrating iterative-convergence independence before quoting a discretisation error is I. B. Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, which requires iterative convergence to be established first, and ASME V&V 20-2009
//   locator:  factor of safety by phi-c reduction on the checked-in tests/corpus/kv-slp-002-griffiths-lane-example1.k2d, one fixed mesh, solved at tolerated residuals 1e-1, 1e-2, 1e-3 (what the strength-reduction search uses when nothing is asked for), 1e-4, 1e-6 and 1e-8; and the elastic strip-load benchmark, whose linear system is solved directly and therefore carries no iterative tolerance at all
//   quantity: factor of safety [-] and, for the Hardening Soil oedometer KV-CST-002, the settlement of the loading step [m], each as a function of the tolerated force residual, with the spread relative to the default run and the SIGN of the deviation on the loose side
//   expected: below the search's own stopping rule, a spread far under the mesh dependence of the same quantity (4-8% in KV-SLP-003), so that the published comparison is a statement about the model and the mesh rather than about the stopping rule; above it, a monotone one-sided error -- a looser rule must report a HIGHER factor of safety, since a trial that stops early counts as equilibrium; for the Hardening Soil family, whose default residual is a hundred times looser, a bounded and stated cost rather than an unknown one
//   band:     as asserted below and MEASURED on this tree, not inherited. Slope factor of safety on this mesh: 1.421021 at 1e-3, 1e-4, 1e-6 and 1e-8 -- BIT-IDENTICAL, spread 0.0000%, inside the reduction search's own resolution of 6.3e-4; 1.449585 at 1e-2 (+2.0%) and 2.069116 at 1e-1 (+45.6%), monotone and always high. Hardening Soil oedometer: 0.018934 m at the default 1e-2 (+0.413% vs the closed form), 0.019045 m at 1e-4 (+1.002%), 0.019044 m at 1e-6 -- so the default costs 0.59% on this quantity and 1e-4 is already converged. Note the default's smaller deviation is cancellation, not accuracy. Asserted: the tight-side spread below 1% and inside the search resolution, monotone unsafe-sided inflation on the loose side with 1e-1 above +30%, 1e-4 and 1e-6 agreeing to 0.05%, and the default HS run within 3% of the closed form

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

    // The strength-reduction search runs its trial solves at 1e-3; the sweep brackets it two
    // orders looser and five orders tighter. Both directions are needed, and the loose side is
    // the one nobody measures: a trial that stops early looks converged, so the slope appears to
    // survive a strength reduction it cannot actually sustain.
    constexpr int kN = 6;
    const double tolerances[kN] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8};
    constexpr int kDefault = 2;   // 1e-3: what the search uses when nothing is asked for
    double fos[kN] = {};
    for (int i = 0; i < kN; ++i) {
        katai::app::PhaseIO io;
        io.config = &pr.initial;
        io.numeric.tolerance = tolerances[i];
        const auto R = katai::app::solve_gravity_le(
            pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, io);
        check(R.ok, "the run converges at tolerance " + std::to_string(tolerances[i]));
        if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return 1; }
        fos[i] = R.fos;
        std::printf("  tolerated residual %.0e -> FoS %s %.6f   (%+.2f%% vs the default rule)\n",
                    tolerances[i], R.fos_lower_bound ? ">" : "=", R.fos,
                    100.0 * (R.fos - 1.421020508) / 1.421020508);
    }

    // The default run must be reproduced exactly by an untouched PhaseIO: the override is a
    // measuring instrument, and an instrument that changes the thing it measures is useless.
    katai::app::PhaseIO io_default;
    io_default.config = &pr.initial;
    const auto R_default = katai::app::solve_gravity_le(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, io_default);
    check(R_default.ok && R_default.fos == fos[kDefault],
          "an untouched PhaseIO reproduces the 1e-3 run bit-for-bit (the override changes "
          "nothing when it is not set)");

    // THE LOOSE SIDE, and it is one-sided. Every relaxation of the stopping rule RAISES the
    // reported factor of safety, because the question the search asks -- did this reduced
    // strength still reach equilibrium? -- is answered "yes" by a solver that was allowed to
    // stop before it got there. The error therefore always points the unsafe way, and it is
    // large long before the tolerance looks absurd.
    check(fos[0] > fos[1] && fos[1] > fos[2],
          "a looser stopping rule reports a HIGHER factor of safety, monotonically");
    std::printf("\n  1e-1 inflates the factor of safety by %+.1f%%, 1e-2 by %+.1f%% -- always "
                "unsafe-sided\n",
                100.0 * (fos[0] - fos[kDefault]) / fos[kDefault],
                100.0 * (fos[1] - fos[kDefault]) / fos[kDefault]);
    check((fos[0] - fos[kDefault]) / fos[kDefault] > 0.30,
          "and at 1e-1 it is not a rounding matter: over 30% too high");

    // THE TIGHT SIDE, which is what the published comparison rests on: below the search's own
    // rule the answer stops moving entirely.
    const double lo = std::fmin(fos[kDefault], std::fmin(fos[3], std::fmin(fos[4], fos[5])));
    const double hi = std::fmax(fos[kDefault], std::fmax(fos[3], std::fmax(fos[4], fos[5])));
    const double spread = (hi - lo) / fos[kDefault];
    std::printf("\n  spread from the default rule down through five orders of magnitude = "
                "%.4f%% of the default run\n", 100.0 * spread);
    std::printf("  mesh dependence of the same quantity (KV-SLP-003) = 4.0%% over a fourfold "
                "refinement, 7.9%% on a finer family\n");

    check(spread < 0.01,
          "at and below the default rule the factor of safety is independent of it to within 1%");
    check(spread < 0.04 / 4.0,
          "and at least four times tighter than the mesh dependence it must not be confused "
          "with");

    // The four runs agreeing bit-for-bit says something further, and it should be said rather
    // than admired: what quantises this factor of safety is not the residual tolerance but the
    // strength-reduction SEARCH. The safety strategy bisects the reduction factor between 0.4
    // and 3.0 for 12 iterations (phase_solver/safety.hpp), so the reported factor carries a
    // resolution of 2.6/2^12 = 6.3e-4, about 0.05% at a factor of 1.42. That is finer than the
    // mesh dependence by two orders of magnitude -- and, worth noting beside the reference,
    // finer than the 0.05 trial increments Griffiths and Lane stepped through.
    std::printf("  resolution of the reported factor: 2.6/2^12 = %.1e (the bisection's own "
                "granularity)\n", 2.6 / 4096.0);
    check(spread * fos[kDefault] < 2.6 / 4096.0 * 2.0,
          "the remaining spread is inside the search's own resolution, so below the default the "
          "stopping rule is not what sets this number");

    // The elastic benchmark needs no sweep, and saying why is part of the record: its system is
    // linear, so it is factorised and solved once. There is no iteration to stop, and therefore
    // no tolerance for its answer to depend on. A study that swept it anyway would be reporting
    // three identical numbers as evidence.
    std::printf("\n  note: the elastic strip-load benchmark (KV-NUM-005) is a LINEAR system --\n"
                "        solved directly, no iteration, so no stopping rule enters its answer.\n");

    // ---------------------------------------------------------------------------------------
    // The case this study was really needed for. The driver gives the hardening and soft-soil
    // families a tolerated residual of 1e-2 -- a hundred times looser than Mohr-Coulomb -- on
    // the grounds that their tangent is continuum rather than consistent. That is a defensible
    // engineering choice, and it was an undeclared one: nothing measured whether an HS answer
    // depends on it. KV-CST-002, the oedometer, is the first HS boundary-value problem in the
    // corpus, and its settlement has a closed form to be judged against.
    std::printf("\n== Hardening Soil at its default 1%% residual, and tighter (KV-CST-002) ==\n");
    m::Project hs;
    const std::string hs_path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
    if (!m::load_project(hs_path, hs, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", hs_path.c_str(), err.c_str());
        return 1;
    }
    const auto HM = katai::app::mesh_from_project(hs);
    check(HM.ok, "oedometer mesh built from the case file");
    if (!HM.ok) return 1;

    // The closed form of the loading step, as in the corpus oracle: with c = 0 the HS stiffness
    // factor is (sigma/p_ref)^m, so -eps_1 = (p_ref^m/Eoed_ref)[sigma^(1-m)]/(1-m) over 50..200.
    const double coef = std::pow(100.0, 0.5) / 30000.0;
    const double want = coef * (std::pow(200.0, 0.5) - std::pow(50.0, 0.5)) / 0.5 * 4.0;
    const double hs_tol[3] = {1e-2, 1e-4, 1e-6};
    double u[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        katai::app::NumericalControls nc;
        nc.tolerance = hs_tol[i];
        const auto res = katai::app::solve_phases(
            hs, HM.mesh, katai::app::initial_phase_from(hs.initial_procedure), nullptr, nullptr,
            nc);
        check(res.size() == 2 && res[1].ok,
              "the oedometer converges at tolerance " + std::to_string(hs_tol[i]));
        if (res.size() != 2 || !res[1].ok) return 1;
        int top = 0;
        double bd = 1e300;
        for (int n = 0; n < res[1].mesh.node_count; ++n) {
            const double d = std::hypot(res[1].mesh.x[n] - 0.5, res[1].mesh.y[n] - 4.0);
            if (d < bd) { bd = d; top = n; }
        }
        u[i] = -res[1].disp[top * 2 + 1];
        std::printf("  tolerated residual %.0e -> settlement %.6f m  (%+.3f%% vs closed form)\n",
                    hs_tol[i], u[i], 100.0 * (u[i] - want) / want);
    }
    const double hs_lo = std::fmin(u[0], std::fmin(u[1], u[2]));
    const double hs_hi = std::fmax(u[0], std::fmax(u[1], u[2]));
    const double hs_spread = (hs_hi - hs_lo) / u[0];
    std::printf("  spread from 1e-2 (the default) to 1e-6 = %.4f%% of the default run\n",
                100.0 * hs_spread);
    check(hs_spread < 0.01,
          "the Hardening Soil answer is independent of its stopping rule to within 1%");
    check(std::fabs(u[0] - want) / want < 0.03,
          "and the default 1% residual still lands within 3% of the closed form");
    // 1e-4 is already converged: tightening another two orders changes nothing. So the cost of
    // the default is the 1e-2 -> 1e-4 step and no more, which is a bounded, quotable number.
    check(std::fabs(u[1] - u[2]) / u[1] < 5e-4,
          "1e-4 and 1e-6 agree to 0.05%: the residual is converged by 1e-4");
    std::printf("  the DEFAULT costs %.3f%% on this quantity, and it is not accuracy: the\n"
                "  converged run sits at %+.3f%% from the closed form while the looser default\n"
                "  sits at %+.3f%% -- closer by cancellation, not by being better.\n",
                100.0 * std::fabs(u[0] - u[2]) / u[2], 100.0 * (u[2] - want) / want,
                100.0 * (u[0] - want) / want);

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
