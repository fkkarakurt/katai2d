// Numerics belong in the file, because a published number is a claim about two things.
//
// "The settlement is 19.0 mm" is a statement about a model AND about the tolerated error, the
// number of load increments and the iteration limit it was solved with. Until now KATAI could be
// TOLD those numbers -- there is a jobs-layer seam, and KV-NUM-007 uses it -- but it could not be
// ASKED for them in a .k2d, so a project handed to a reviewer carried the model and dropped the
// numerics. The reviewer then re-runs it at whatever this build happens to choose, gets a
// different number, and has no way to tell whether the model changed or the stopping rule did.
// PLAXIS puts these on the phase (Reference Manual, "Numerical control parameters"), and so does
// this, for the reason the manual gives them: they are part of the calculation, not part of the
// program.
//
// What this case checks is deliberately narrow and deliberately hostile: that the value in the
// file REACHES the solver. A control that is read, validated, echoed in the GUI and then quietly
// not used produces a run indistinguishable from one that honoured it -- which is exactly how
// the strength-reduction search came to have a hard-coded tolerance while a whole verification
// case believed it was sweeping one (see the note at the top of KV-NUM-007). So every control
// here is asserted twice: the file must give the SAME answer as the seam (it is the same
// control) and a DIFFERENT answer from the default (it is being read at all).
//
// verify: KV-NUM-008
//   oracle:   independent_path
//   source:   the same numerical control applied through two independent routes -- the .k2d file (model::Phase::tolerance / load_steps / max_iterations, io v7) and the jobs-layer seam (app::NumericalControls, the route KV-NUM-007 measures with) -- must reach the same solver and produce the same run; the controls themselves are the phase-level numerical control parameters of PLAXIS 2D 2025.1 Reference Manual (Tolerated error, Max iterations), whose per-phase presence in the input is what makes a calculation reproducible by a third party
//   locator:  tests/corpus/kv-cst-002-hs-oedometer.k2d (Hardening Soil, whose answer is known from KV-NUM-007 to move with the tolerance) solved four ways per control: default, control set in the FILE, the same control passed through the seam, and both set at once with different values
//   quantity: settlement of the oedometer top [m] and the file round trip of the three control fields
//   expected: file == seam bit-for-bit; file != default (the control is read); seam wins when both are set (documented precedence); and the three fields survive a write/read round trip
//   band:     exact -- these are identity checks, not approximations. Measured on this tree: default 0.018934049 m; tolerance 1e-6 from the file 0.019044410 m, identical to the seam to 0.0e+00 and 0.58% away from the default, so the control is unmistakably read; on the staged phase alone 0.019026773 m, which differs again and is what "per phase" means; 4 load increments 0.019461550 m; an iteration limit of 8 0.018934204 m (the solver cuts back and arrives by another path)

#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
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

// The settlement of the oedometer top -- the same quantity KV-CST-002 and KV-NUM-007 read.
double settlement(const m::Project& pr, const katai::mesh::Mesh& mesh,
                  const katai::app::NumericalControls& nc, bool* ok_out) {
    const auto res = katai::app::solve_phases(
        pr, mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, nullptr, nc);
    *ok_out = res.size() == 2 && res[1].ok;
    if (!*ok_out) return 0.0;
    int top = 0;
    double best = 1e300;
    for (int n = 0; n < res[1].mesh.node_count; ++n) {
        const double d = std::hypot(res[1].mesh.x[n] - 0.5, res[1].mesh.y[n] - 4.0);
        if (d < best) { best = d; top = n; }
    }
    return -res[1].disp[top * 2 + 1];
}

// The file's controls are PER PHASE; the seam's apply to the whole run. To compare the two
// routes as the same control, the file has to say on every phase what the seam says once --
// including the initial phase, which is a phase like any other and where a run usually spends
// its first and most delicate solve.
m::Project all_phases(const m::Project& base, double tol, int steps, int iters) {
    m::Project pr = base;
    const auto set = [&](m::Phase& ph) {
        if (tol > 0.0) ph.tolerance = tol;
        if (steps > 0) ph.load_steps = steps;
        if (iters > 0) ph.max_iterations = iters;
    };
    set(pr.initial);
    for (auto& ph : pr.phases) set(ph);
    return pr;
}

}  // namespace

int main() {
    std::printf("== numerical controls carried by the file (KV-NUM-008) ==\n");

    m::Project base;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
    if (!m::load_project(path, base, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        return 1;
    }
    check(!base.phases.empty(), "the case has a staged phase to control");
    if (base.phases.empty()) return 1;

    const auto M = katai::app::mesh_from_project(base);
    check(M.ok, "oedometer mesh built from the case file");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }

    bool ok = false;
    const double u_default = settlement(base, M.mesh, {}, &ok);
    check(ok, "the untouched case solves");
    if (!ok) return 1;
    std::printf("  default (the material class chooses)      settlement = %.9f m\n", u_default);

    // --- 1. The tolerated error ------------------------------------------------------------
    // KV-NUM-007 measured what this control is worth on this very problem: 0.59% between the
    // Hardening Soil default of 1e-2 and a converged 1e-6. That is the margin this check needs --
    // large enough that a dropped control cannot hide inside it.
    const m::Project from_file = all_phases(base, 1e-6, 0, 0);
    const double u_file = settlement(from_file, M.mesh, {}, &ok);
    check(ok, "the case with a tolerated error in the FILE solves");
    katai::app::NumericalControls seam;
    seam.tolerance = 1e-6;
    const double u_seam = settlement(base, M.mesh, seam, &ok);
    check(ok, "the same tolerance through the seam solves");
    std::printf("  tolerance 1e-6 from the file              settlement = %.9f m\n", u_file);
    std::printf("  tolerance 1e-6 through the seam           settlement = %.9f m\n", u_seam);
    std::printf("  difference between the two routes         %.3e m\n", std::fabs(u_file - u_seam));
    check(u_file == u_seam, "the file and the seam are the same control, bit for bit");
    std::printf("  difference from the default run           %.3e m (%.2f%%)\n",
                std::fabs(u_file - u_default), 100.0 * std::fabs(u_file - u_default) / u_default);
    check(std::fabs(u_file - u_default) / u_default > 1e-3,
          "and the file's tolerance changes the answer, so it is genuinely being read");

    // PER PHASE means per phase. Tightening only the staged step, and leaving the initial phase
    // to the material class, is a different calculation from tightening both -- and being able
    // to say which is the reason the controls live on the phase rather than on the project.
    m::Project one_phase = base;
    one_phase.phases[0].tolerance = 1e-6;
    const double u_one = settlement(one_phase, M.mesh, {}, &ok);
    check(ok, "the case with the tolerance on ONE phase solves");
    std::printf("  tolerance 1e-6 on the staged phase only   settlement = %.9f m\n", u_one);
    check(u_one != u_file,
          "and it differs from tightening every phase: the control really is per phase");

    // --- 2. The load increments -------------------------------------------------------------
    // A Hardening Soil phase is given 40 increments by default. Four is a coarser path to the
    // same load, and a path-dependent model does not arrive at quite the same place -- which is
    // the point: if the number were ignored, it would.
    const m::Project steps_file = all_phases(base, 0.0, 4, 0);
    const double u_steps_file = settlement(steps_file, M.mesh, {}, &ok);
    check(ok, "the case with 4 load increments in the FILE solves");
    katai::app::NumericalControls steps_seam;
    steps_seam.steps = 4;
    const double u_steps_seam = settlement(base, M.mesh, steps_seam, &ok);
    check(ok, "the same 4 increments through the seam solve");
    std::printf("  4 load increments: file %.9f m, seam %.9f m (default has 40: %.9f m)\n",
                u_steps_file, u_steps_seam, u_default);
    check(u_steps_file == u_steps_seam, "load increments: file == seam, bit for bit");
    check(u_steps_file != u_default, "and 4 increments is not 40, so the count is read");

    // --- 3. The iteration limit -------------------------------------------------------------
    // Eight iterations per increment is less than this model wants, so the solver cuts the
    // increment back and works its way through anyway. The answer that comes out is therefore
    // reached by a different path -- an honest way to see the limit take effect without asking
    // for a run that fails. (Two would not solve at all: the initial phase needs more than that,
    // which is itself worth knowing before setting the limit in a real project.)
    const m::Project iter_file = all_phases(base, 0.0, 0, 8);
    const double u_iter_file = settlement(iter_file, M.mesh, {}, &ok);
    check(ok, "the case with an iteration limit of 8 in the FILE solves");
    katai::app::NumericalControls iter_seam;
    iter_seam.max_iterations = 8;
    const double u_iter_seam = settlement(base, M.mesh, iter_seam, &ok);
    check(ok, "the same limit through the seam solves");
    std::printf("  max 8 iterations:  file %.9f m, seam %.9f m\n", u_iter_file, u_iter_seam);
    check(u_iter_file == u_iter_seam, "iteration limit: file == seam, bit for bit");
    check(u_iter_file != u_default, "and the limit changes the path, so it is read");

    // --- 4. Precedence, stated and tested ---------------------------------------------------
    // When both are set the SEAM wins. The seam exists so that a given file can be re-run at
    // other numerics (that is what KV-NUM-007 does to a checked-in corpus case); a file that
    // always won would make that study impossible to perform.
    const m::Project both = all_phases(base, 1e-2, 0, 0);   // the file asks for the loose default
    katai::app::NumericalControls tight;
    tight.tolerance = 1e-6;               // the caller asks for a converged run
    const double u_both = settlement(both, M.mesh, tight, &ok);
    check(ok, "the case with a control in BOTH places solves");
    std::printf("  file 1e-2 + seam 1e-6 -> %.9f m (seam run was %.9f m)\n", u_both, u_seam);
    check(u_both == u_seam,
          "an explicit override beats the file, so a checked-in case can still be re-measured");

    // --- 5. The file round trip -------------------------------------------------------------
    // A control the writer drops is a control the reviewer never receives.
    m::Project trip;
    m::Project written = base;
    written.phases[0].tolerance = 3.5e-4;
    written.phases[0].load_steps = 17;
    written.phases[0].max_iterations = 33;
    written.initial.tolerance = 1e-7;     // the initial phase carries them too
    const std::string text = m::project_to_json(written);
    check(m::project_from_json(text, trip, &err, nullptr), "the written project reads back");
    check(trip.phases[0].tolerance == 3.5e-4 && trip.phases[0].load_steps == 17 &&
              trip.phases[0].max_iterations == 33 && trip.initial.tolerance == 1e-7,
          "all three controls survive the round trip, on the initial phase as well");
    check(m::project_to_json(base).find("\"tol\"") == std::string::npos,
          "a project that sets no controls writes no control keys (the corpus stays as it is)");
    const katai::io::ValidationReport rep = katai::io::validate_project(written);
    check(rep.ok(), "a project with controls set validates");

    // A tolerated error of 1 accepts a residual the size of the load itself. The validator must
    // refuse it rather than let a run report equilibrium it never reached.
    m::Project absurd = base;
    absurd.phases[0].tolerance = 1.0;
    check(!katai::io::validate_project(absurd).ok(),
          "a tolerated error of 1.0 (100%) is refused, not silently accepted");

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
