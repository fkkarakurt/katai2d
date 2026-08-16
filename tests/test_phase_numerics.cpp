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

// verify: KV-NUM-009
//   oracle:   closed_form
//   source:   the Hardening Soil oedometric stiffness law as published in the PLAXIS 2D Material Models Manual, E_oed = E_oed^ref ((c cos(phi) + sigma_1 sin(phi))/(c cos(phi) + p_ref sin(phi)))^m, integrated over one-dimensional primary loading -- the same closed form KV-CST-002 is measured against; the error-decomposition procedure is the solution-verification one this program already applies to meshes (Roache 1994; Celik et al. 2008, ASME J. Fluids Eng. 130(7):078001), here applied on the axis a path-dependent model actually discretises
//   locator:  tests/corpus/kv-cst-002-hs-oedometer.k2d, swept on three axes independently -- mesh density 0.5/0.25/0.125 m, load increments 10/20/40/80/160 at a converged tolerance with the seating phase pinned, and the stress range walked over 50-100, 100-200 and 200-400 kPa
//   quantity: settlement of the oedometer top [m] on each sweep, and the observed order of convergence of the load-step sweep [-]
//   expected: the mesh contributes nothing (a weightless column has a uniform strain field, which is exact in the element space); the load path dominates; its observed order is NOT stable across triplets of one sweep, so no GCI band may be quoted from it; and the residual left over is a model deviation, identified by growing with stress level and changing sign near p_ref rather than staying a constant fraction
//   band:     no band is published for this case, deliberately, and the measurements are why. Mesh: 85 -> 1105 nodes changes the answer by 5e-15 relative, which is round-off. Load path: 10 -> 160 increments moves it +3.47% -> +0.54% (a 16x refinement worth 2.9 percentage points against the mesh's zero), and the observed order from the three overlapping triplets of that one sweep is 2.204, 0.469, 1.132 -- a factor of 4.7 apart, so the triplets are not in an asymptotic range and a Richardson extrapolation built on any of them would be an invented number. Tolerance: 1e-6 and 1e-8 agree to six figures, so the sweep above is tolerance-converged, while the Hardening Soil default of 1e-2 is not and makes the same sweep non-monotone. Residual: -0.2425% over 50-100 kPa, +0.6419% over 100-200, +1.0952% over 200-400 -- near zero at the reference pressure and growing away from it with a sign change, the signature of a cap calibrated at p_ref, not of a discretisation. Ceiling: refining the seating phase to 160 increments at 1e-6 does not converge at all, because the tolerated error is an absolute residual and the increment it must satisfy keeps shrinking

// verify: KV-NUM-010
//   oracle:   closed_form
//   source:   Terzaghi (1943) one-dimensional consolidation, degree-of-consolidation series U(Tv) = 1 - sum_j (2/M^2) exp(-M^2 Tv), M = (2j+1) pi/2, Tv = cv t / H_dr^2, cv = k Eoed / gamma_w -- the same closed form KV-CON-002 is measured against; the estimator is the solution-verification procedure of Roache (1994) and Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, applied on the TIME axis, whose order is fixed in advance by the scheme: katai/analysis/consolidation.hpp integrates fully implicitly (alpha = 1), i.e. backward Euler, whose global error is O(dt)
//   locator:  tests/corpus/kv-con-002-terzaghi-column.k2d with the phase ended at Tv = 0.5, swept over 30/60/120/240 time steps at the file's own 0.4 m tri6 mesh, and over 0.8/0.4/0.2 m at a pinned 240 steps
//   quantity: degree of consolidation U at Tv = 0.5 [-]; the observed order of the time refinement [-]; and the GCI band on the file's own 120 steps [%]
//   expected: the observed order is 1, because backward Euler is first order and that is decided before the run; it is STABLE across overlapping triplets, unlike the load path of KV-NUM-009; Richardson extrapolation of three time steps recovers Terzaghi's closed form; and the band contains the true error, which can be checked here rather than trusted because the exact answer is known
//   band:     +/- 0.2480% on the file's own 120 steps (GCI at the observed order, Fs = 1.25), and the actual error there is 0.1965%, so the band contains it. Observed order 0.9850 from 120/60/30 and 0.9924 from 240/120/60 -- the two agree to 0.007, which is what an asymptotic range looks like and is exactly what KV-NUM-009's load path could not produce (2.204 / 0.469 / 1.132). Richardson U(dt->0) = 0.763961996 against the series' 0.763950331, +0.00153%. The mesh, again, is not the axis: 169 -> 1884 nodes moves U by 1.7e-8

// verify: KV-NUM-011
//   oracle:   closed_form
//   source:   1D SH site response of a damped elastic shear column on a rigid base at fundamental-mode resonance (Kramer 1996, Geotechnical Earthquake Engineering, ch. 7): |u_surf| = (4/pi) A/(w_1^2 2 xi), w_1 = 2 pi Vs/(4H) -- the same closed form KV-DYN-002 is measured against; the estimator is Roache (1994) / Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, applied on the TIME axis. The order Newmark alone justifies is 2: katai/analysis/dynamics.hpp integrates with gamma = 1/2, beta = 1/4 (average acceleration, no numerical damping)
//   locator:  tests/corpus/kv-dyn-002-resonant-column.k2d, swept on two axes -- the DURATION over 10/20/40 cycles at a fixed 0.0025 s step, and the time step over 800/1600/3200 steps at 40 cycles where the resonant buildup is finished -- plus two mesh densities
//   quantity: peak surface displacement |u_surf| of the resonant phase [m], and the observed order of the time refinement [-]
//   expected: the duration is an axis of its own, because a resonant amplitude is built up rather than imposed; the time axis is a proper refinement axis and extrapolates onto the closed form; and its observed order is NOT the scheme's 2
//   band:     the prediction FAILED and that is the finding. Observed order 3.147 at 40 cycles, 3.147 and 2.998 at 160 cycles, and 3.097 / 2.795 on a grid deliberately INCOMMENSURATE with the period (T/dt = 18.43, 36.86, 73.72, 147.45), which eliminates sampling phase-lock as the explanation -- the order is about 3, reproducibly, where average-acceleration Newmark alone justifies 2. What the quantity is explains it: a PEAK of a resonant response, not the response at an instant, and a derived quantity need not inherit its scheme's order. So "the algebra fixes the order before the run", true for KV-STR-003's q h^2/12 and for KV-NUM-010's backward Euler, is NOT a general rule, and this is the case that bounds it. The axis is nonetheless clean: Richardson lands on the closed form to -0.036% at 40 cycles and -0.031% at 160, from two independent grids agreeing to 1e-5. Duration: the shipped 20 cycles is 0.21% short of steady state (envelope 1 - exp(-xi w t) = 0.9981) and 40 cycles is converged, so at the file's own settings the buildup is worth about as much as the time step. Mesh: 459 -> 6389 nodes moves |u| by 4.5e-6

#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/math/grid_convergence.hpp>
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

    // --- 6. WHERE THE NUMERICAL ERROR OF THIS CASE ACTUALLY LIVES (KV-NUM-009) ---------------
    // Everything above shows the controls are read. This asks the question underneath it: of the
    // +1% this case sits from its closed form, how much is the mesh, how much the stopping rule,
    // how much the load path -- and how much is not numerical at all. Until now the corpus has
    // only ever banded the MESH, which for a path-dependent model is the wrong axis to start on.
    std::printf("\n== where the numerical error lives (KV-NUM-009) ==\n");

    // (a) THE MESH CONTRIBUTES NOTHING HERE, AND THAT IS MEASURABLE RATHER THAN ARGUABLE. The
    //     column is weightless, so sigma_1 is the surcharge and the strain field is UNIFORM; a
    //     uniform field lies exactly in the element space, so refinement has nothing to improve.
    //     Three densities over a 13x node count say so. This is the estimator's `Exact` branch
    //     meeting a real case, and it is why a mesh band for this case would be a fiction.
    double u_mesh[3] = {0, 0, 0};
    int nodes[3] = {0, 0, 0};
    const double sizes[3] = {0.5, 0.25, 0.125};
    bool mesh_ok = true;
    for (int i = 0; i < 3; ++i) {
        m::Project pr = base;
        pr.mesh.elem_size = sizes[i];
        const auto Mi = katai::app::mesh_from_project(pr);
        if (!Mi.ok) { mesh_ok = false; break; }
        nodes[i] = Mi.mesh.node_count;
        u_mesh[i] = settlement(pr, Mi.mesh, {}, &ok);
        if (!ok) { mesh_ok = false; break; }
        std::printf("  elem %5.3f m  %5d nodes  settlement = %.12e\n", sizes[i], nodes[i],
                    u_mesh[i]);
    }
    check(mesh_ok, "the case solves at all three mesh densities");
    if (mesh_ok) {
        const double d1 = std::fabs(u_mesh[1] - u_mesh[0]) / u_mesh[0];
        const double d2 = std::fabs(u_mesh[2] - u_mesh[0]) / u_mesh[0];
        std::printf("  relative change over a %.0fx node count: %.2e and %.2e\n",
                    (double)nodes[2] / nodes[0], d1, d2);
        check(nodes[2] > 10 * nodes[0], "the refinement is a real one (>10x the nodes)");
        check(d1 < 1e-8 && d2 < 1e-8,
              "the mesh contributes nothing: a uniform field is already exact in the element space");
    }

    // (b) THE LOAD PATH IS THE DOMINANT DISCRETISATION -- AND IT IS NOT A RICHARDSON PARAMETER.
    //     Refining the increments moves this answer far more than any mesh does. But the
    //     estimator needs phi(d) = phi_exact + C d^p, and a path-dependent integration with its
    //     own adaptive substepping does not supply one: the observed order computed from
    //     successive triplets of the SAME sweep disagrees with itself by a factor of several.
    //     Quoting a GCI from it would dress that disagreement up as a confidence interval, so
    //     this case reports the measured SPREAD instead and says why.
    //     The seating phase is pinned throughout: refining IT is a separate effect, see (d).
    const int path_steps[5] = {10, 20, 40, 80, 160};
    double u_path[5] = {0, 0, 0, 0, 0};
    bool path_ok = true;
    for (int i = 0; i < 5; ++i) {
        m::Project pr = base;
        pr.initial.load_steps = 40;            // PINNED
        pr.initial.tolerance = 1e-6;
        pr.initial.max_iterations = 500;
        pr.phases[0].load_steps = path_steps[i];
        pr.phases[0].tolerance = 1e-6;         // converged: 1e-8 agrees to six figures
        pr.phases[0].max_iterations = 500;
        u_path[i] = settlement(pr, M.mesh, {}, &ok);
        if (!ok) { path_ok = false; break; }
        std::printf("  %4d load increments  settlement = %.12e\n", path_steps[i], u_path[i]);
    }
    check(path_ok, "the case solves at every step count in the sweep");
    if (path_ok) {
        const double spread = std::fabs(u_path[0] - u_path[4]) / u_path[4];
        std::printf("  spread over a 16x refinement: %.3f%% (the mesh gave %.1e)\n",
                    100.0 * spread, mesh_ok ? std::fabs(u_mesh[2] - u_mesh[0]) / u_mesh[0] : 0.0);
        check(spread > 5e-3, "the load path moves the answer by more than half a percent");

        // Three overlapping triplets of one monotone-looking sweep, three different orders.
        double p_obs[3] = {0, 0, 0};
        bool orders_ok = true;
        for (int i = 0; i < 3; ++i) {
            const double f = u_path[i + 2], mid = u_path[i + 1], c = u_path[i];
            const double e21 = mid - f, e32 = c - mid;
            if (e21 == 0.0 || e32 / e21 <= 0.0) { orders_ok = false; break; }
            p_obs[i] = std::log(std::fabs(e32 / e21)) / std::log(2.0);
        }
        check(orders_ok, "each triplet is at least monotone, so an order can be computed at all");
        if (orders_ok) {
            std::printf("  observed order from triplets %d/%d/%d, %d/%d/%d, %d/%d/%d: "
                        "%.3f, %.3f, %.3f\n",
                        path_steps[2], path_steps[1], path_steps[0], path_steps[3], path_steps[2],
                        path_steps[1], path_steps[4], path_steps[3], path_steps[2], p_obs[0],
                        p_obs[1], p_obs[2]);
            double lo = p_obs[0], hi = p_obs[0];
            for (int i = 1; i < 3; ++i) { lo = std::fmin(lo, p_obs[i]); hi = std::fmax(hi, p_obs[i]); }
            std::printf("  they span a factor of %.1f -- NOT an asymptotic range\n", hi / lo);
            check(hi / lo > 2.0,
                  "the observed order is not stable, so no GCI band may be quoted from the path");
        }
    }

    // (c) WHAT SURVIVES ALL OF IT IS NOT NUMERICAL. The residual left after an exact mesh, a
    //     converged tolerance and a 16x refined path is still ~0.5% -- and its SIGNATURE says
    //     what it is. A pure offset in the fitted stiffness would show the same relative
    //     deviation on every stress range. This one is near zero close to p_ref = 100 kPa and
    //     grows away from it, changing sign: the mark of a cap whose alpha/beta were calibrated
    //     at the reference pressure (hs_calibrate_cap, as PLAXIS derives them by simulating an
    //     oedometer). It is a MODEL deviation from the closed form, not a mesh or step artefact,
    //     and reporting it as numerical uncertainty would be wrong in both directions.
    const double ranges[3][2] = {{50.0, 100.0}, {100.0, 200.0}, {200.0, 400.0}};
    double dev[3] = {0, 0, 0};
    bool range_ok = true;
    for (int i = 0; i < 3; ++i) {
        m::Project pr = base;
        pr.loads[0].qy1 = pr.loads[0].qy2 = -ranges[i][0];
        pr.loads[1].qy1 = pr.loads[1].qy2 = -(ranges[i][1] - ranges[i][0]);
        pr.initial.load_steps = 40; pr.initial.tolerance = 1e-6; pr.initial.max_iterations = 500;
        pr.phases[0].load_steps = 160; pr.phases[0].tolerance = 1e-6;
        pr.phases[0].max_iterations = 500;
        const double got = settlement(pr, M.mesh, {}, &ok);
        if (!ok) { range_ok = false; break; }
        // The closed form of the HS oedometric law over this range (c = 0, so sin(phi) cancels).
        const double kH = 4.0, kEoed = 30000.0, kPref = 100.0, kM = 0.5;
        const double want = std::pow(kPref, kM) / kEoed *
                            (std::pow(ranges[i][1], 1.0 - kM) - std::pow(ranges[i][0], 1.0 - kM)) /
                            (1.0 - kM) * kH;
        dev[i] = (got - want) / want;
        std::printf("  %3.0f -> %3.0f kPa  settlement %.9f vs closed form %.9f  %+.4f%%\n",
                    ranges[i][0], ranges[i][1], got, want, 100.0 * dev[i]);
    }
    check(range_ok, "the case solves over each stress range");
    if (range_ok) {
        check(dev[0] < dev[1] && dev[1] < dev[2],
              "the residual grows with the stress level: it is not a constant offset");
        check(dev[0] < 0.0 && dev[2] > 0.0,
              "and it changes sign near p_ref, which is where the cap was calibrated");
    }

    // (d) A DECLARED CEILING ON PATH REFINEMENT, because a study that cannot be repeated is not
    //     a study. The tolerated error is an ABSOLUTE force residual, so shrinking the increment
    //     does not shrink what each increment must achieve. On this weightless column -- whose
    //     confining stress starts near zero, where the Hardening Soil stiffness is at its
    //     smallest -- refining the SEATING phase to 160 increments at a converged tolerance stops
    //     converging altogether rather than getting better. That is why (b) pins it at 40. If
    //     this check ever fails because the run now succeeds, the ceiling has moved and the
    //     sentence above needs rewriting, which is the point of pinning it.
    m::Project seat_fine = base;
    seat_fine.initial.load_steps = 160;
    seat_fine.initial.tolerance = 1e-6;
    seat_fine.initial.max_iterations = 500;
    settlement(seat_fine, M.mesh, {}, &ok);
    std::printf("  seating phase at 160 increments, tol 1e-6: %s\n", ok ? "converged" : "does NOT converge");
    check(!ok, "refining the seating phase past the ceiling fails openly instead of drifting");

    // --- 7. THE SAME QUESTION WHERE THE AXIS *IS* CLEAN (KV-NUM-010) -------------------------
    // KV-NUM-009 refused a band. A refusal is only worth something if the same procedure, applied
    // to a case whose axis IS a proper discretisation parameter, produces one -- otherwise the
    // refusal is indistinguishable from the machinery not working. Terzaghi consolidation is that
    // case, and its order is known before any run: consolidation.hpp integrates FULLY IMPLICITLY
    // (alpha = 1), which is backward Euler, whose global error is O(dt). The prediction is p = 1,
    // decided by the algebra, exactly as KV-STR-003's peak moment was decided at p = 2.
    std::printf("\n== the same procedure where the axis is clean (KV-NUM-010) ==\n");

    m::Project tz;
    const std::string tz_path = std::string(KATAI_CORPUS_DIR) + "/kv-con-002-terzaghi-column.k2d";
    if (!m::load_project(tz_path, tz, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", tz_path.c_str(), err.c_str());
        return 1;
    }
    check(!tz.phases.empty(), "the consolidation case has a phase to control");

    // The end of the phase is moved to Tv = 0.5, where U is still moving fast enough to measure a
    // step error against. Every run then ends at the SAME physical time, so what is compared is
    // one quantity rather than one label.
    const double tzH = 12.0, tzE = 1000.0, tzK = 0.1, tzQ = 10.0;
    const double cv = tzK * tzE / katai::app::kGammaWater;   // nu = 0 -> Eoed = E
    const double s_inf = tzQ * tzH / tzE;
    const double Tv_end = 0.5;
    const double t_end = Tv_end * tzH * tzH / cv;
    double u_exact = 1.0;                                    // Terzaghi's series at Tv_end
    for (int j = 0; j < 80; ++j) {
        const double Mj = (2 * j + 1) * 3.14159265358979323846 / 2.0;
        u_exact -= (2.0 / (Mj * Mj)) * std::exp(-Mj * Mj * Tv_end);
    }
    std::printf("  Tv = %.2f, Terzaghi U = %.9f (closed form, series)\n", Tv_end, u_exact);

    const auto degree_of_consolidation = [&](int nsteps, double elem, int* nodes) -> double {
        m::Project pr = tz;
        pr.phases[0].duration = t_end;
        pr.phases[0].time_steps = nsteps;
        if (elem > 0.0) pr.mesh.elem_size = elem;
        const auto Mt = katai::app::mesh_from_project(pr);
        if (!Mt.ok) return -1.0;
        if (nodes) *nodes = Mt.mesh.node_count;
        const auto r = katai::app::solve_phases(
            pr, Mt.mesh, katai::app::initial_phase_from(pr.initial_procedure));
        if (r.size() != 2 || !r[1].ok || r[1].consol_settlement.empty()) return -1.0;
        return r[1].consol_settlement.back() / s_inf;
    };

    // (a) The time axis, four densities so the order can be checked for STABILITY and not just
    //     computed once -- which is precisely what KV-NUM-009's load path failed.
    const int nst[4] = {30, 60, 120, 240};
    double U[4] = {0, 0, 0, 0};
    bool time_ok = true;
    for (int i = 0; i < 4; ++i) {
        U[i] = degree_of_consolidation(nst[i], 0.0, nullptr);
        if (U[i] < 0.0) { time_ok = false; break; }
        std::printf("  %4d time steps  U = %.9f  (%+.4f%% vs Terzaghi)\n", nst[i], U[i],
                    100.0 * (U[i] - u_exact) / u_exact);
    }
    check(time_ok, "the consolidation case solves at every time-step count");

    if (time_ok) {
        // The window comes from the INTEGRATOR, not from the elements. The default OrderPolicy
        // brackets [0.5, 4.0] because that is what a tri6 mesh can deliver; a backward-Euler time
        // axis is first order, so an observed order far from 1 here is evidence that the triplet
        // is not asymptotic rather than evidence about the discretisation, and the assumed order
        // to fall back on is 1 rather than 2.
        katai::math::OrderPolicy time_axis;
        time_axis.p_assumed = 1.0;
        time_axis.p_min = 0.75;
        time_axis.p_max = 1.5;

        // The file's own 120 steps as the finest of the reported triplet, so the band published
        // is the band for the run the corpus actually ships.
        katai::math::GridTriplet g;
        g.h1 = t_end / 120.0; g.h2 = t_end / 60.0; g.h3 = t_end / 30.0;
        g.phi1 = U[2];        g.phi2 = U[1];       g.phi3 = U[0];
        const auto e = katai::math::grid_convergence_band(g, time_axis);
        std::printf("  triplet 120/60/30: %s, observed order p = %.4f (backward Euler predicts 1)\n",
                    katai::math::convergence_kind_name(e.kind), e.p);
        std::printf("  Richardson U(dt->0) = %.9f  (%+.5f%% vs Terzaghi)\n", e.phi_extrapolated,
                    100.0 * (e.phi_extrapolated - u_exact) / u_exact);
        std::printf("  band on the file's own 120 steps: +/- %.4f%% (%s)\n", 100.0 * e.band,
                    e.band_basis.c_str());
        check(e.ok && e.kind == katai::math::ConvergenceKind::MonotonicConvergence,
              "the time refinement converges monotonically");
        check(std::fabs(e.p - 1.0) < 0.1,
              "the observed order is the one backward Euler was known to have: 1");
        check(e.asymptotic, "the triplet is inside the asymptotic range, so the order may be quoted");
        check(std::fabs(e.phi_extrapolated - u_exact) / u_exact < 5e-4,
              "Richardson recovers Terzaghi's closed form from three time steps alone");

        // The order is STABLE, which is the property KV-NUM-009's load path lacked. Same sweep,
        // the other overlapping triplet.
        katai::math::GridTriplet g2;
        g2.h1 = t_end / 240.0; g2.h2 = t_end / 120.0; g2.h3 = t_end / 60.0;
        g2.phi1 = U[3];        g2.phi2 = U[2];        g2.phi3 = U[1];
        const auto e2 = katai::math::grid_convergence_band(g2, time_axis);
        std::printf("  triplet 240/120/60: p = %.4f -- the two triplets agree to %.3f\n", e2.p,
                    std::fabs(e2.p - e.p));
        check(std::fabs(e2.p - e.p) < 0.05,
              "and it is STABLE across triplets: this axis really is asymptotic");

        // The band has to CONTAIN the error it claims to bound, and here that can be checked
        // rather than trusted, because the exact answer is known.
        const double actual = std::fabs(U[2] - u_exact) / u_exact;
        std::printf("  actual error at 120 steps %.4f%% vs band %.4f%%\n", 100.0 * actual,
                    100.0 * e.band);
        check(actual < e.band, "the published band contains the true error");
    }

    // (b) And the mesh, again, is not where the error is. Three densities at a pinned 240 steps.
    //     Two corpus cases in a row whose error lives on an axis no mesh sweep would have found.
    double Um[3] = {0, 0, 0};
    int nm[3] = {0, 0, 0};
    const double tz_sizes[3] = {0.8, 0.4, 0.2};
    bool tz_mesh_ok = true;
    for (int i = 0; i < 3; ++i) {
        Um[i] = degree_of_consolidation(240, tz_sizes[i], &nm[i]);
        if (Um[i] < 0.0) { tz_mesh_ok = false; break; }
        std::printf("  elem %.1f m  %5d nodes  U = %.9f\n", tz_sizes[i], nm[i], Um[i]);
    }
    check(tz_mesh_ok, "the consolidation case solves at all three mesh densities");
    if (tz_mesh_ok) {
        const double dm = std::fabs(Um[2] - Um[0]) / Um[0];
        std::printf("  relative change over a %.0fx node count: %.2e\n", (double)nm[2] / nm[0], dm);
        check(dm < 1e-6,
              "the mesh contributes essentially nothing to an integral quantity like U either");
    }

    // --- 8. WHERE "THE ALGEBRA DECIDES THE ORDER" STOPS BEING TRUE (KV-NUM-011) --------------
    // KV-STR-003 knew its order was 2 before any run (the moment overshoot is q h^2/12) and
    // KV-NUM-010 knew its order was 1 (backward Euler). It is tempting to generalise that, and
    // this case is the counterexample that stops it. Newmark here is gamma = 1/2, beta = 1/4 --
    // average acceleration, no numerical damping, SECOND order -- so the prediction was p = 2.
    // The measurement says 3, reproducibly. The order belongs to the QUANTITY as much as to the
    // scheme, and the published quantity here is the PEAK of a resonant response, not the
    // response at an instant.
    std::printf("\n== when the scheme's order is not the quantity's order (KV-NUM-011) ==\n");

    m::Project dy;
    const std::string dy_path = std::string(KATAI_CORPUS_DIR) + "/kv-dyn-002-resonant-column.k2d";
    if (!m::load_project(dy_path, dy, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", dy_path.c_str(), err.c_str());
        return 1;
    }
    check(!dy.phases.empty(), "the dynamic case has a phase to control");

    const double dyE = 208000.0, dyNu = 0.3, dyGam = 19.62, dyH = 20.0, dyA = 1.0, dyXi = 0.05;
    const double Gs = dyE / (2.0 * (1.0 + dyNu)), rho = dyGam / katai::app::kGammaWater * 1.0;
    const double f1 = std::sqrt(Gs / (dyGam / 9.81)) / (4.0 * dyH);
    const double w1 = 2.0 * 3.14159265358979323846 * f1;
    const double u_ex = (4.0 / 3.14159265358979323846) * dyA / (w1 * w1 * 2.0 * dyXi);
    (void)rho;
    std::printf("  f_1 = %.4f Hz, steady-state |u_surf| = %.9f m (closed form)\n", f1, u_ex);

    const auto peak_u = [&](int steps, double duration, double elem) -> double {
        m::Project pr = dy;
        pr.phases[0].duration = duration;
        pr.phases[0].time_steps = steps;
        if (elem > 0.0) pr.mesh.elem_size = elem;
        const auto Md = katai::app::mesh_from_project(pr);
        if (!Md.ok) return -1.0;
        const auto r = katai::app::solve_phases(
            pr, Md.mesh, katai::app::initial_phase_from(pr.initial_procedure));
        if (r.size() < 2 || !r[1].ok) return -1.0;
        return r[1].max_disp;
    };

    // (a) THE DURATION IS AN AXIS OF ITS OWN, and at the shipped settings it costs about as much
    //     as the time step. A resonant amplitude is BUILT UP, not imposed: the envelope goes like
    //     1 - exp(-xi w t), so the file's 20 cycles is still 1 - exp(-2 pi 0.05 20) = 0.9981 of
    //     the way there. Sweeping dt without knowing that would attribute the shortfall to the
    //     integrator. dt is held at 0.0025 s throughout so only the duration moves.
    const double cycles[3] = {10.0, 20.0, 40.0};
    double u_dur[3] = {0, 0, 0};
    bool dur_ok = true;
    for (int i = 0; i < 3; ++i) {
        const double dur = cycles[i] / f1;
        u_dur[i] = peak_u((int)std::llround(dur / 0.0025), dur, 0.0);
        if (u_dur[i] < 0.0) { dur_ok = false; break; }
        std::printf("  %4.0f cycles (%5.1f s)  |u| = %.9f  (%+.4f%%)  envelope %.6f\n", cycles[i],
                    dur, u_dur[i], 100.0 * (u_dur[i] - u_ex) / u_ex,
                    1.0 - std::exp(-2.0 * 3.14159265358979323846 * dyXi * cycles[i]));
    }
    check(dur_ok, "the dynamic case solves at every duration");
    if (dur_ok) {
        check(std::fabs(u_dur[1] - u_ex) / u_ex > 1.5e-3,
              "at the shipped 20 cycles the resonant amplitude is still measurably short");
        check(std::fabs(u_dur[2] - u_ex) / u_ex < 5e-4,
              "by 40 cycles the buildup is done, so a time sweep there measures the time step");
    }

    // (b) THE TIME AXIS, at a duration where the buildup is finished so it is the only thing
    //     moving. It IS a proper refinement axis -- monotone, extrapolable, and the extrapolated
    //     value lands on the closed form. What it is NOT is second order.
    const int dsteps[3] = {800, 1600, 3200};
    const double dyn_dur = 40.0 / f1;   // 16 s
    double u_dt[3] = {0, 0, 0};
    bool dt_ok = true;
    for (int i = 0; i < 3; ++i) {
        u_dt[i] = peak_u(dsteps[i], dyn_dur, 0.0);
        if (u_dt[i] < 0.0) { dt_ok = false; break; }
        std::printf("  %5d steps (dt = %.6f s, %.0f per cycle)  |u| = %.9f  (%+.4f%%)\n",
                    dsteps[i], dyn_dur / dsteps[i], (1.0 / f1) / (dyn_dur / dsteps[i]), u_dt[i],
                    100.0 * (u_dt[i] - u_ex) / u_ex);
    }
    check(dt_ok, "the dynamic case solves at every time-step count");
    if (dt_ok) {
        katai::math::OrderPolicy dyn_axis;
        dyn_axis.p_assumed = 2.0;   // what Newmark alone would justify
        dyn_axis.p_min = 0.75;
        dyn_axis.p_max = 4.0;
        katai::math::GridTriplet g;
        g.h1 = dyn_dur / dsteps[2]; g.h2 = dyn_dur / dsteps[1]; g.h3 = dyn_dur / dsteps[0];
        g.phi1 = u_dt[2];           g.phi2 = u_dt[1];           g.phi3 = u_dt[0];
        const auto e = katai::math::grid_convergence_band(g, dyn_axis);
        std::printf("  %s, observed order p = %.4f -- Newmark alone would give 2\n",
                    katai::math::convergence_kind_name(e.kind), e.p);
        std::printf("  Richardson |u|(dt->0) = %.9f  (%+.5f%% vs the closed form)\n",
                    e.phi_extrapolated, 100.0 * (e.phi_extrapolated - u_ex) / u_ex);
        std::printf("  band on this triplet: +/- %.4f%% (%s)\n", 100.0 * e.band,
                    e.band_basis.c_str());
        check(e.ok && e.kind == katai::math::ConvergenceKind::MonotonicConvergence,
              "the time refinement converges monotonically: it is a proper axis");
        check(std::fabs(e.phi_extrapolated - u_ex) / u_ex < 1e-3,
              "and it extrapolates onto the closed-form resonant amplitude");
        // THE POINT. Asserted as an inequality against the predicted order, because the finding
        // is that the prediction is wrong, and a test that merely recorded 3.1 would go quiet the
        // day the mechanism changed. Measured 3.147 here; 3.147 and 2.998 at 160 cycles; and
        // 3.097 / 2.795 on a grid deliberately INCOMMENSURATE with the period, which eliminates
        // sampling phase-lock as the explanation. What remains is that this quantity is a peak of
        // a resonant response, and a peak does not have to inherit the scheme's order.
        check(e.p > 2.3,
              "the observed order is NOT the scheme's 2: the quantity has an order of its own");
    }

    // (c) And the mesh, for the third corpus case in a row, is not the axis.
    const double u_m1 = peak_u(3200, dyn_dur, 0.8);
    const double u_m2 = peak_u(3200, dyn_dur, 0.4);
    check(u_m1 > 0.0 && u_m2 > 0.0, "the dynamic case solves at both mesh densities");
    if (u_m1 > 0.0 && u_m2 > 0.0) {
        std::printf("  elem 0.8 m |u| = %.9f   elem 0.4 m |u| = %.9f   relative change %.2e\n",
                    u_m1, u_m2, std::fabs(u_m2 - u_m1) / u_m1);
        check(std::fabs(u_m2 - u_m1) / u_m1 < 1e-4, "the mesh is not where this error lives either");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
