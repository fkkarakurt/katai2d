// Numerics belong in the file, because a published number is a claim about two things.
//
// "The settlement is 19.0 mm" is a statement about a model AND about the tolerated error, the
// number of load increments and the iteration limit it was solved with. Until now KATAI could be
// TOLD those numbers -- there is a jobs-layer seam, and KV-NUM-007 uses it -- but it could not be
// ASKED for them in a .k2d, so a project handed to a reviewer carried the model and dropped the
// numerics. The reviewer then re-runs it at whatever this build happens to choose, gets a
// different number, and has no way to tell whether the model changed or the stopping rule did.
// KATAI's input contract puts these on the phase (docs/k2d-format.md, `tol` / `loadsteps` /
// `maxiter` / `substol`), for a plain reason: they are part of the calculation, not part of the
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
//   source:   the same numerical control applied through two independent routes -- the .k2d file (model::Phase::tolerance / load_steps / max_iterations since io v7, and substep_tolerance since io v15) and the jobs-layer seam (app::NumericalControls, the route KV-NUM-007 measures with) -- must reach the same solver and produce the same run; the controls themselves are the phase-level numerical control parameters of KATAI 2D input contract (docs/k2d-format.md, tol / loadsteps / maxiter / substol), whose per-phase presence in the input is what makes a calculation reproducible by a third party
//   locator:  tests/corpus/kv-cst-002-hs-oedometer.k2d (Hardening Soil, whose answer is known from KV-NUM-007 to move with the tolerance) solved four ways per control: default, control set in the FILE, the same control passed through the seam, and both set at once with different values
//   quantity: settlement of the oedometer top [m] and the file round trip of the four control fields
//   expected: file == seam bit-for-bit on every control; the control demonstrably reaches the solver; seam wins when both are set (documented precedence); and the four fields survive a write/read round trip
//   band:     exact -- these are identity checks, not approximations. Measured on this tree (2026-08-24, with the local convergence criteria binding): default 0.018647102 m; tolerance 1e-6 from the file 0.018646941 m, identical to the seam to 0.0e+00; on the staged phase alone 0.018643754 m, which differs again and is what "per phase" means; 4 load increments 0.018643176 m against 40's 0.018647102 m. TWO OF THIS CASE'S GUARDS DECAYED AND WERE REWRITTEN, both for the same reason and both recorded rather than quietly repaired. (1) The tolerated error used to be proved READ by showing the answer moved when it changed; it moves by 0.0009% now, because the default run already stands on the converged answer, so the proof is taken where the control lands instead -- the run reports the tolerance it ran under (1e-6 against the default's 1e-2), it demonstrably MET it, and reaching it cost 473 iterations against 194. (2) The ITERATION LIMIT was pinned at a threshold of 3; it is 5 now. The threshold is RECORDED AND CHECKED rather than pinned-and-asserted: the recorded value has to straddle the boundary (it converges, one below it refuses), which is two solves and exactly the proof a scan would give, and the scan runs only when that straddle stops holding -- and then it reports the value it found. Searching from scratch every run cost up to a hundred two-phase solves and made this file the slowest test in the suite by a factor of two, which is a real price for a property that changes about once a year. There turned out to be TWO thresholds and the old check conflated them: at 5 the run stops REFUSING but survives by cutting increments back, so it walks a different load path and lands 0.1796% away; only from 6 does it reproduce the default bit for bit. A guard that proves a control is read by pointing at a difference stops proving anything when the difference is the defect being fixed -- which has now happened three times on this case

// verify: KV-NUM-009
//   oracle:   closed_form
//   source:   Schanz, T., Vermeer, P.A. & Bonnier, P.G. (1999). The Hardening Soil model: formulation and verification. Beyond 2000 in Computational Geotechnics, Balkema, 281-296 -- the Hardening Soil oedometric stiffness law, integrated over one-dimensional primary loading -- the same closed form KV-CST-002 is measured against; the error-decomposition procedure is the solution-verification one this program already applies to meshes (Roache 1994; Celik et al. 2008, ASME J. Fluids Eng. 130(7):078001), here applied on the axis a path-dependent model actually discretises
//   locator:  E_oed = E_oed^ref ((c cos(phi) + sigma_1 sin(phi))/(c cos(phi) + p_ref sin(phi)))^m (compression-positive sigma_1), integrated over the vertical stress range of each sweep; tests/corpus/kv-cst-002-hs-oedometer.k2d, swept on three axes independently -- mesh density 0.5/0.25/0.125 m, load increments 10/20/40/80/160 at a converged tolerance with the seating phase pinned, and the stress range walked over 50-100, 100-200 and 200-400 kPa -- and then the SAME calibrated material walked as a one-dimensional oedometer at the stress point (katai::core::hs_integrate, built through the registry entry the driver uses), where none of those axes exist
//   quantity: settlement of the oedometer top [m] on each sweep, and the same settlement computed at the stress point [m]
//   expected: the mesh contributes nothing (a weightless column has a uniform strain field, which is exact in the element space); the load path is now below its own noise, so no order can be computed from it and no GCI quoted; the deviation from the closed form is ONE-SIGNED and closes as the stress rises; and the boundary-value run reproduces the constitutive routine, which is what makes the remainder a model deviation rather than an FE error
//   band:     no band is published for this case, deliberately, and the measurements are why. Mesh: 85 -> 1105 nodes changes the answer by 5.6e-16 relative, which is round-off. Load path: 10 -> 160 increments moves it by 0.020% -- until 2026-08-20 this axis was worth 2.9 percentage points, and it was not the load path: the stress-point integrator held the stress-dependent moduli fixed at the state each INCREMENT started from, a first-order error in the increment that no substep tolerance can see (0.9.0 N-1; the same fixture's material point moved 4.5% over a 32x outer refinement at a fixed integration tolerance of 1e-6 before the fix, and 6e-5 after). What is left is not a convergent sequence, so no order may be computed from it. Tolerance: 1e-4 and 1e-6 agree to 0.043%. Deviation: -1.3335% over 50-100 kPa, -0.9873% over 100-200, -0.6088% over 200-400 -- one-signed, largest just below p_ref, closing as the stress rises. The same three ranges at the STRESS POINT, no FE at all: -1.3490%, -0.9518%, -0.9130%, the worst of the three 0.305 pp from its FE counterpart, so at most a third of a percentage point of the deviation can be the finite element method and the rest is the model's distance from the idealised power law. Backends: the case that used to split between PARDISO and Eigen (+0.6419% against +0.4583% over 100-200) now agrees to 15 significant figures on all three ranges. Ceiling: THE CEILING HAS BEEN LIFTED. Refining the seating phase to 160 increments at 1e-6 used to stop converging altogether rather than getting better -- the tolerated error is an absolute residual, so shrinking the increment does not shrink what each one must achieve, and on a confining stress starting near zero the increments that could not achieve it were the ones whose stress points sat furthest from their own material law. Requiring those points to settle (0.9.0 N-2) carries the refined path through: 0.018640386 m, 0.0360% from the file's own 40 increments. The check that pinned the ceiling is now the one the study needed all along and could not ask while the run refused -- refining the seating path fourfold must not move the answer

// verify: KV-NUM-010
//   oracle:   closed_form
//   source:   Terzaghi (1943) one-dimensional consolidation, degree-of-consolidation series U(Tv) = 1 - sum_j (2/M^2) exp(-M^2 Tv), M = (2j+1) pi/2, Tv = cv t / H_dr^2, cv = k Eoed / gamma_w -- the same closed form KV-CON-002 is measured against; the estimator is the solution-verification procedure of Roache (1994) and Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, applied on the TIME axis, whose order is fixed in advance by the scheme: katai/analysis/consolidation.hpp integrates fully implicitly (alpha = 1), i.e. backward Euler, whose global error is O(dt)
//   locator:  tests/corpus/kv-con-002-terzaghi-column.k2d with the phase ended at Tv = 0.5, swept over 30/60/120/240 time steps at the file's own 0.4 m tri6 mesh, and over 0.8/0.4/0.2 m at a pinned 240 steps
//   quantity: degree of consolidation U at Tv = 0.5 [-]; the observed order of the time refinement [-]; and the GCI band on the file's own 120 steps [%]
//   expected: the observed order is 1, because backward Euler is first order and that is decided before the run; it is STABLE across overlapping triplets, unlike the load path of KV-NUM-009; Richardson extrapolation of three time steps recovers Terzaghi's closed form; and the band contains the true error, which can be checked here rather than trusted because the exact answer is known
//   band:     +/- 0.2480% on the file's own 120 steps (GCI at the observed order, Fs = 1.25), and the actual error there is 0.1965%, so the band contains it. Observed order 0.9850 from 120/60/30 and 0.9924 from 240/120/60 -- the two agree to 0.007, which is what an asymptotic range looks like and is exactly what KV-NUM-009's load path cannot produce -- that axis is now below its own noise floor and has no order to quote at all. Richardson U(dt->0) = 0.763961996 against the series' 0.763950331, +0.00153%. The mesh, again, is not the axis: 169 -> 1884 nodes moves U by 1.7e-8

// verify: KV-NUM-011
//   oracle:   closed_form
//   source:   1D SH site response of a damped elastic shear column on a rigid base at fundamental-mode resonance (Kramer 1996, Geotechnical Earthquake Engineering, ch. 7): |u_surf| = (4/pi) A/(w_1^2 2 xi), w_1 = 2 pi Vs/(4H) -- the same closed form KV-DYN-002 is measured against; the estimator is Roache (1994) / Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001, applied on the TIME axis. The order Newmark alone justifies is 2: katai/analysis/dynamics.hpp integrates with gamma = 1/2, beta = 1/4 (average acceleration, no numerical damping)
//   locator:  tests/corpus/kv-dyn-002-resonant-column.k2d, swept on two axes -- the DURATION over 10/20/40 cycles at a fixed 0.0025 s step, and the time step over 800/1600/3200 steps at 40 cycles where the resonant buildup is finished -- plus two mesh densities
//   quantity: peak surface displacement |u_surf| of the resonant phase [m], and the observed order of the time refinement [-]
//   expected: the duration is an axis of its own, because a resonant amplitude is built up rather than imposed; the time axis is a proper refinement axis and extrapolates onto the closed form; and its observed order is NOT the scheme's 2
//   band:     the prediction FAILED and that is the finding. Observed order 3.147 at 40 cycles, 3.147 and 2.998 at 160 cycles, and 3.097 / 2.795 on a grid deliberately INCOMMENSURATE with the period (T/dt = 18.43, 36.86, 73.72, 147.45), which eliminates sampling phase-lock as the explanation -- the order is about 3, reproducibly, where average-acceleration Newmark alone justifies 2. What the quantity is explains it: a PEAK of a resonant response, not the response at an instant, and a derived quantity need not inherit its scheme's order. So "the algebra fixes the order before the run", true for KV-STR-003's q h^2/12 and for KV-NUM-010's backward Euler, is NOT a general rule, and this is the case that bounds it. The axis is nonetheless clean: Richardson lands on the closed form to -0.036% at 40 cycles and -0.031% at 160, from two independent grids agreeing to 1e-5. Duration: the shipped 20 cycles is 0.21% short of steady state (envelope 1 - exp(-xi w t) = 0.9981) and 40 cycles is converged, so at the file's own settings the buildup is worth about as much as the time step. Mesh: 459 -> 6389 nodes moves |u| by 4.5e-6

#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/registry.hpp>
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
                  const katai::app::NumericalControls& nc, bool* ok_out,
                  katai::core::NewtonResult::Convergence* conv_out = nullptr,
                  int* iters_out = nullptr) {
    const auto res = katai::app::solve_phases(
        pr, mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, nullptr, nc);
    *ok_out = res.size() == 2 && res[1].ok;
    if (res.size() == 2) {
        if (conv_out) *conv_out = res[1].convergence;
        if (iters_out) *iters_out = res[1].iterations;
    }
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
m::Project all_phases(const m::Project& base, double tol, int steps, int iters,
                      double substol = 0.0) {
    m::Project pr = base;
    const auto set = [&](m::Phase& ph) {
        if (tol > 0.0) ph.tolerance = tol;
        if (steps > 0) ph.load_steps = steps;
        if (iters > 0) ph.max_iterations = iters;
        if (substol > 0.0) ph.substep_tolerance = substol;
    };
    set(pr.initial);
    for (auto& ph : pr.phases) set(ph);
    return pr;
}

// THE FOUR CASES IN THIS FILE RUN SEPARATELY, and the reason is wall clock rather than tidiness.
// With `ctest -j 6` a suite cannot finish faster than its single longest test, and this file was
// that test by a wide margin -- measured 2683 s on 2026-08-26, which IS the whole suite's 2683 s.
// The four verification cases share nothing but a corpus file, so they are four ctest entries
// over one executable: `test_phase_numerics <008|009|010|011>`, no argument meaning all four (the
// way it ran before, kept so a developer can still get the whole story in one run).
//
// KV-NUM-008 and KV-NUM-009 share the oedometer fixture below -- the same corpus case, its mesh,
// and the settlement the untouched file produces. Building it twice costs one extra solve, which
// is the price of the two running side by side.
struct Oedometer {
    m::Project base;
    katai::mesh::Mesh mesh;
    double u_default = 0.0;
    bool ok = false;
};

Oedometer oedometer_fixture() {
    Oedometer O;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
    if (!m::load_project(path, O.base, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        ++g_failures;
        return O;
    }
    check(!O.base.phases.empty(), "the case has a staged phase to control");
    if (O.base.phases.empty()) return O;

    const auto M = katai::app::mesh_from_project(O.base);
    check(M.ok, "oedometer mesh built from the case file");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return O; }
    O.mesh = M.mesh;

    bool solved = false;
    O.u_default = settlement(O.base, O.mesh, {}, &solved);
    check(solved, "the untouched case solves");
    if (!solved) return O;
    std::printf("  default (the material class chooses)      settlement = %.9f m\n", O.u_default);
    O.ok = true;
    return O;
}

void case_008(const Oedometer& O) {
    std::printf("== numerical controls carried by the file (KV-NUM-008) ==\n");
    if (!O.ok) return;
    const m::Project& base = O.base;
    const katai::mesh::Mesh& mesh = O.mesh;
    const double u_default = O.u_default;
    std::string err;
    bool ok = false;
    (void)u_default;

    // --- 1. The tolerated error ------------------------------------------------------------
    // This control used to be worth 0.59% on this problem between the Hardening Soil default of
    // 1e-2 and a converged 1e-6, and the check below leaned on that margin. It is worth 0.0009%
    // now: the local convergence criteria bind (0.9.0 N-2), so the default run already stands on
    // the converged answer. The margin is gone; what replaced it is below.
    const m::Project from_file = all_phases(base, 1e-6, 0, 0);
    const double u_file = settlement(from_file, mesh, {}, &ok);
    check(ok, "the case with a tolerated error in the FILE solves");
    katai::app::NumericalControls seam;
    seam.tolerance = 1e-6;
    const double u_seam = settlement(base, mesh, seam, &ok);
    check(ok, "the same tolerance through the seam solves");
    std::printf("  tolerance 1e-6 from the file              settlement = %.9f m\n", u_file);
    std::printf("  tolerance 1e-6 through the seam           settlement = %.9f m\n", u_seam);
    std::printf("  difference between the two routes         %.3e m\n", std::fabs(u_file - u_seam));
    check(u_file == u_seam, "the file and the seam are the same control, bit for bit");
    std::printf("  difference from the default run           %.3e m (%.2f%%)\n",
                std::fabs(u_file - u_default), 100.0 * std::fabs(u_file - u_default) / u_default);
    // The proof that the control is READ used to be that the answer moved. It no longer does:
    // with the local convergence criteria binding (0.9.0 N-2) the default run already stands on
    // the converged answer, and four decades of extra tolerance move it by 1.6e-7 m -- 0.0009%.
    // That is the second time on this case that a guard proving a control is read BY POINTING AT
    // A DIFFERENCE stopped proving anything the moment the difference was the defect being fixed
    // (the iteration limit below is the first).
    //
    // So it is proved where the control actually lands instead of where its consequence used to
    // show: the solve REPORTS the tolerance it ran under, and it spends more iterations reaching
    // it. Both are properties of the run rather than of how wrong the looser run happened to be,
    // so neither decays the next time the solver gets better.
    bool ok_c = false;
    katai::core::NewtonResult::Convergence conv_file, conv_default;
    int it_file = 0, it_default = 0;
    settlement(from_file, mesh, {}, &ok_c, &conv_file, &it_file);
    settlement(base, mesh, {}, &ok_c, &conv_default, &it_default);
    std::printf("  tolerance reported back: file %.3e, default %.3e; iterations %d vs %d\n",
                conv_file.tolerated, conv_default.tolerated, it_file, it_default);
    check(conv_file.tolerated == 1e-6 && conv_default.tolerated > 1e-6,
          "the run reports the tolerance it ran under, and the file's is the one it used");
    check(conv_file.force_error <= 1e-6,
          "...and it actually MET that tolerance, so the number is not merely carried");
    check(it_file > it_default,
          "...and reaching it cost iterations, so the control changed the calculation");

    // PER PHASE means per phase. Tightening only the staged step, and leaving the initial phase
    // to the material class, is a different calculation from tightening both -- and being able
    // to say which is the reason the controls live on the phase rather than on the project.
    m::Project one_phase = base;
    one_phase.phases[0].tolerance = 1e-6;
    const double u_one = settlement(one_phase, mesh, {}, &ok);
    check(ok, "the case with the tolerance on ONE phase solves");
    std::printf("  tolerance 1e-6 on the staged phase only   settlement = %.9f m\n", u_one);
    check(u_one != u_file,
          "and it differs from tightening every phase: the control really is per phase");

    // --- 2. The load increments -------------------------------------------------------------
    // A Hardening Soil phase is given 40 increments by default. Four is a coarser path to the
    // same load, and a path-dependent model does not arrive at quite the same place -- which is
    // the point: if the number were ignored, it would.
    const m::Project steps_file = all_phases(base, 0.0, 4, 0);
    const double u_steps_file = settlement(steps_file, mesh, {}, &ok);
    check(ok, "the case with 4 load increments in the FILE solves");
    katai::app::NumericalControls steps_seam;
    steps_seam.steps = 4;
    const double u_steps_seam = settlement(base, mesh, steps_seam, &ok);
    check(ok, "the same 4 increments through the seam solve");
    std::printf("  4 load increments: file %.9f m, seam %.9f m (default has 40: %.9f m)\n",
                u_steps_file, u_steps_seam, u_default);
    check(u_steps_file == u_steps_seam, "load increments: file == seam, bit for bit");
    check(u_steps_file != u_default, "and 4 increments is not 40, so the count is read");

    // --- 3. The iteration limit -------------------------------------------------------------
    // This control is read where it BITES, and a converged run is one it stops biting on. Until
    // 2026-08-20 the proof here was that a limit of 8 changed the answer; with the constitutive
    // integration under an error tolerance (0.9.0 N-1) it no longer does -- every increment of
    // this case now converges in three iterations, so 8, 12 and no limit at all agree bit for
    // bit. A guard that proves a control is read BY POINTING AT A DIFFERENCE stops proving
    // anything the moment the difference is the defect being fixed.
    //
    // So the pair is measured AT THE THRESHOLD instead, which proves more than the old check
    // did: one below it the run REFUSES and names the budget as the reason (the KV-NUM-012
    // contract -- a fraction of the load that is an iteration limit must not be reported as a
    // capacity), and at it the run converges to the default answer bit for bit. One number
    // apart, two different outcomes: the control reaches the solver, and above the threshold it
    // does not move the answer, which is what a converged run should do.
    //
    // The threshold itself is now FOUND rather than written down. It was 3 until the local
    // convergence criteria began to bind (0.9.0 N-2) -- the second time a number pinned in this
    // block moved because the solver got better. What is asserted is that a threshold exists,
    // that it is small enough to still be a threshold and not a budget, and that the two
    // outcomes straddle it; the measured value is printed rather than pinned.
    // HOW the threshold is found matters for what this test COSTS. Walking limit = 2, 3, 4 ... and
    // then again from the threshold up to 96, each step a full two-phase solve, is up to a hundred
    // solves and it made this the slowest test in the suite by a factor of two -- 2597 s of a
    // 2597 s wall clock, so the whole suite was as long as this one file. The search is kept but
    // it is no longer the FIRST thing tried: the recorded value is CHECKED, which costs two solves
    // and proves exactly the same straddle, and the scan runs only if the recorded value has
    // stopped being the boundary. Then it reports the new one loudly. A guard that cannot decay
    // silently was the whole point; paying a hundred solves for it every run was not.
    constexpr int kRecordedThreshold = 5;   // 2026-08-24; it was 3 before the local criteria bound
    const auto converges_at = [&](int limit, double* u_out) {
        const double u = settlement(all_phases(base, 0.0, 0, limit), mesh, {}, &ok);
        if (u_out) *u_out = u;
        return ok;
    };
    int threshold = 0;
    double u_at_threshold = 0.0;
    const bool at_recorded = converges_at(kRecordedThreshold, &u_at_threshold);
    const bool below_recorded = converges_at(kRecordedThreshold - 1, nullptr);
    if (at_recorded && !below_recorded) {
        threshold = kRecordedThreshold;
    } else {
        std::printf("  NOTE: the recorded iteration-limit threshold %d no longer straddles the "
                    "boundary (converges at it: %s, one below: %s) -- scanning for the new one\n",
                    kRecordedThreshold, at_recorded ? "yes" : "no",
                    below_recorded ? "yes" : "no");
        for (int limit = 2; limit <= 24 && threshold == 0; ++limit)
            if (converges_at(limit, &u_at_threshold)) threshold = limit;
    }
    std::printf("  iteration-limit threshold: %d\n", threshold);
    check(threshold >= 3 && threshold <= 24,
          "there is an iteration limit below which this case cannot converge, and it is a "
          "threshold rather than a budget");
    if (threshold == 0) return;

    // One below the threshold refuses. When the recorded value held, this is already known from
    // the probe above and is not re-solved; only a moved threshold pays for it again.
    if (threshold != kRecordedThreshold) {
        settlement(all_phases(base, 0.0, 0, threshold - 1), mesh, {}, &ok);
    } else {
        ok = below_recorded;
    }
    check(!ok, "one below the threshold, the FILE's limit makes the run refuse rather than drift");
    katai::app::NumericalControls starved_seam;
    starved_seam.max_iterations = threshold - 1;
    settlement(base, mesh, starved_seam, &ok);
    check(!ok, "and the same limit through the seam refuses too: both routes reach the solver");

    // The FILE's run at the threshold is the probe above -- solving the same calculation twice
    // would prove nothing the first solve did not.
    const double u_iter_file = u_at_threshold;
    check(threshold > 0 && u_iter_file > 0.0, "at the threshold the case in the FILE solves");
    katai::app::NumericalControls iter_seam;
    iter_seam.max_iterations = threshold;
    const double u_iter_seam = settlement(base, mesh, iter_seam, &ok);
    check(ok, "the same limit through the seam solves");
    std::printf("  at the threshold:  file %.9f m, seam %.9f m (default %.9f m)\n",
                u_iter_file, u_iter_seam, u_default);
    check(u_iter_file == u_iter_seam, "iteration limit: file == seam, bit for bit");

    // There are TWO thresholds here and the old check conflated them, which only became visible
    // once the first one moved. At `threshold` the run stops REFUSING -- but it survives by
    // cutting increments back, so it walks a different load path and lands 0.2% away. The answer
    // becomes the default's only at a HIGHER limit, where no increment is cut back at all. Both
    // are found, and what is asserted is the ordering and the fact that the second exists: a
    // limit high enough is a limit that does not enter the answer, which is the property that
    // makes it a patience setting rather than an accuracy one.
    // Recorded-then-searched, for the same reason and at the same saving as the threshold above:
    // this scan could run to 96 full solves on its own.
    constexpr int kRecordedSettled = 6;   // 2026-08-24, one above the threshold
    const auto reproduces_default = [&](int limit) {
        double u = 0.0;
        return converges_at(limit, &u) && u == u_default;
    };
    int settled = 0;
    if (reproduces_default(kRecordedSettled) &&
        (kRecordedSettled <= threshold || !reproduces_default(kRecordedSettled - 1))) {
        settled = kRecordedSettled;
    } else {
        std::printf("  NOTE: the recorded settling limit %d is no longer the first that reproduces "
                    "the default -- scanning for the new one\n", kRecordedSettled);
        for (int limit = threshold; limit <= 96 && settled == 0; ++limit)
            if (reproduces_default(limit)) settled = limit;
    }
    std::printf("  refuses below %d; converges at %d (%.4f%% from the default, by cutting back); "
                "reproduces the default from %d\n",
                threshold, threshold, 100.0 * std::fabs(u_iter_file - u_default) / u_default,
                settled);
    check(settled >= threshold && settled <= 96,
          "a high enough iteration limit stops entering the answer at all: above it the run "
          "reproduces the default bit for bit");
    check(std::fabs(u_at_threshold - u_default) / u_default < 0.01,
          "and between the two thresholds the run is cutting back rather than drifting: it "
          "lands close, not anywhere");

    // --- 3b. The FOURTH control: the constitutive integration tolerance ----------------------
    // The three above govern the equilibrium iteration. This one governs the material law: how
    // accurately a stress point is walked along it INSIDE an increment. It was an environment
    // variable until .k2d v15, which meant the half of a published claim that says "and this is
    // how the constitutive path was integrated" could not be written down at all.
    //
    // It is asserted the same hostile way as the others -- file == seam, and both different from
    // the default -- because the failure this whole case exists to catch is a control that is
    // read, stored, echoed and then not used. Here that failure would be especially quiet: a
    // dropped integration tolerance produces a perfectly convergent run with a slightly wrong
    // stress path, which no residual reports and no plot shows.
    // ALL THREE RUNS USE FOUR LOAD INCREMENTS instead of the file's forty, and that is a
    // measurement, not a shortcut. What is asserted here is that a control REACHES the solver --
    // identity between two routes and a difference from the default -- and neither claim needs
    // the converged settlement. The full-increment version of this block cost 912 s, because a
    // LOOSER integration is expensive on this case: the equilibrium iteration grinds against the
    // integration noise it can no longer resolve (the same effect that set the default at 1e-5 in
    // the first place). Paying fifteen minutes per run to re-learn that would be a poor trade for
    // a check about plumbing.
    constexpr int kSubSteps = 4;
    const double u_sub_base = settlement(all_phases(base, 0.0, kSubSteps, 0), mesh, {}, &ok);
    check(ok, "the four-increment baseline solves");
    const m::Project sub_file = all_phases(base, 0.0, kSubSteps, 0, 1.0e-3);
    const double u_sub_file = settlement(sub_file, mesh, {}, &ok);
    check(ok, "the case with a substep tolerance in the FILE solves");
    katai::app::NumericalControls sub_seam;
    sub_seam.substep_tolerance = 1.0e-3;
    const double u_sub_seam = settlement(all_phases(base, 0.0, kSubSteps, 0), mesh, sub_seam, &ok);
    check(ok, "the same substep tolerance through the seam solves");
    std::printf("  substep tolerance 1e-3 (%d increments):  file %.9f m, seam %.9f m "
                "(default 1e-5 %.9f m, apart by %.4f%%)\n",
                kSubSteps, u_sub_file, u_sub_seam, u_sub_base,
                100.0 * std::fabs(u_sub_file - u_sub_base) / u_sub_base);
    check(u_sub_file == u_sub_seam,
          "substep tolerance: file == seam, bit for bit -- one control, two routes");
    check(std::fabs(u_sub_file - u_sub_base) / u_sub_base > 1e-4,
          "and a looser constitutive integration changes the answer, so the file's value is "
          "genuinely reaching the material routine rather than the environment's default");

    // --- 4. Precedence, stated and tested ---------------------------------------------------
    // When both are set the SEAM wins. The seam exists so that a given file can be re-run at
    // other numerics (that is what KV-NUM-007 does to a checked-in corpus case); a file that
    // always won would make that study impossible to perform.
    const m::Project both = all_phases(base, 1e-2, 0, 0);   // the file asks for the loose default
    katai::app::NumericalControls tight;
    tight.tolerance = 1e-6;               // the caller asks for a converged run
    const double u_both = settlement(both, mesh, tight, &ok);
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
    written.phases[0].substep_tolerance = 2.5e-6;   // .k2d v15, the fourth control
    written.initial.tolerance = 1e-7;     // the initial phase carries them too
    const std::string text = m::project_to_json(written);
    check(m::project_from_json(text, trip, &err, nullptr), "the written project reads back");
    check(trip.phases[0].tolerance == 3.5e-4 && trip.phases[0].load_steps == 17 &&
              trip.phases[0].max_iterations == 33 && trip.initial.tolerance == 1e-7 &&
              trip.phases[0].substep_tolerance == 2.5e-6,
          "all four controls survive the round trip, on the initial phase as well");
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

}

// --- 6. WHERE THE NUMERICAL ERROR OF THIS CASE ACTUALLY LIVES (KV-NUM-009) -------------------
// KV-NUM-008 shows the controls are read. This asks the question underneath it: of the +1% this
// case sits from its closed form, how much is the mesh, how much the stopping rule, how much the
// load path -- and how much is not numerical at all. Until KV-NUM-009 the corpus had only ever
// banded the MESH, which for a path-dependent model is the wrong axis to start on.
void case_009(const Oedometer& O) {
    std::printf("\n== where the numerical error lives (KV-NUM-009) ==\n");
    if (!O.ok) return;
    const m::Project& base = O.base;
    const katai::mesh::Mesh& mesh = O.mesh;
    const double u_default = O.u_default;
    std::string err;
    bool ok = false;
    (void)err;

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

    // (b) THE LOAD PATH WAS THE DOMINANT DISCRETISATION, AND IS NOT ANY MORE. Until 2026-08-20
    //     a 16x refinement of the increments moved this answer by 2.9 percentage points, and this
    //     block asserted that it did -- the axis had to be real for the rest of the section to
    //     mean anything. It was real, and it was not the load path: the stress-point integrator
    //     held the stress-dependent moduli (E_ur, E_i, q_a, q_f) fixed at the state each
    //     INCREMENT started from, which is a first-order error in the increment size that no
    //     substep tolerance can see. Reading them at each substep instead (0.9.0 N-1) collapses
    //     the same sweep to 0.020%, a factor of 145.
    //     What that leaves is an axis below its own noise: the five settlements no longer form a
    //     convergent sequence at all, so no order can be computed from them and no GCI may be
    //     quoted -- the same refusal as before, on a better reason. The measured spread is
    //     reported, and it is now two orders of magnitude smaller than the model deviation (c)
    //     measures, which is what makes that deviation attributable.
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
        u_path[i] = settlement(pr, mesh, {}, &ok);
        if (!ok) { path_ok = false; break; }
        std::printf("  %4d load increments  settlement = %.12e\n", path_steps[i], u_path[i]);
    }
    check(path_ok, "the case solves at every step count in the sweep");
    if (path_ok) {
        const double spread = std::fabs(u_path[0] - u_path[4]) / u_path[4];
        std::printf("  spread over a 16x refinement: %.3f%% (the mesh gave %.1e)\n",
                    100.0 * spread, mesh_ok ? std::fabs(u_mesh[2] - u_mesh[0]) / u_mesh[0] : 0.0);
        check(spread < 1e-3,
              "the load path is no longer a percentage-point axis: a 16x refinement moves it "
              "less than 0.1%");

        // Is what is left a convergent sequence, or noise? The estimator needs
        // phi(d) = phi_exact + C d^p, which requires successive differences of one sign and a
        // stable ratio. Count how many of the three overlapping triplets even keep their sign.
        int monotone = 0;
        for (int i = 0; i < 3; ++i) {
            const double e21 = u_path[i + 1] - u_path[i + 2], e32 = u_path[i] - u_path[i + 1];
            if (e21 != 0.0 && e32 / e21 > 0.0) ++monotone;
        }
        std::printf("  of the three overlapping triplets, %d keep a consistent sign\n", monotone);
        check(monotone < 3,
              "the sequence is inside its own noise, so no order may be computed and no GCI "
              "quoted from the load path");
    }

    // (c) WHAT SURVIVES IS THE MODEL, AND IT IS SHOWN BY LEAVING THE FE OUT. The residual left
    //     after an exact mesh, a converged tolerance and a 16x refined path is about one per
    //     cent, one-signed, and largest just below p_ref.
    //
    //     Until 2026-08-20 this block argued from a SIGNATURE -- the deviation grew with stress
    //     level and changed sign near p_ref, therefore it was a cap calibrated at p_ref. That
    //     signature is gone, because it was the frozen-modulus error described in (b) rather than
    //     the model: the deviation does not change sign, and its magnitude FALLS as the stress
    //     rises above p_ref. An argument from a shape is only as good as the shape.
    //
    //     So the claim is made the other way now, by measuring the same law where none of the
    //     numerics under discussion exist. The SAME calibrated material -- built through the same
    //     registry entry the driver uses, so alpha and beta are the run's own -- is walked as a
    //     one-dimensional oedometer at the STRESS POINT: no mesh, no load path, no equilibrium
    //     iteration, no linear solver. If the boundary-value run and the constitutive routine
    //     land in the same place, then whatever separates BOTH of them from the closed form is
    //     the model's distance from the idealised power law, and the FE is not on trial for it.
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
        const double got = settlement(pr, mesh, {}, &ok);
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

    // The same law at the stress point. hs_integrate is the routine the FE calls at every Gauss
    // point; here it is called directly, on the material the driver would build from this file.
    const katai::core::ModelEntry* hs_entry = katai::core::find_model("HardeningSoil");
    check(hs_entry != nullptr && !base.materials.empty(),
          "the case's material resolves through the constitutive registry");
    double dev_pt[3] = {0, 0, 0};
    bool point_ok = hs_entry != nullptr && !base.materials.empty();
    if (point_ok) {
        const katai::core::MaterialModel mm =
            hs_entry->build(katai::app::to_material_params(base.materials[0]));
        const double kH = 4.0, kEoed = 30000.0, kPref = 100.0, kM = 0.5;
        const double de = 2.0e-5;
        for (int i = 0; i < 3 && point_ok; ++i) {
            const double sa = ranges[i][0], sb = ranges[i][1];
            Eigen::Vector3d sig(0.5, 0.5, 0.5);          // primary loading from near zero, as the
            double gp = 0.0, pp = 0.5, eps = 0.0;        // column itself is loaded from rest
            double eps_a = 0.0, eps_b = 0.0;
            bool got_a = false, got_b = false;
            for (int k = 0; k < 400000 && !got_b; ++k) {
                const katai::core::HsIntegrated r =
                    katai::core::hs_integrate(mm.hs, sig, gp, pp, Eigen::Vector3d(de, 0.0, 0.0));
                if (!r.stress.allFinite()) break;
                const double s0 = sig(0), s1 = r.stress(0);
                if (s1 > s0) {   // linear in strain between samples, so the mark is not a step artefact
                    if (!got_a && s1 >= sa) { eps_a = eps + de * (sa - s0) / (s1 - s0); got_a = true; }
                    if (got_a && !got_b && s1 >= sb) {
                        eps_b = eps + de * (sb - s0) / (s1 - s0); got_b = true;
                    }
                }
                sig = r.stress; gp = r.gamma_p; pp = r.pp; eps += de;
            }
            if (!got_b) { point_ok = false; break; }
            const double got = kH * (eps_b - eps_a);
            const double want = std::pow(kPref, kM) / kEoed *
                                (std::pow(sb, 1.0 - kM) - std::pow(sa, 1.0 - kM)) / (1.0 - kM) * kH;
            dev_pt[i] = (got - want) / want;
            std::printf("  %3.0f -> %3.0f kPa  STRESS POINT %.9f  %+.4f%%   (FE %+.4f%%, "
                        "apart by %.3f pp)\n",
                        sa, sb, got, 100.0 * dev_pt[i], 100.0 * dev[i],
                        100.0 * std::fabs(dev[i] - dev_pt[i]));
        }
    }
    check(point_ok, "the same law integrates at the stress point over each range");

    if (range_ok && point_ok) {
        check(dev[0] < 0.0 && dev[1] < 0.0 && dev[2] < 0.0 &&
              dev_pt[0] < 0.0 && dev_pt[1] < 0.0 && dev_pt[2] < 0.0,
              "the deviation is ONE-SIGNED on both routes: the model is stiffer than the "
              "idealised law at every stress level, and does not cross at p_ref");
        check(dev[0] < dev[1] && dev[1] < dev[2],
              "and the gap CLOSES as the stress rises above p_ref, rather than growing");
        double worst = 0.0;
        for (int i = 0; i < 3; ++i) worst = std::fmax(worst, std::fabs(dev[i] - dev_pt[i]));
        std::printf("  FE against the stress point, worst of the three ranges: %.3f pp\n",
                    100.0 * worst);
        check(worst < 5e-3,
              "the boundary-value run reproduces the constitutive routine to within half a "
              "percentage point, so at most that much of the deviation can be the FE");
    }

    // (d) THE CEILING ON PATH REFINEMENT IS STILL THERE, AND ASKING ONE POINT SAID OTHERWISE.
    //     This block used to refine the SEATING phase to 160 increments and assert that it
    //     converged, "where it used to refuse" -- reading that success as proof that the local
    //     convergence criteria (0.9.0 N-2) had lifted a ceiling. That reading came from a single
    //     point on an axis, and it did not survive its own neighbours. Measured 2026-08-27 by
    //     running the case file itself through the CLI on both linear-solver backends:
    //
    //       seating steps    Eigen                    PARDISO
    //          40            ok, 3.093e-08            ok, 3.093e-08     <- identical, every digit
    //          80            ok, 1.478e-08            ok, 1.478e-08     <- identical, every digit
    //         120            REFUSED at 74% of load   REFUSED at 83% of load
    //         160            REFUSED at 68% of load   ok, 6.570e-07 (of 1e-6), 390 s
    //
    //     BOTH backends refuse at 120. The ceiling was not lifted; it sits between 80 and 120,
    //     and 160 is the one rung where PARDISO gets through -- non-monotonically, a coarser path
    //     failing where a finer one passes, at 66% of the tolerance and for twenty times the run
    //     time of a rung that converges thirty times further inside it. Nothing about that is a
    //     property of this program: below the ceiling the two backends agree to every printed
    //     digit, and above it the answer is decided by round-off. The portable composition is
    //     what said so -- this is what §11.7 exists for, and the first thing it caught.
    //
    //     SO THE CHECK ASKS THE QUESTION THE STUDY ACTUALLY NEEDED, at a rung where the answer
    //     exists on both compositions: refining the seating path does not move the answer. That
    //     is path independence, which is what (b) pinning the count at 40 was standing in for.
    //
    //     WHAT THIS DOES NOT ASSERT, DELIBERATELY: that 120 refuses. Asserting a refusal is the
    //     trap this block already fell into once -- the old sentence had to end "if this check
    //     ever fails because the run now succeeds, the ceiling has moved", and then a single
    //     backend moved it and the sentence was rewritten the wrong way. The asymmetry is
    //     honest: if the ceiling drops below 80 this check fails, and if it rises the table
    //     above goes stale without failing anything. Re-measuring it is four CLI runs on the
    //     case file with `initial.loadsteps` set, and it belongs to N-3, whose subject is
    //     choosing the step size rather than being handed one.
    m::Project seat_fine = base;
    seat_fine.initial.load_steps = 80;
    seat_fine.initial.tolerance = 1e-6;
    seat_fine.initial.max_iterations = 500;
    const double u_seat_fine = settlement(seat_fine, mesh, {}, &ok);
    check(ok, "the seating phase converges at 80 increments, the last rung below the ceiling");
    if (ok) {
        const double drift = std::fabs(u_seat_fine - u_default) / u_default;
        std::printf("  seating phase at 80 increments, tol 1e-6: %.9f m (%.4f%% from the "
                    "file's own count)\n", u_seat_fine, 100.0 * drift);
        check(drift < 0.01,
              "and refining the seating path does not move the answer: the seating "
              "count is not what sets it");
    }

}

// --- 7. THE SAME QUESTION WHERE THE AXIS *IS* CLEAN (KV-NUM-010) -----------------------------
void case_010() {
    std::string err;
    bool ok = false;
    (void)ok;
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
        ++g_failures;
        return;
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

}

// --- 8. WHERE "THE ALGEBRA DECIDES THE ORDER" STOPS BEING TRUE (KV-NUM-011) ------------------
void case_011() {
    std::string err;
    bool ok = false;
    (void)ok;
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
        ++g_failures;
        return;
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

}

}  // namespace

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "all";
    if (which == "008" || which == "009" || which == "all") {
        const Oedometer O = oedometer_fixture();
        if (which != "009") case_008(O);
        if (which != "008") case_009(O);
    }
    if (which == "010" || which == "all") case_010();
    if (which == "011" || which == "all") case_011();
    if (which != "all" && which != "008" && which != "009" && which != "010" && which != "011") {
        std::printf("FAIL: unknown case '%s' (use 008, 009, 010, 011 or nothing)\n",
                    which.c_str());
        return 1;
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
