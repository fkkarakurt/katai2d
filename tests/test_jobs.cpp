// The job object (katai/jobs/job.hpp) pinned end to end, using the checked-in corpus
// .k2d files as fixtures (they are the contract: a job runs what a file says):
//
//   (1) EQUIVALENCE -- the Job path reproduces the direct driver path BIT-identically
//       on the same thread (same mesh, same solve; the job layer may not change the
//       answer). This is the seed of Stage D's front-end-equivalence gate.
//   (2) PROGRESS -- the callback reports (index, total, name) for every phase, in
//       order, with no UI assumption.
//   (3) CANCELLATION -- cooperative; a cancel issued from the progress callback stops
//       the announced phase from running, the completed prefix is kept and the state is
//       Cancelled with an honest message. A request that arrives while a phase is being
//       solved stops it at its next step or iteration (a Plastic and a Safety phase):
//       the phase is dropped whole, a Safety phase is not reported as a collapse, and
//       the next run on the same thread is not touched by it. A stress-point walk polls
//       too, between its substeps: a Hardening Soil cap calibration stopped half way
//       leaves no calibration behind -- the material, calibrated again, meets its targets.
//   (4) TIMINGS -- one entry per completed phase, named, non-negative, sum > 0.
//   (4a) LIVE -- what the job says while it runs, read from another thread: the steps only grow, and
//       the last phase's curve ends where its result ends.
//   (5) CONCURRENCY -- two jobs on two threads in one process both finish Done and
//       match their serial runs (no globals, no static solver state). The comparison
//       uses a 1e-10 relative band rather than bit-equality: MKL's internal threading
//       may partition differently under concurrent load, which can move the last bits
//       without moving the answer.
//
// verify: none -- job-layer mechanics (equivalence, progress, cancellation, timings,
// concurrency); the physics NUMBERS behind these runs are declared by the corpus cases
// themselves (KV-FND-008, KV-EXC-001) and asserted in test_input_corpus.
#include <katai/jobs/job.hpp>
#include <katai/io/project_io.hpp>
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/interrupt.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <thread>
#include <vector>

namespace m = katai::model;
namespace jobs = katai::jobs;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

m::Project load_corpus(const char* file) {
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/" + file;
    if (!m::load_project(path, pr, &err))
        std::fprintf(stderr, "FATAL: cannot load %s (%s)\n", path.c_str(), err.c_str());
    return pr;
}

// ------------------------------------------------------------------- (1) equivalence --
void test_equivalence() {
    std::printf("-- job path == driver path, bit-identical (kv-fnd-008) --\n");
    const m::Project pr = load_corpus("kv-fnd-008-strip-load.k2d");

    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "direct path meshed");
    const auto direct = katai::app::solve_gravity_le(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    check(direct.ok, "direct path solved");

    jobs::Job job(pr);
    check(job.state() == jobs::JobState::Pending, "job starts Pending");
    check(job.run(), "job ran");
    check(job.state() == jobs::JobState::Done, "job state Done");
    check(job.results().size() == 1, "single-phase run produced one result");
    if (!direct.ok || job.results().size() != 1) return;

    const auto& R = job.results()[0];
    const double ddisp = (R.disp - direct.disp).lpNorm<Eigen::Infinity>();
    std::printf("   max|disp(job) - disp(direct)| = %.3e\n", ddisp);
    check(ddisp == 0.0, "displacements bit-identical");
    check(R.max_disp == direct.max_disp, "max_disp bit-identical");
}

// ------------------------------------------------------- (2) progress + (4) timings --
void test_progress_and_timings() {
    std::printf("-- progress callback + per-phase timings (kv-exc-001) --\n");
    jobs::Job job(load_corpus("kv-exc-001-staged-excavation.k2d"));
    std::vector<std::string> seen;
    job.set_on_phase([&seen](int cur, int total, const std::string& name) {
        seen.push_back(std::to_string(cur) + "/" + std::to_string(total) + " " + name);
    });
    check(job.run(), "staged job ran");
    check(seen.size() == 2, "callback fired once per phase");
    for (const auto& s : seen) std::printf("   announced: %s\n", s.c_str());
    check(seen.size() == 2 && seen[0] == "0/2 Initial phase" && seen[1] == "1/2 Excavate",
          "announcements carry (index, total, name) in order");

    const auto& t = job.timings();
    check(t.size() == 2, "one timing per completed phase");
    double sum = 0.0;
    for (const auto& pt : t) {
        std::printf("   timing: %-14s %.4f s\n", pt.name.c_str(), pt.seconds);
        check(pt.seconds >= 0.0, "timing non-negative");
        sum += pt.seconds;
    }
    check(sum > 0.0, "the run took measurable time");
}

// --------------------------------------------------------------------------- (4a) live --
void test_live_progress() {
    std::printf("-- what a running job says, read from another thread (kv-exc-001) --\n");
    jobs::Job job(load_corpus("kv-exc-001-staged-excavation.k2d"));
    check(job.live().step == 0 && job.live().curve.empty(), "nothing to say before the run");
    std::atomic<bool> done{false};
    // A new phase starts its count again, and nothing else may lower it. Falls are counted rather than
    // judged by the phase at the moment: the count resets just before the phase is announced, and a reader
    // under load can land between the two.
    int falls = 0;
    int last_step = 0;
    std::atomic<int> phase_seen{0};
    job.set_on_phase([&phase_seen](int, int, const std::string&) { ++phase_seen; });
    std::thread runner([&] { job.run(); done = true; });
    while (!done) {
        const auto now = job.live();
        if (now.step < last_step) ++falls;
        last_step = now.step;
        std::this_thread::yield();
    }
    runner.join();
    check(job.state() == jobs::JobState::Done, "the watched job finished");
    check(falls < phase_seen.load(), "the steps of a phase only grow");
    const auto end = job.live();
    std::printf("   %d steps, %zu curve points, factor %.4f, movement %.4f m\n", end.step,
                end.curve.size(), end.factor, end.movement);
    check(end.step > 0 && (size_t)end.step == end.curve.size(), "one curve point per accepted step");
    if (job.results().empty() || end.curve.empty()) return;
    const auto& R = job.results().back();
    double top = 0.0;
    for (const auto& p : end.curve) top = std::max(top, p[1]);
    check(std::fabs(top - R.load_factor) <= 1e-12, "the curve's highest factor is the factor the phase reports");
}

// ------------------------------------------------------------------ (3) cancellation --
void test_cancellation() {
    std::printf("-- cooperative cancellation at a phase boundary (kv-exc-001) --\n");
    jobs::Job job(load_corpus("kv-exc-001-staged-excavation.k2d"));
    job.set_on_phase([&job](int cur, int, const std::string&) {
        if (cur == 1) job.request_cancel();   // cancel the announced second phase
    });
    check(!job.run(), "cancelled run reports failure to complete");
    check(job.state() == jobs::JobState::Cancelled, "state is Cancelled");
    check(job.results().size() == 1, "the completed prefix (initial phase) is kept");
    check(!job.results().empty() && job.results()[0].ok, "the completed phase is intact");
    std::printf("   message: %s\n", job.message().c_str());
    check(job.message() == "cancelled after 1 of 2 phase(s)", "honest cancellation message");

    jobs::Job early(load_corpus("kv-fnd-008-strip-load.k2d"));
    early.request_cancel();
    check(!early.run() && early.state() == jobs::JobState::Cancelled &&
              early.results().empty(),
          "a cancel before run() stops everything");
}

// A request that arrives INSIDE a phase. The predicate counts its polls from the moment the
// phase is announced and answers yes on the fourth: the boundary poll right after the
// announcement is the first, so the stop can only come from a poll in the solver's own loops.
// Cancelled at phase boundaries alone, the phase would finish and be in the results.
void test_cancellation_inside_a_phase() {
    std::printf("-- cancellation inside a phase (kv-exc-001) --\n");
    m::Project pr = load_corpus("kv-exc-001-staged-excavation.k2d");
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "meshed");
    const auto init = katai::app::initial_phase_from(pr.initial_procedure);

    const auto stopped_in = [&](const m::Project& p, int phase, size_t& kept, bool& kept_ok) {
        int armed = -1, polls = 0;
        const auto on_phase = [&](int cur, int, const std::string&) { if (cur == phase) armed = cur; };
        const std::function<bool()> cancelled = [&] { return armed >= 0 && ++polls >= 4; };
        const auto out = katai::app::solve_phases(p, M.mesh, init, on_phase, cancelled);
        kept = out.size();
        kept_ok = true;
        for (const auto& r : out) kept_ok = kept_ok && r.ok;
        return polls;
    };

    size_t kept = 0;
    bool kept_ok = false;
    int polls = stopped_in(pr, 1, kept, kept_ok);
    std::printf("   Plastic: %d polls, %zu phase(s) kept\n", polls, kept);
    check(polls == 4 && kept == 1 && kept_ok,
          "a Plastic phase stops inside its steps; only the initial phase is kept");

    // The excavation's soil is elastic and has no strength to reduce: the Safety case runs on
    // Mohr-Coulomb ground (c = 10 kPa, phi = 30 deg), so that the strength reduction iterates.
    m::Project with_safety = pr;
    for (auto& mat : with_safety.materials) {
        mat.model = m::SoilModel::MohrCoulomb;
        mat.c = 10.0;
        mat.phi = 30.0;
    }
    m::Phase safety = pr.phases.back();
    safety.name = "Safety";
    safety.type = m::PhaseType::Safety;
    with_safety.phases.push_back(safety);
    polls = stopped_in(with_safety, 2, kept, kept_ok);
    std::printf("   Safety: %d polls, %zu phase(s) kept\n", polls, kept);
    check(polls == 4 && kept == 2 && kept_ok,
          "a Safety phase stops inside its strength reduction and is not reported as a collapse");

    const m::Project column = load_corpus("kv-con-002-terzaghi-column.k2d");
    const auto MC = katai::app::mesh_from_project(column);
    {
        int armed = -1, cpolls = 0;
        const auto on_phase = [&](int cur, int, const std::string&) { if (cur == 1) armed = cur; };
        const std::function<bool()> cancelled = [&] { return armed >= 0 && ++cpolls >= 4; };
        const auto out = katai::app::solve_phases(
            column, MC.mesh, katai::app::initial_phase_from(column.initial_procedure), on_phase,
            cancelled);
        std::printf("   Consolidation: %d polls, %zu phase(s) kept\n", cpolls, out.size());
        check(cpolls == 4 && out.size() == 1 && out[0].ok,
              "a consolidation phase stops inside its time steps (kv-con-002)");
    }

    const auto again = katai::app::solve_phases(with_safety, M.mesh, init);
    std::printf("   plain run: %zu phase(s), last ok %d (%s)\n", again.size(),
                again.empty() ? 0 : (int)again.back().ok,
                again.empty() ? "" : again.back().message.c_str());
    check(again.size() == 3 && again.back().ok,
          "the next run on the same thread runs to the end (the request did not outlive its run)");
}

// A cancel that reaches inside a Hardening Soil cap calibration (its walks poll the thread's
// interrupt hook, katai/materials/interrupt.hpp) must leave nothing in the calibration cache: the
// alpha and beta of a search stopped half way are not the material's, and every later build of
// it in the process read them back. A material no other test calibrates, so the first
// calibration is the stopped one.
void test_cancel_inside_calibration() {
    std::printf("-- a cancel inside a Hardening Soil cap calibration --\n");
    katai::core::HardeningSoilParams p;
    p.p_ref = 100.0; p.E50_ref = 3.17e4; p.Eur_ref = 9.51e4; p.Eoed_ref = 3.17e4;
    p.m = 0.5; p.nu_ur = 0.2; p.friction = 31.0 * 3.14159265358979323846 / 180.0; p.Rf = 0.9;
    const double k0nc = 1.0 - std::sin(p.friction);
    struct Stop {};
    int polls = 0;
    katai::core::interrupt_hook() = {[](void* n) { if (++*static_cast<int*>(n) > 50) throw Stop{}; }, &polls};
    katai::core::HardeningSoilParams half = p;
    bool stopped = false;
    try { katai::core::hs_calibrate_cap(half, k0nc); } catch (const Stop&) { stopped = true; }
    katai::core::interrupt_hook() = {};
    double e_half = 0.0, k0_half = 0.0;
    katai::core::hs_oedometer_probe(half, e_half, k0_half);
    katai::core::HardeningSoilParams full = p;
    katai::core::hs_calibrate_cap(full, k0nc);
    double e = 0.0, k0 = 0.0;
    katai::core::hs_oedometer_probe(full, e, k0);
    std::printf("   stopped after %d polls: K0 %.6f; calibrated again: K0 %.6f (target %.6f), Eoed %.2f (%.0f)\n",
                polls, k0_half, k0, k0nc, e, p.Eoed_ref);
    check(stopped && std::fabs(k0_half - k0nc) > 1e-3,
          "the calibration was stopped inside, at a state that is not the material's");
    check(std::fabs(k0 - k0nc) < 1e-4 && std::fabs(e - p.Eoed_ref) < 1e-3 * p.Eoed_ref,
          "calibrated again, the material meets its targets (nothing half done was kept)");
}

// ------------------------------------------------------- validator enforcement --
void test_validation_refusal() {
    std::printf("-- run() executes the input contract first --\n");
    m::Project bad = load_corpus("kv-fnd-008-strip-load.k2d");
    bad.mesh.order = 9;   // indefensible on sight
    jobs::Job job(bad);
    check(!job.run() && job.state() == jobs::JobState::Failed,
          "an invalid project is refused, not solved");
    check(!job.report().ok(), "the validation report carries the error");
    check(job.results().empty(), "nothing was solved");
    std::printf("   message: %s\n", job.message().c_str());
    check(job.message().find("mesh.order") != std::string::npos,
          "the refusal names the field path");
}

// ------------------------------------------------------------------- (5) concurrency --
void test_concurrency() {
    std::printf("-- two jobs concurrently in one process --\n");
    const m::Project pa = load_corpus("kv-fnd-008-strip-load.k2d");
    const m::Project pb = load_corpus("kv-exc-001-staged-excavation.k2d");

    jobs::Job serial_a(pa), serial_b(pb);
    check(serial_a.run() && serial_b.run(), "serial reference runs");

    jobs::Job conc_a(pa), conc_b(pb);
    std::thread ta([&conc_a] { conc_a.run(); });
    std::thread tb([&conc_b] { conc_b.run(); });
    ta.join();
    tb.join();
    check(conc_a.state() == jobs::JobState::Done && conc_b.state() == jobs::JobState::Done,
          "both concurrent jobs finished Done");
    if (conc_a.results().empty() || conc_b.results().empty()) return;

    const auto rel = [](double got, double ref) {
        return std::fabs(got - ref) / std::fmax(std::fabs(ref), 1e-30);
    };
    const double da = rel(conc_a.results()[0].max_disp, serial_a.results()[0].max_disp);
    const double db = rel(conc_b.results().back().max_disp, serial_b.results().back().max_disp);
    std::printf("   rel|max_disp conc - serial|: A = %.3e  B = %.3e\n", da, db);
    check(da < 1e-10 && db < 1e-10, "concurrent results match serial (1e-10 relative)");
}

}  // namespace

int main() {
    std::printf("The job object: headless runs with progress, cancellation, timings\n\n");
    test_equivalence();
    test_progress_and_timings();
    test_live_progress();
    test_cancellation();
    test_cancellation_inside_a_phase();
    test_cancel_inside_calibration();
    test_validation_refusal();
    test_concurrency();
    if (g_failures == 0) {
        std::printf("\nOK: the job layer runs what the file says and cannot change the answer\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
