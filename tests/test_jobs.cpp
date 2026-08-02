// The job object (katai/jobs/job.hpp) pinned end to end, using the checked-in corpus
// .k2d files as fixtures (they are the contract: a job runs what a file says):
//
//   (1) EQUIVALENCE -- the Job path reproduces the direct driver path BIT-identically
//       on the same thread (same mesh, same solve; the job layer may not change the
//       answer). This is the seed of Stage D's front-end-equivalence gate.
//   (2) PROGRESS -- the callback reports (index, total, name) for every phase, in
//       order, with no UI assumption.
//   (3) CANCELLATION -- cooperative, phase-granular; a cancel issued from the progress
//       callback stops the announced phase from running, the completed prefix is kept
//       and the state is Cancelled with an honest message.
//   (4) TIMINGS -- one entry per completed phase, named, non-negative, sum > 0.
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

#include <cmath>
#include <cstdio>
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
    test_cancellation();
    test_validation_refusal();
    test_concurrency();
    if (g_failures == 0) {
        std::printf("\nOK: the job layer runs what the file says and cannot change the answer\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
