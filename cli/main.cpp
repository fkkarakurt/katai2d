// katai -- the command-line front end (layer 4, a composition root). It programs
// against the published facade exclusively, and every solve goes through
// katai::api::Job -- the SAME layer-2 work the GUI's solve button submits, so the
// front end cannot change the answer. tests/test_cli.cpp pins that across processes:
// the .res this program writes is bit-identical to an in-process Job run.
//
// Subcommands:
//   katai solve <file.k2d> [--out <file.res>]   run every phase; progress on stderr
//   katai validate <file.k2d>                    print the input-contract report
//   katai info                                   the file-format versions this build speaks
//
// Exit codes (documented here, pinned by test_cli):
//   0  success (validate: no errors; solve: every phase converged)
//   2  usage error
//   3  a file cannot be read, parsed or written
//   4  the input contract refused the project: validation errors, or reader ERROR
//      notes (a forward-version file whose enums this build had to clamp) -- the
//      reader-notes rule: display may proceed, a solve must not
//   5  the solve failed (a phase did not converge or the engine refused it)

#include <katai/api/katai.hpp>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace api = katai::api;

namespace {

int usage() {
    std::fprintf(stderr,
                 "katai -- KATAI 2D command line\n"
                 "usage:\n"
                 "  katai solve <file.k2d> [--out <file.res>]\n"
                 "  katai validate <file.k2d>\n"
                 "  katai info\n"
                 "exit codes: 0 ok, 2 usage, 3 file error, 4 input contract refused, "
                 "5 solve failed\n");
    return 2;
}

void print_issues(const std::vector<api::Issue>& issues) {
    for (const api::Issue& i : issues)
        std::printf("%s: %s: %s\n",
                    i.severity == api::Severity::Error ? "error" : "warning",
                    i.path.c_str(), i.message.c_str());
}

// What the run did that the file does not literally say: a line clipped to the soil, a
// fallback taken, an object refused. The input contract is checked before a mesh exists,
// so these are the only report a user gets on everything that is decided against the mesh
// -- printing them is what keeps a discarded load from passing as a clean run.
void print_diagnostics(const api::SolveResult& R) {
    for (const api::Diagnostic& d : R.diagnostics) {
        const char* sev = d.severity == api::DiagnosticSeverity::Refusal  ? "refusal"
                        : d.severity == api::DiagnosticSeverity::Warning  ? "warning"
                                                                          : "note   ";
        std::printf("  %s %s  %s%s%s\n", sev, d.code.c_str(),
                    d.subject.empty() ? "" : d.subject.c_str(), d.subject.empty() ? "" : ": ",
                    d.message.c_str());
    }
}

bool any_error_note(const std::vector<api::Issue>& notes) {
    for (const api::Issue& i : notes)
        if (i.severity == api::Severity::Error) return true;
    return false;
}

// Load + surface the reader's forward-version notes. Returns 0, or the exit code.
int load(const std::string& path, api::Project& pr, std::vector<api::Issue>& notes) {
    std::string err;
    if (!api::load_project(path, pr, &err, &notes)) {
        std::fprintf(stderr, "error: cannot read %s: %s\n", path.c_str(), err.c_str());
        return 3;
    }
    print_issues(notes);
    return 0;
}

int cmd_validate(const std::string& path) {
    api::Project pr;
    std::vector<api::Issue> notes;
    if (const int rc = load(path, pr, notes)) return rc;
    const api::ValidationReport rep = api::validate_project(pr);
    std::fputs(rep.to_string().c_str(), stdout);
    if (!rep.ok() || any_error_note(notes)) {
        std::printf("REFUSED: the project fails the input contract (%d error(s), %d warning(s))\n",
                    rep.errors, rep.warnings);
        return 4;
    }
    std::printf("OK: %s satisfies the input contract (%d warning(s))\n", path.c_str(),
                rep.warnings);
    return 0;
}

int cmd_solve(const std::string& path, const std::string& out) {
    api::Project pr;
    std::vector<api::Issue> notes;
    if (const int rc = load(path, pr, notes)) return rc;
    if (any_error_note(notes)) {
        std::printf("REFUSED: the file carries reader ERROR notes (see above); "
                    "a solve must not proceed on clamped input\n");
        return 4;
    }

    api::Job job(pr);
    job.set_on_phase([](int cur, int total, const std::string& name) {
        std::fprintf(stderr, "phase %d/%d: %s\n", cur + 1, total, name.c_str());
    });
    if (!job.run()) {
        if (!job.report().ok()) {
            std::fputs(job.report().to_string().c_str(), stdout);
            std::printf("REFUSED: %s\n", job.message().c_str());
            return 4;
        }
        // Every phase that ran, including the one that stopped the job: its diagnostics name
        // the object at fault, so they must survive the failure path.
        const std::vector<api::SolveResult>& partial = job.results();
        for (size_t i = 0; i < partial.size(); ++i) {
            if (partial[i].diagnostics.empty()) continue;
            std::printf("phase %zu:\n", i + 1);
            print_diagnostics(partial[i]);
        }
        std::printf("FAILED: %s\n", job.message().c_str());
        return 5;
    }

    const std::vector<api::SolveResult>& res = job.results();
    double total_s = 0.0;
    for (const api::PhaseTiming& t : job.timings()) total_s += t.seconds;
    for (size_t i = 0; i < res.size(); ++i) {
        const api::SolveResult& R = res[i];
        std::printf("phase %zu/%zu: ok  max|u| = %.6e m", i + 1, res.size(), R.max_disp);
        if (R.fos >= 0.0)
            std::printf("  FoS %s %.3f", R.fos_lower_bound ? ">" : "=", R.fos);
        std::printf("\n");
        print_diagnostics(R);
    }
    std::printf("solved %zu phase(s) in %.2f s\n", res.size(), total_s);

    if (!out.empty()) {
        std::string err;
        const std::uint64_t hash = api::fnv1a64(api::project_to_json(pr));
        if (!api::save_results(out, hash, res, &err)) {
            std::fprintf(stderr, "error: cannot write %s: %s\n", out.c_str(), err.c_str());
            return 3;
        }
        std::printf("results written to %s\n", out.c_str());
    }
    return 0;
}

int cmd_info() {
    std::printf("%s %s (%s)\n", api::kAppName, api::kVersion, api::kVersionDate);
    std::printf("  .k2d project file version: %d\n", api::kProjectFileVersion);
    std::printf("  .res results file version: %u\n", (unsigned)api::kResultsFileVersion);
    std::printf("  linear solver backend: %s\n", api::backend_name());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return usage();
    const std::string cmd = argv[1];
    if (cmd == "info") return argc == 2 ? cmd_info() : usage();
    if (cmd == "validate") return argc == 3 ? cmd_validate(argv[2]) : usage();
    if (cmd == "solve") {
        std::string file, out;
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) out = argv[++i];
            else if (argv[i][0] == '-' || !file.empty()) return usage();
            else file = argv[i];
        }
        return file.empty() ? usage() : cmd_solve(file, out);
    }
    return usage();
}
