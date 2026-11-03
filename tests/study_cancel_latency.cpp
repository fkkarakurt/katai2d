// Measurement instrument for cancelling a run (a `study_` target: EXCLUDE_FROM_ALL, not wired into
// CTest -- its numbers are wall-clock times). A run is stopped at its next poll of the request
// (katai/analysis/cancel.hpp), so the longest a user waits after pressing Cancel is the longest
// stretch of the run between two polls. This runs a project once with a request that is never set
// and records the time of every poll: the largest gaps, and the phase each falls in, are the
// latency a Cancel can meet there. The meshing before the first phase is timed apart. With a time
// T the request is set at T seconds into the run, and the time the run then takes to return is
// the latency itself, the unwinding included.
//
// usage: study_cancel_latency <project.k2d> [gaps to list] [T]
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace m = katai::model;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: study_cancel_latency <project.k2d> [gaps to list]\n");
        return 2;
    }
    const int list = argc > 2 ? std::atoi(argv[2]) : 8;
    const double cancel_at = argc > 3 ? std::atof(argv[3]) : -1.0;
    m::Project pr;
    std::string err;
    if (!m::load_project(argv[1], pr, &err)) {
        std::printf("cannot read %s: %s\n", argv[1], err.c_str());
        return 1;
    }
    const auto t0 = clk::now();
    const auto M = katai::app::mesh_from_project(pr);
    const double t_mesh = std::chrono::duration<double>(clk::now() - t0).count();
    if (!M.ok) {
        std::printf("meshing failed: %s\n", M.message.c_str());
        return 1;
    }
    struct Poll { double t; int phase; };
    std::vector<Poll> polls;
    std::vector<std::string> names;
    int phase = -1;
    const auto t1 = clk::now();
    auto now = [&] { return std::chrono::duration<double>(clk::now() - t1).count(); };
    const auto progress = [&](int cur, int, const std::string& name) {
        phase = cur;
        names.resize((size_t)std::max<int>((int)names.size(), cur + 1));
        names[(size_t)cur] = name;
        polls.push_back({now(), cur});   // a phase boundary is a poll too
    };
    const auto request = [&] {
        const double t = now();
        polls.push_back({t, phase});
        return cancel_at >= 0.0 && t >= cancel_at;
    };
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure),
                                              progress, request);
    const double t_end = now();
    polls.push_back({t_end, phase});
    std::printf("%s: %d nodes, meshing %.2f s (no poll), run %.2f s, %zu phase(s), %zu polls\n", argv[1],
                M.mesh.node_count, t_mesh, t_end, res.size(), polls.size());
    if (cancel_at >= 0.0 && t_end > cancel_at)
        std::printf("  cancel requested at %.2f s: the run returned %.3f s later, in phase %d\n", cancel_at,
                    t_end - cancel_at, phase);
    struct Gap { double len, at; int phase; };
    std::vector<Gap> gaps;
    for (size_t k = 1; k < polls.size(); ++k)
        gaps.push_back({polls[k].t - polls[k - 1].t, polls[k - 1].t, polls[k].phase});
    std::sort(gaps.begin(), gaps.end(), [](const Gap& a, const Gap& b) { return a.len > b.len; });
    for (int k = 0; k < list && k < (int)gaps.size(); ++k) {
        const Gap& g = gaps[(size_t)k];
        std::printf("  gap %8.3f s  from t = %8.2f s  in phase %d (%s)\n", g.len, g.at, g.phase,
                    g.phase >= 0 && g.phase < (int)names.size() ? names[(size_t)g.phase].c_str() : "-");
    }
    for (size_t k = 0; k < res.size(); ++k)
        std::printf("  phase %zu: %s\n", k, res[k].ok ? "ok" : res[k].message.substr(0, 100).c_str());
    return 0;
}
