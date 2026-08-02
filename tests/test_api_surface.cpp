// The published facade pinned for SUFFICIENCY: this translation unit includes exactly
// ONE project header -- <katai/api/katai.hpp> -- and standard headers, and still
// performs the whole front-end workflow: construct a project naming every kind of
// vocabulary (component types, enums, mesh settings, initial procedure), save it,
// load it back with reader notes, validate it, run it as a Job with progress, read
// the results against a closed form, and watch an invalid project be REFUSED by the
// job itself with its field-path report. If a name this file needs is missing from
// the facade, the facade is insufficient and this file stops compiling -- which is
// the point (Stage D, section 7.3: total over the input, closed over the
// implementation).
//
// verify: none -- api-surface sufficiency and refusal mechanics; the physics number
// checked here (drained confined column u = -qH/M') is a smoke sanity, and the real
// oracle record lives in the corpus (test_input_corpus) and the engine suites.
#include <katai/api/katai.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace api = katai::api;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// A drained, weightless, laterally confined column under a top surcharge: the
// simplest project that still names materials, polygons, BC enums, loads, mesh
// settings and the initial procedure through the facade.
constexpr double kE = 1.0e4, kNu = 0.3, kH = 10.0, kW = 2.0, kQ = 50.0;

api::Project make_column() {
    api::Project pr;
    pr.name = "api-surface column";
    pr.x_min = 0.0; pr.x_max = kW;
    pr.y_min = 0.0; pr.y_max = kH;
    pr.has_water = false;
    pr.initial_procedure = api::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    api::Material s;
    s.name = "Weightless elastic";
    s.model = api::SoilModel::LinearElastic;
    s.drainage = api::Drainage::Drained;
    s.E = kE; s.nu = kNu;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);

    api::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)api::BCType::FullyFixed, (int)api::BCType::NormallyFixed,
                 (int)api::BCType::Free, (int)api::BCType::NormallyFixed};
    pr.polygons.push_back(P);

    api::Load L;
    L.kind = api::LoadKind::Distributed;
    L.name = "Surcharge";
    L.x1 = 0; L.y1 = kH; L.x2 = kW; L.y2 = kH;
    L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    return pr;
}

void test_workflow() {
    std::printf("-- construct -> save -> load -> validate -> run -> read --\n");
    const api::Project built = make_column();

    // Files, through the facade (written into the ctest working directory).
    const std::string path = "api_surface_roundtrip.k2d";
    std::string err;
    check(api::save_project(built, path, &err), "save_project through the facade");
    api::Project pr;
    std::vector<api::Issue> notes;
    check(api::load_project(path, pr, &err, &notes), "load_project through the facade");
    check(notes.empty(), "reader reports no notes");
    check(api::project_to_json(pr) == api::project_to_json(built),
          "round trip is JSON-identical");
    std::remove(path.c_str());

    // Validation, through the facade.
    const api::ValidationReport rep = api::validate_project(pr);
    check(rep.ok() && rep.issues.empty(), "the column validates clean");

    // The run, through the facade.
    api::Job job(pr);
    int announced = 0;
    job.set_on_phase([&announced](int, int, const std::string&) { ++announced; });
    check(job.run(), "the job ran");
    check(job.state() == api::JobState::Done, "job state Done");
    check(announced == 1, "progress announced the single phase");
    check(job.timings().size() == 1, "one timing recorded");

    // Results, through the facade: u_top = -q H / M', M' = E(1-nu)/((1+nu)(1-2nu)).
    const double Mp = kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));
    const double u_exact = -kQ * kH / Mp;
    const api::SolveResult& R = job.results().back();
    double u_top = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > kH - 1e-6) u_top = std::fmin(u_top, R.disp[n * 2 + 1]);
    std::printf("   u_top = %.5e m (closed form %.5e, %+.2f%%)\n", u_top, u_exact,
                100.0 * (u_top - u_exact) / u_exact);
    check(std::fabs(u_top - u_exact) < 0.02 * std::fabs(u_exact),
          "drained confined settlement -qH/M' within 2% (smoke sanity)");
}

void test_refusal() {
    std::printf("-- an invalid project is refused by the job itself --\n");
    api::Project bad = make_column();
    bad.materials[0].E = -1.0;   // indefensible on sight
    api::Job job(bad);
    check(!job.run(), "run() reports failure");
    check(job.state() == api::JobState::Failed, "state is Failed");
    check(!job.report().ok() && job.report().errors >= 1, "the validation report carries the error");
    bool at_path = false;
    for (const api::Issue& i : job.report().issues)
        if (i.severity == api::Severity::Error && i.path == "materials[0].E") at_path = true;
    check(at_path, "the error is reported at its exact field path");
    check(job.results().empty() && job.mesh().mesh.node_count == 0,
          "nothing was meshed or solved");
    std::printf("   message: %s\n", job.message().c_str());
    check(job.message().find("materials[0].E") != std::string::npos,
          "the job's message names the field path");
}

}  // namespace

int main() {
    std::printf("The published facade: one header, the whole front-end workflow\n\n");
    test_workflow();
    test_refusal();
    if (g_failures == 0) {
        std::printf("\nOK: the api surface suffices, and it cannot solve an invalid project\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
