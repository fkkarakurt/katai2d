// Points for curves (Project::curves, .k2d v30): the node nearest to each point has its displacement
// recorded at every step of every phase -- against Sum-Mstage in a Plastic phase, Sum-Msf in a
// Safety phase, time in a consolidation phase.
//   (a) a linear elastic block loaded in five steps by a strip load and, beside it, a strip driven
//       down by a prescribed displacement: every point's displacement is Sum-Mstage times its
//       displacement at the end (linear elasticity), the curve ends at the node's displacement in
//       the phase's result, and the point on the driven strip moves with the strip (Sum-Mstage
//       times the prescribed 10 mm) -- read from the free unknowns alone it stood still;
//   (b) the same for the phase's path: its largest displacement is the largest of the result,
//       the driven strip's included;
//   (c) the Terzaghi column of the corpus: the curve of the point at its top is the settlement-time
//       curve of the column (the largest settlement of the model, which its top's nodes share to
//       within the discretisation), and ends at its node's displacement;
//   (d) the Griffiths-Lane slope, a Safety phase: the point at the crest starts at no displacement
//       and its curve runs beside the path, never above the largest displacement;
//   (e) the points survive the project file and the curves the results file; a point outside the
//       model is warned about.
//
// verify: KV-NUM-017
//   oracle:   closed_form
//   source:   linearity of the elastic boundary value problem: the displacement under a load (and a prescribed displacement) scaled by a factor is the displacement under the whole load scaled by it (KATAI 2D, docs/k2d-format.md, `curves`)
//   locator:  a 10 x 5 m linear elastic block (E = 10 MPa, nu = 0.3), base fixed, sides on rollers, tri6 at 0.5 m; a 2 m strip load of 100 kPa at x = 0..2 and a 2 m strip driven down 10 mm at x = 6..8, applied in one Plastic phase of five steps; points at (1, 5), (7, 5) and (4, 2.5)
//   quantity: the displacement of each point at each step [m]
//   expected: u_k = Sum-Mstage_k * u_end; the driven point's u_y = -0.01 Sum-Mstage_k
//   band:     1e-9 relative (an elastic solve is exact to rounding); the column's curve within 1e-3 of its settlement curve and never above it
#include <katai/io/project_io.hpp>
#include <katai/io/results_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
// The displacement of point i at step k (component c).
double cu(const katai::app::SolveResult& R, size_t k, size_t i, int c) {
    return R.curve_u[(k * R.curve_nodes.size() + i) * 2 + (size_t)c];
}
m::CurvePoint point(const char* name, double x, double y) {
    m::CurvePoint C;
    C.name = name; C.x = x; C.y = y;
    return C;
}
}  // namespace

int main() {
    std::printf("== Points for curves (KV-NUM-017) ==\n");
    std::string err;
    {
        // (a), (b) the elastic block
        m::Project pr;
        pr.x_min = 0; pr.x_max = 10; pr.y_min = 0; pr.y_max = 5;
        pr.has_water = false;
        pr.mesh.elem_size = 0.5;
        m::Material soil;
        soil.name = "Elastic"; soil.model = m::SoilModel::LinearElastic;
        soil.E = 1.0e4; soil.nu = 0.3; soil.gamma_unsat = soil.gamma_sat = 0.0;
        pr.materials.push_back(soil);
        m::SoilPolygon P;
        P.name = "Block"; P.material = 0;
        P.x = {0, 10, 10, 0}; P.y = {0, 0, 5, 5};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                     (int)m::BCType::HorizontallyFixed};
        pr.polygons.push_back(P);
        m::Load q;
        q.kind = m::LoadKind::Distributed; q.name = "Strip";
        q.x1 = 0; q.y1 = 5; q.x2 = 2; q.y2 = 5; q.qy1 = q.qy2 = -100.0;
        pr.loads.push_back(q);
        m::PrescribedDisp D;
        D.name = "Driven"; D.x1 = 6; D.y1 = 5; D.x2 = 8; D.y2 = 5; D.set_uy = true; D.uy = -0.01;
        pr.disps.push_back(D);
        pr.curves = {point("A", 1, 5), point("B", 7, 5), point("C", 4, 2.5)};
        pr.initial.load_active = {0};
        pr.initial.disp_active = {0};
        m::Phase ph;
        ph.name = "Load"; ph.type = m::PhaseType::Plastic;
        ph.load_active = {1}; ph.disp_active = {1}; ph.load_steps = 5;
        pr.phases = {ph};
        const auto M = katai::app::mesh_from_project(pr);
        const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
        if (res.size() != 2 || !res[1].ok) {
            std::printf("FATAL: %s\n", res.empty() ? "no result" : res.back().message.c_str());
            return 1;
        }
        const auto& R = res[1];
        const size_t n = R.path_factor.size(), np = R.curve_nodes.size();
        std::printf("   %zu steps, %zu points; the driven point ends at u_y = %.9f m\n", n, np,
                    np == 3 && R.curve_u.size() == n * 6 ? cu(R, n - 1, 1, 1) : 0.0);
        bool shaped = np == 3 && n >= 4 && R.curve_u.size() == n * 2 * np;
        check(shaped, "(a) three points, one displacement pair each at every step of the path");
        if (shaped) {
            double worst = 0.0, worst_end = 0.0, worst_driven = 0.0;
            for (size_t i = 0; i < np; ++i)
                for (int c = 0; c < 2; ++c) {
                    const double end = cu(R, n - 1, i, c);
                    const double node = R.disp(2 * R.curve_nodes[i] + c);
                    worst_end = std::max(worst_end, std::fabs(end - node) / std::max(std::fabs(node), 1e-12));
                    for (size_t k = 0; k < n; ++k)
                        worst = std::max(worst, std::fabs(cu(R, k, i, c) - R.path_factor[k] * end) /
                                                    std::max(std::fabs(end), 1e-12));
                }
            for (size_t k = 0; k < n; ++k)
                worst_driven = std::max(worst_driven, std::fabs(cu(R, k, 1, 1) + 0.01 * R.path_factor[k]) / 0.01);
            std::printf("   worst |u_k - Mstage_k u_end| / |u_end| = %.2e; end against the result %.2e; "
                        "driven point against Mstage x 10 mm %.2e\n", worst, worst_end, worst_driven);
            check(worst < 1e-9, "(a) every point's displacement is Sum-Mstage times its displacement at the end");
            check(worst_end < 1e-12, "(a) each curve ends at its node's displacement in the phase's result");
            check(worst_driven < 1e-12, "(a) the point on the driven strip moves with the strip");
            double umax = 0.0;
            for (int k = 0; k < (int)(R.disp.size() / 2); ++k)
                umax = std::max(umax, std::hypot(R.disp(2 * k), R.disp(2 * k + 1)));
            check(std::fabs(R.path_disp.back() - umax) <= 1e-12 * umax && umax >= 0.01,
                  "(b) the path's largest displacement is the result's, the driven strip's included");
        }

        // (e) the project file and the results file
        m::Project back;
        check(m::project_from_json(m::project_to_json(pr), back, &err) && back.curves.size() == 3 &&
                  back.curves[1].name == "B" && back.curves[1].x == 7.0 && back.curves[1].y == 5.0,
              "(e) the points come back from the project file");
        const std::string file = "test_curve_points.res";
        std::vector<katai::app::SolveResult> rback;
        check(katai::app::save_results(file, 0x5151, res, &err) && katai::app::load_results(file, 0x5151, rback, &err) &&
                  rback.size() == 2 && rback[1].curve_nodes == R.curve_nodes && rback[1].curve_u == R.curve_u,
              "(e) the curves come back from the results file exactly");
        std::remove(file.c_str());
        m::Project outside = pr;
        outside.curves.push_back(point("Far", 30, 5));
        const auto rep = katai::io::validate_project(outside);
        bool warned = false;
        for (const auto& is : rep.issues)
            warned = warned || (is.severity == katai::io::Severity::Warning && is.path.rfind("curves[3]", 0) == 0);
        check(warned, "(e) a point outside the model is warned about");
    }
    {
        // (c) the Terzaghi column
        m::Project pr;
        if (!m::load_project(std::string(KATAI_CORPUS_DIR) + "/kv-con-002-terzaghi-column.k2d", pr, &err)) {
            std::printf("FATAL: %s\n", err.c_str());
            return 1;
        }
        double ytop = pr.polygons[0].y[0], xmid = 0.0;
        for (size_t k = 0; k < pr.polygons[0].y.size(); ++k) {
            ytop = std::max(ytop, pr.polygons[0].y[k]);
            xmid += pr.polygons[0].x[k] / (double)pr.polygons[0].x.size();
        }
        pr.curves = {point("Top", xmid, ytop)};
        const auto M = katai::app::mesh_from_project(pr);
        const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
        const katai::app::SolveResult* C = nullptr;
        for (const auto& r : res)
            if (r.ok && r.consol_time.size() > 2) C = &r;
        // The settlement curve is the largest settlement of the model at each step; the top's nodes
        // settle alike to within the discretisation (a 2-D mesh of a 1-D problem), so the top
        // point's curve lies on it, never above it.
        bool same = C && C->curve_nodes.size() == 1 && C->curve_u.size() == 2 * C->consol_time.size();
        double worst = 0.0, above = 0.0, end_gap = 1.0;
        if (same) {
            for (size_t k = 1; k < C->consol_time.size(); ++k) {
                const double s = std::fabs(cu(*C, k, 0, 1)), smax = C->consol_settlement[k];
                worst = std::max(worst, (smax - s) / smax);
                above = std::max(above, s - smax);
            }
            const size_t last = C->consol_time.size() - 1;
            end_gap = std::fabs(cu(*C, last, 0, 1) - C->disp(2 * C->curve_nodes[0] + 1));
        }
        std::printf("   column: %zu time steps, the top's curve below the settlement curve by at most %.2e of it; "
                    "its end against the node's result %.2e m\n",
                    C ? C->consol_time.size() : (size_t)0, worst, end_gap);
        check(same && worst <= 1e-3 && above <= 0.0,
              "(c) the column's top point follows its settlement-time curve, never above it");
        check(same && end_gap <= 1e-15, "(c) ...and ends at its node's displacement in the phase's result");
        check(same && cu(*C, 0, 0, 1) == 0.0 && cu(*C, C->consol_time.size() - 1, 0, 1) < 0.0,
              "(c) ...starting at no displacement and settling");
    }
    {
        // (d) the slope, a Safety phase
        m::Project slope;
        if (!m::load_project(std::string(KATAI_CORPUS_DIR) + "/kv-slp-001-griffiths-lane-slope.k2d", slope, &err)) {
            std::printf("FATAL: %s\n", err.c_str());
            return 1;
        }
        slope.initial_procedure = m::InitialProcedure::GravityLoading;
        m::Phase safety;
        safety.name = "Safety";
        safety.type = m::PhaseType::Safety;
        slope.phases.push_back(safety);
        slope.curves = {point("Crest", 50, 35)};
        const auto MS = katai::app::mesh_from_project(slope);
        const auto res = katai::app::solve_phases(slope, MS.mesh, katai::app::initial_phase_from(slope.initial_procedure));
        const bool ran = res.size() == 2 && res[1].ok;
        const auto& S = res.back();
        const size_t n = S.path_factor.size();
        bool beside = ran && S.curve_nodes.size() == 1 && S.curve_u.size() == 2 * n && n > 2;
        double last = 0.0;
        if (beside) {
            beside = cu(S, 0, 0, 0) == 0.0 && cu(S, 0, 0, 1) == 0.0;
            for (size_t k = 0; k < n; ++k)
                beside = beside && std::hypot(cu(S, k, 0, 0), cu(S, k, 0, 1)) <= S.path_disp[k] * (1.0 + 1e-12);
            last = std::hypot(cu(S, n - 1, 0, 0), cu(S, n - 1, 0, 1));
        }
        std::printf("   slope: FoS %.4f, %zu steps; the crest moves %.4e m, the largest %.4e m\n", S.fos, n, last,
                    n ? S.path_disp.back() : 0.0);
        check(beside && last > 0.0, "(d) the crest's curve starts at rest and never exceeds the largest displacement");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
