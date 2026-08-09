// Solution verification of the Giroud rigid-footing benchmark (PLAXIS 2D Validation Manual V8,
// Section 2.1): the footing force on three NESTED meshes, and the discretisation uncertainty
// that follows from them.
//
// Why nested. The mesher rebuilds an unstructured mesh from scratch at every density, so three
// independently generated meshes differ by more than h^p and the observed order comes out at
// values no element can deliver (measured: 0.24 / 3.01 / 6.98 on an earlier study). Uniform
// ("red") refinement -- every triangle into four similar children -- gives a ratio of exactly 2
// and preserves every angle, which is what the Richardson/GCI procedure assumes.
//
// The quantity is the footing force at the imposed settlement: F = 2 |sum R_y| over the footing
// nodes (half model, so twice the half-model reaction). The reference values are the manual's:
// 15.15 kN/m analytic (Giroud 1972), 15.24 kN/m PLAXIS 2D.
//
// Not a ctest gate: it is a study, run when the published record needs its number re-measured.
// Build: cmake --build <dir> --target study_giroud_convergence

#include <katai/analysis/results.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/math/grid_convergence.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace gc = katai::math;
namespace m = katai::model;

namespace {

// The manual's own numbers for Section 2.1.
constexpr double kAnalytic = 15.15;   // kN/m, Giroud (1972): F = 2 (1 + nu) G B s / rho
constexpr double kPlaxis = 15.24;     // kN/m, as published in the same section

// F = 2 |sum R_y| over the nodes the prescribed displacement acts on (half model).
double footing_force(const katai::app::SolveResult& r, const m::PrescribedDisp& d) {
    if (r.reaction.size() != (Eigen::Index)(2 * r.mesh.node_count)) return 0.0;
    const double ex = d.x2 - d.x1, ey = d.y2 - d.y1, l2 = ex * ex + ey * ey;
    double ry = 0.0;
    for (int n = 0; n < r.mesh.node_count; ++n) {
        const double t = std::clamp(((r.mesh.x[n] - d.x1) * ex + (r.mesh.y[n] - d.y1) * ey) / l2,
                                    0.0, 1.0);
        if (std::hypot(r.mesh.x[n] - (d.x1 + t * ex), r.mesh.y[n] - (d.y1 + t * ey)) > 1e-6) continue;
        ry += r.reaction[n * 2 + 1];
    }
    return 2.0 * std::fabs(ry);
}

double mesh_area(const katai::mesh::Mesh& mesh) {
    double a = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int i = mesh.node_of(e, 0), j = mesh.node_of(e, 1), k = mesh.node_of(e, 2);
        a += 0.5 * std::fabs((mesh.x[j] - mesh.x[i]) * (mesh.y[k] - mesh.y[i]) -
                             (mesh.x[k] - mesh.x[i]) * (mesh.y[j] - mesh.y[i]));
    }
    return a;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d";
    std::string err;
    m::Project pr;
    if (!m::load_project(path, pr, &err)) {
        std::printf("cannot load %s: %s\n", path.c_str(), err.c_str());
        return 1;
    }
    if (pr.disps.empty()) { std::printf("the case has no prescribed displacement\n"); return 1; }

    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("mesh: %s\n", M.message.c_str()); return 1; }

    std::printf("== Giroud rigid footing, solution verification (nested refinement) ==\n");
    std::printf("   file: %s\n", path.c_str());
    std::printf("   references: analytic %.2f kN/m (Giroud 1972), PLAXIS 2D %.2f kN/m\n\n",
                kAnalytic, kPlaxis);

    // Coarse -> fine: the file's own mesh, then two uniform refinements of it.
    std::vector<katai::mesh::Mesh> meshes;
    meshes.push_back(M.mesh);
    meshes.push_back(katai::mesh::refine_uniform(meshes[0]));
    meshes.push_back(katai::mesh::refine_uniform(meshes[1]));

    std::vector<double> h(3, 0.0), F(3, 0.0);
    for (int i = 0; i < 3; ++i) {
        const auto res = katai::app::solve_phases(pr, meshes[i],
                                                  katai::app::initial_phase_from(pr.initial_procedure),
                                                  nullptr, nullptr, {});
        if (res.empty() || !res.back().ok) {
            std::printf("mesh %d did not solve: %s\n", i,
                        res.empty() ? "no result" : res.back().message.c_str());
            return 1;
        }
        F[i] = footing_force(res.back(), pr.disps[0]);
        h[i] = gc::mesh_representative_size(mesh_area(meshes[i]), meshes[i].element_count);
        std::printf("   mesh %d: %6d elements, %6d nodes, h = %.5f m  ->  F = %.5f kN/m "
                    "(%+.2f%% analytic, %+.2f%% PLAXIS)\n",
                    i, meshes[i].element_count, meshes[i].node_count, h[i], F[i],
                    100.0 * (F[i] - kAnalytic) / kAnalytic,
                    100.0 * (F[i] - kPlaxis) / kPlaxis);
    }

    gc::GridTriplet t;
    t.h1 = h[2]; t.phi1 = F[2];     // fine
    t.h2 = h[1]; t.phi2 = F[1];     // medium
    t.h3 = h[0]; t.phi3 = F[0];     // coarse
    const gc::ConvergenceEstimate e = gc::grid_convergence_band(t);

    std::printf("\n   refinement ratios r21 = %.4f, r32 = %.4f\n", e.r21, e.r32);
    std::printf("   observed order p = %.4f (%s)\n", e.p, gc::convergence_kind_name(e.kind));
    if (!e.message.empty()) std::printf("   note: %s\n", e.message.c_str());
    std::printf("   Richardson extrapolation F(h->0) = %.5f kN/m (%+.2f%% analytic, %+.2f%% PLAXIS)\n",
                e.phi_extrapolated, 100.0 * (e.phi_extrapolated - kAnalytic) / kAnalytic,
                100.0 * (e.phi_extrapolated - kPlaxis) / kPlaxis);
    std::printf("   asymptotic ratio %.4f, extrapolation quotable: %s\n", e.asymptotic_ratio,
                e.asymptotic ? "yes" : "no");
    std::printf("   NUMERICAL UNCERTAINTY of the fine mesh: +/- %.4f%% (%s)\n",
                100.0 * e.band, e.band_basis.c_str());
    std::printf("\n   => F(fine) = %.4f kN/m +/- %.4f%%; distance to PLAXIS %+.2f%%, "
                "to the analytic value %+.2f%%\n",
                F[2], 100.0 * e.band, 100.0 * (F[2] - kPlaxis) / kPlaxis,
                100.0 * (F[2] - kAnalytic) / kAnalytic);
    return 0;
}
