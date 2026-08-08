// The first published number to carry a numerical-uncertainty band, and the machinery that
// gives it one.
//
// Every benchmark in this repository has so far been reported at ONE mesh density -- "+1.2% on
// the file's own 1 m tri6 mesh". A reader cannot tell from that how much of the deviation is
// the MODEL and how much is the DISCRETISATION, and neither can the author. This case answers
// that for the strip-load benchmark by solving the same file at three systematically refined
// densities and putting the triplet through the Grid Convergence Index (katai/math/
// grid_convergence.hpp, verified against manufactured data in KV-NUM-004).
//
// The strip load is the right case to do first because its answer is known in closed form
// (Boussinesq integrated over a uniform strip), so the study can be checked in a way that a
// case without an analytic solution cannot: the Richardson-extrapolated value -- the estimate
// of what this formulation would give on an infinitely fine mesh -- must sit CLOSER to the
// closed form than the finest mesh's own value does. That is the whole claim of the method,
// and here it is measured rather than assumed.
//
// The three meshes are NESTED: the mesher is visited once and the finer two are made by
// splitting every triangle at its edge midpoints (katai::mesh::refine_uniform, KV-NUM-006).
// That is not a convenience -- it is what makes the study mean anything. Three separate visits
// to the mesher give three unrelated meshes whose differences carry the mesher's irregularity
// as well as the discretisation, and measured on this tree that contamination was enough to
// put the same study outside the asymptotic range. The representative cell size h is measured
// from each mesh as sqrt(total area / element count), not taken from the requested element
// size, because the achieved density is what the solution actually saw (Celik et al., step 1).
//
// verify: KV-NUM-005
//   oracle:   closed_form
//   source:   classical elasticity: Boussinesq point solution integrated over a uniform strip on a half-plane; discretisation-error estimation per P. J. Roache (1994), ASME J. Fluids Eng. 116(3):405-413 and I. B. Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001 (Grid Convergence Index, factor of safety 1.25 for a three-mesh study), the solution-verification framework of ASME V&V 20-2009 / V&V 10-2006
//   locator:  sigma_z = (q/pi)[alpha + sin(alpha) cos(theta1 + theta2)], theta_i = atan((x -+ a)/z) from the strip edges, alpha = theta1 - theta2; h = sqrt(total area / element count); p = |ln|eps32/eps21| + q(p)|/ln(r21); phi_ext = (r21^p phi1 - phi2)/(r21^p - 1); GCI_fine = 1.25 e_a/(r21^p - 1) (stated in full)
//   quantity: vertical stress sigma_z at 2 m and 4 m depth on the strip centre line and the settlement under the strip centre, interpolated with the element shape functions at the exact probe point, from tests/corpus/kv-fnd-008-strip-load.k2d solved on a NESTED family (mesher visited once at elem_size 1.4 m, then refined uniformly twice: 1190 / 4760 / 19040 elements, h = 0.8199 / 0.4100 / 0.2050 m, refinement factor exactly 2) [kPa; m]
//   expected: the closed form above at each stress probe; a reported band that covers the distance from the fine mesh to the limit it estimates; and, where the triplet is in the asymptotic range, a Richardson-extrapolated value closer to the closed form than the finest mesh's own
//   band:     as asserted below and MEASURED on this tree, not inherited -- sigma_z(2 m): p = 3.65 monotonic, band +/-0.003% from the observed order, fine mesh -0.067% from the closed form (so the deviation is NOT discretisation), extrapolated -0.06%; sigma_z(4 m): oscillatory, band +/-1.066% from the assumed order p = 2 with Fs = 3 widened to the observed spread, fine mesh -0.398%; settlement: p = 3.90 monotonic, band +/-0.001%. The generic assertions are: refinement factor exactly 2, band covers the distance to the limit, band below 5%, fine mesh within 4% of the closed form

#include <katai/fem/elements/point_location.hpp>
#include <katai/fem/elements/tri15.hpp>
#include <katai/fem/elements/tri6.hpp>
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

namespace m = katai::model;
namespace gc = katai::math;

namespace {

constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;

void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

// Boussinesq strip load on a weightless elastic half-plane (independent of E and nu): the same
// closed form KV-FND-008 already verifies against, evaluated on the strip centre line.
// q is the applied pressure (positive), a the half width, z the depth below the surface.
double boussinesq_strip_sigma_z(double q, double a, double z) {
    const double th1 = std::atan(a / z), th2 = std::atan(-a / z);
    const double alpha = th1 - th2;
    return (q / kPi) * (alpha + std::sin(alpha) * std::cos(th1 + th2));
}

// A nodal field read at an arbitrary physical point, with the element's own shape functions --
// NOT at the nearest node. On a convergence study the nearest node moves with every mesh, and
// that sampling jitter would be indistinguishable from discretisation error; interpolating at
// the exact probe point removes it.
template <class E>
bool interpolate_at(const katai::mesh::Mesh& mesh, const std::vector<double>& nodal, double px,
                    double py, double& out) {
    const auto loc = katai::core::ploc::locate_point(mesh, px, py);
    if (!loc.found) return false;
    typename E::NodeCoords X;
    for (int k = 0; k < E::kNodeCount; ++k) {
        X(k, 0) = mesh.x[mesh.node_of(loc.element, k)];
        X(k, 1) = mesh.y[mesh.node_of(loc.element, k)];
    }
    const auto lc = katai::core::ploc::physical_to_local<E>(X, px, py);
    if (!lc.converged) return false;
    const auto N = E::shape_functions(lc.xi, lc.eta);
    out = 0.0;
    for (int k = 0; k < E::kNodeCount; ++k) out += N(k) * nodal[mesh.node_of(loc.element, k)];
    return true;
}

struct MeshRun {
    bool ok = false;
    double h = 0.0;             // representative cell size, measured from the mesh
    int elements = 0, nodes = 0;
    std::vector<double> sigma_z;  // one per probe depth [kPa], sign as reported (compression < 0)
    double u_y = 0.0;             // settlement under the strip centre [m]
};

// Solve the case file ON A GIVEN MESH and read the probes. The mesh is passed in rather than
// derived from the project, because the three meshes of this study are related to each other by
// refinement -- not by three separate visits to the mesher.
MeshRun run_on(const m::Project& pr, const katai::mesh::Mesh& mesh,
               const std::vector<double>& probe_depth, double y_surface) {
    MeshRun r;
    const auto& M = mesh;

    double area = 0.0;
    for (int e = 0; e < M.element_count; ++e) {
        const int a = M.node_of(e, 0), b = M.node_of(e, 1), c = M.node_of(e, 2);
        area += 0.5 * std::fabs((M.x[b] - M.x[a]) * (M.y[c] - M.y[a]) -
                                (M.x[c] - M.x[a]) * (M.y[b] - M.y[a]));
    }
    r.h = gc::mesh_representative_size(area, M.element_count);
    r.elements = M.element_count;
    r.nodes = M.node_count;

    const auto res = katai::app::solve_phases(pr, M,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res.back().ok) {
        std::printf("      (solve: %s)\n", res.empty() ? "no phases" : res.back().message.c_str());
        return r;
    }
    // Nodal sigma_yy from the recovered effective-stress field.
    std::vector<double> syy(M.node_count, 0.0);
    for (int n = 0; n < M.node_count && n < (int)res.back().stress.stress.size(); ++n)
        syy[n] = res.back().stress.stress[n](1);

    const double x_centre = 20.0;   // the strip runs x = 18..22 in the case file
    for (double d : probe_depth) {
        double v = 0.0;
        const bool got = M.nodes_per_element == 15
                             ? interpolate_at<katai::core::Tri15Element>(M, syy, x_centre,
                                                                        y_surface - d, v)
                             : interpolate_at<katai::core::Tri6Element>(M, syy, x_centre,
                                                                       y_surface - d, v);
        if (!got) return r;
        r.sigma_z.push_back(v);
    }
    // The settlement under the strip centre, from the primary unknown. Recovered stress is a
    // DERIVED field -- one differentiation and a nodal averaging away from the solution -- so it
    // converges more slowly and more raggedly than the displacement it came from. Carrying both
    // is the point: they answer different questions about the same three meshes.
    std::vector<double> uy(M.node_count, 0.0);
    for (int n = 0; n < M.node_count && 2 * n + 1 < (int)res.back().disp.size(); ++n)
        uy[n] = res.back().disp[2 * n + 1];
    const bool got_u = M.nodes_per_element == 15
                           ? interpolate_at<katai::core::Tri15Element>(M, uy, x_centre,
                                                                      y_surface, r.u_y)
                           : interpolate_at<katai::core::Tri6Element>(M, uy, x_centre,
                                                                     y_surface, r.u_y);
    if (!got_u) return r;
    r.ok = true;
    return r;
}

}  // namespace

int main() {
    std::printf("== mesh convergence of the strip-load benchmark (KV-NUM-005) ==\n");

    m::Project base;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-fnd-008-strip-load.k2d";
    if (!m::load_project(path, base, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        return 1;
    }
    const double q = 100.0, a = 2.0, y_surface = base.y_max;   // the file's own strip
    const std::vector<double> probes = {2.0, 4.0};

    // The three meshes are NESTED: the mesher is visited once, and the finer two are made by
    // splitting every triangle at its edge midpoints (katai::mesh::refine_uniform, verified in
    // KV-NUM-006). That fixes the refinement factor at exactly 2, preserves every angle, and
    // keeps every node of the coarse mesh -- so the differences between these three solutions
    // are the element SIZE and nothing else. Three visits to the mesher instead would give
    // three unrelated meshes, and their differences would carry the mesher's irregularity as
    // well; measured on this tree, that irregularity was large enough to put the same study
    // outside the asymptotic range (docs/validation/numerical-uncertainty.md).
    m::Project pr = base;
    pr.mesh.elem_size = 1.4;   // the coarse level; the finest is this refined twice
    const auto M0 = katai::app::mesh_from_project(pr);
    check(M0.ok, "base mesh built from the case file");
    if (!M0.ok) { std::printf("      (%s)\n", M0.message.c_str()); return 1; }
    const katai::mesh::Mesh m1 = katai::mesh::refine_uniform(M0.mesh);
    const katai::mesh::Mesh m2 = katai::mesh::refine_uniform(m1);
    const katai::mesh::Mesh* meshes[3] = {&m2, &m1, &M0.mesh};   // fine, medium, coarse

    MeshRun runs[3];
    for (int i = 0; i < 3; ++i) {
        runs[i] = run_on(pr, *meshes[i], probes, y_surface);
        std::printf("  level %d -> %d elements, %d nodes, h = %.4f m%s\n", 2 - i,
                    runs[i].elements, runs[i].nodes, runs[i].h, runs[i].ok ? "" : "  (FAILED)");
        check(runs[i].ok, "mesh " + std::to_string(i + 1) + " solves and probes");
    }
    if (!runs[0].ok || !runs[1].ok || !runs[2].ok) {
        std::printf("\n%d CHECK(S) FAILED\n", ++g_failures);
        return 1;
    }
    check(runs[0].h < runs[1].h && runs[1].h < runs[2].h,
          "the three meshes are ordered fine to coarse by measured h");

    // One quantity of interest at a time. `exact` is the closed form where one exists, NaN
    // where none does -- the study is still worth running there, it just cannot be scored
    // against an analytic answer.
    const auto study = [&](const std::string& name, const char* unit, double phi1, double phi2,
                           double phi3, double exact) {
        std::printf("\n-- %s --\n", name.c_str());
        gc::GridTriplet t;
        t.h1 = runs[0].h; t.h2 = runs[1].h; t.h3 = runs[2].h;
        t.phi1 = phi1; t.phi2 = phi2; t.phi3 = phi3;
        const bool have_exact = std::isfinite(exact) && std::fabs(exact) > 0.0;
        const auto dev = [&](double v) { return 100.0 * (v - exact) / std::fabs(exact); };
        if (have_exact) std::printf("   closed form           %.6g %s\n", exact, unit);
        std::printf("   fine   h = %.4f m  %.6g %s", t.h1, t.phi1, unit);
        if (have_exact) std::printf("  (%+.2f%%)", dev(t.phi1));
        std::printf("\n   medium h = %.4f m  %.6g %s", t.h2, t.phi2, unit);
        if (have_exact) std::printf("  (%+.2f%%)", dev(t.phi2));
        std::printf("\n   coarse h = %.4f m  %.6g %s", t.h3, t.phi3, unit);
        if (have_exact) std::printf("  (%+.2f%%)", dev(t.phi3));
        std::printf("\n");

        const gc::ConvergenceEstimate e = gc::grid_convergence_band(t);
        std::printf("   r21 = %.4f, r32 = %.4f, observed order p = %.4f (%s)\n", e.r21, e.r32, e.p,
                    gc::convergence_kind_name(e.kind));
        if (!e.message.empty()) std::printf("   note: %s\n", e.message.c_str());
        check(e.ok, name + ": the triplet yields an estimate");
        if (!e.ok) return;
        check(std::fabs(e.r21 - 2.0) < 1e-9 && std::fabs(e.r32 - 2.0) < 1e-9,
              name + ": the refinement factor is exactly 2 (nested refinement)");
        std::printf("   extrapolated          %.6g %s", e.phi_extrapolated, unit);
        if (have_exact) std::printf("  (%+.2f%%)", dev(e.phi_extrapolated));
        std::printf("\n   GCI from observed order = %.3f%%   asymptotic ratio %.4f\n",
                    100.0 * e.gci_fine, e.asymptotic_ratio);
        std::printf("   REPORTED BAND = +/- %.3f%%  [%s]\n", 100.0 * e.band,
                    e.band_basis.c_str());

        const double lo = std::fmin(t.phi1, std::fmin(t.phi2, t.phi3));
        const double hi = std::fmax(t.phi1, std::fmax(t.phi2, t.phi3));
        const double spread = (hi - lo) / std::fabs(t.phi1);
        std::printf("   spread across the three meshes = %.3f%% of the fine value\n",
                    100.0 * spread);

        // Whatever basis it was computed on, the published band has to cover the distance from
        // the fine mesh to the limit it estimates. That is what a band means.
        const double to_limit = std::fabs(t.phi1 - e.phi_extrapolated) / std::fabs(t.phi1);
        check(e.band >= to_limit * 0.999,
              name + ": the reported band covers the distance to the extrapolated limit");
        check(e.band < 0.05, name + ": the reported band is below 5%");

        if (e.asymptotic) {
            // In the asymptotic range the extrapolation is a better answer than any of the
            // three runs -- that is the claim of the method, and where a closed form exists
            // it IS checked rather than assumed.
            check(e.band == e.gci_fine, name + ": in the asymptotic range the band is the GCI");
            if (have_exact) {
                const double d_fine = std::fabs(t.phi1 - exact) / std::fabs(exact);
                const double d_ext = std::fabs(e.phi_extrapolated - exact) / std::fabs(exact);
                check(d_ext < d_fine,
                      name + ": asymptotic, so extrapolation moves TOWARD the closed form");
            }
        } else {
            // Not in the asymptotic range -- a finding about the STUDY, not a defect in the
            // runs. An unstructured mesher rebuilds the mesh from scratch at every density, so
            // successive meshes are not nested and the mesh-to-mesh irregularity can dominate
            // the h^p trend at practical densities. The contract is that this is named, that
            // the extrapolated value is not quoted, and that the band falls back to the
            // assumed order with the larger factor of safety.
            check(!e.band_basis.empty() && e.band_basis.find("assumed order") != std::string::npos,
                  name + ": the band falls back to the assumed order, and says so");
            check(e.safety_used == 3.0, name + ": with Roache's assumed-order factor of safety");
            // Outside the asymptotic range the three values are three answers rather than a
            // hierarchy of errors, so the band may not be narrower than their spread.
            check(e.band >= spread * 0.999,
                  name + ": and it covers the spread across the three meshes");
        }
        if (have_exact) {
            const double d_fine = std::fabs(t.phi1 - exact) / std::fabs(exact);
            check(d_fine < 0.04, name + ": the fine mesh is within 4% of the closed form");
            // The claim that makes the band worth publishing: the model deviation from the
            // closed form is inside the numerical uncertainty, or the band is honest about
            // not covering it. Both are stated; only silence would be wrong.
            std::printf("   model deviation on the fine mesh = %.3f%%, %s the reported band\n",
                        100.0 * d_fine, d_fine <= e.band ? "INSIDE" : "OUTSIDE");
        }
    };

    for (size_t k = 0; k < probes.size(); ++k) {
        char label[96];
        std::snprintf(label, sizeof(label), "sigma_z at %.1f m depth (recovered stress)",
                      probes[k]);
        study(label, "kPa", runs[0].sigma_z[k], runs[1].sigma_z[k], runs[2].sigma_z[k],
              -boussinesq_strip_sigma_z(q, a, probes[k]));   // compression negative
    }
    // No closed form: a uniform strip on an elastic HALF-PLANE has an unbounded surface
    // settlement, so the finite domain's value is a property of this model, not of Boussinesq.
    // The study is still the useful one -- it is the primary unknown, and its band is what a
    // settlement prediction would have to carry.
    study("settlement under the strip centre (primary unknown)", "m", runs[0].u_y, runs[1].u_y,
          runs[2].u_y, std::nan(""));

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
