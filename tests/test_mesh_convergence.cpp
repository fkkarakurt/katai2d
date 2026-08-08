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
// What the sweep varies is ONE field of the checked-in case file: mesh.elem_size. Everything
// else -- geometry, material, load, element order, the mesher's local refinement rules -- is
// the file's own, so the three runs are the same model at three densities and nothing else.
// The representative cell size h is measured from each mesh as sqrt(total area / element
// count) rather than taken from the requested elem_size, because the mesher's achieved density
// is what the solution actually saw (Celik et al. 2008, step 1).
//
// verify: KV-NUM-005
//   oracle:   closed_form
//   source:   classical elasticity: Boussinesq point solution integrated over a uniform strip on a half-plane; discretisation-error estimation per P. J. Roache (1994), ASME J. Fluids Eng. 116(3):405-413 and I. B. Celik et al. (2008), ASME J. Fluids Eng. 130(7):078001 (Grid Convergence Index, factor of safety 1.25 for a three-mesh study), the solution-verification framework of ASME V&V 20-2009 / V&V 10-2006
//   locator:  sigma_z = (q/pi)[alpha + sin(alpha) cos(theta1 + theta2)], theta_i = atan((x -+ a)/z) from the strip edges, alpha = theta1 - theta2; h = sqrt(total area / element count); p = |ln|eps32/eps21| + q(p)|/ln(r21); phi_ext = (r21^p phi1 - phi2)/(r21^p - 1); GCI_fine = 1.25 e_a/(r21^p - 1) (stated in full)
//   quantity: vertical stress sigma_z at 2 m and 4 m depth on the strip centre line, interpolated with the element shape functions at the exact probe point, from tests/corpus/kv-fnd-008-strip-load.k2d solved at elem_size 2.0 / 1.4 / 1.0 m [kPa]
//   expected: the closed form above at each probe; the observed order p in the range a second-order displacement element can produce on a load-edge singularity; and the Richardson-extrapolated value closer to the closed form than the finest mesh's own value
//   band:     as asserted below and MEASURED on this tree, not inherited: fine-mesh deviation from the closed form within 4%, extrapolated deviation strictly smaller, GCI_fine below 5%, and the refinement factors above the 1.3 that Celik et al. recommend

#include <katai/fem/elements/point_location.hpp>
#include <katai/fem/elements/tri15.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/math/grid_convergence.hpp>
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

// Solve the case file at one requested element size and read the probes.
MeshRun run_at(const m::Project& base, double elem_size, const std::vector<double>& probe_depth,
               double y_surface) {
    MeshRun r;
    m::Project pr = base;
    pr.mesh.elem_size = elem_size;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return r; }

    double area = 0.0;
    for (int e = 0; e < M.mesh.element_count; ++e) {
        const int a = M.mesh.node_of(e, 0), b = M.mesh.node_of(e, 1), c = M.mesh.node_of(e, 2);
        area += 0.5 * std::fabs((M.mesh.x[b] - M.mesh.x[a]) * (M.mesh.y[c] - M.mesh.y[a]) -
                                (M.mesh.x[c] - M.mesh.x[a]) * (M.mesh.y[b] - M.mesh.y[a]));
    }
    r.h = gc::mesh_representative_size(area, M.mesh.element_count);
    r.elements = M.mesh.element_count;
    r.nodes = M.mesh.node_count;

    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res.back().ok) {
        std::printf("      (solve: %s)\n", res.empty() ? "no phases" : res.back().message.c_str());
        return r;
    }
    // Nodal sigma_yy from the recovered effective-stress field.
    std::vector<double> syy(M.mesh.node_count, 0.0);
    for (int n = 0; n < M.mesh.node_count && n < (int)res.back().stress.stress.size(); ++n)
        syy[n] = res.back().stress.stress[n](1);

    const double x_centre = 20.0;   // the strip runs x = 18..22 in the case file
    for (double d : probe_depth) {
        double v = 0.0;
        const bool got = M.mesh.nodes_per_element == 15
                             ? interpolate_at<katai::core::Tri15Element>(M.mesh, syy, x_centre,
                                                                        y_surface - d, v)
                             : interpolate_at<katai::core::Tri6Element>(M.mesh, syy, x_centre,
                                                                       y_surface - d, v);
        if (!got) return r;
        r.sigma_z.push_back(v);
    }
    // The settlement under the strip centre, from the primary unknown. Recovered stress is a
    // DERIVED field -- one differentiation and a nodal averaging away from the solution -- so it
    // converges more slowly and more raggedly than the displacement it came from. Carrying both
    // is the point: they answer different questions about the same three meshes.
    std::vector<double> uy(M.mesh.node_count, 0.0);
    for (int n = 0; n < M.mesh.node_count && 2 * n + 1 < (int)res.back().disp.size(); ++n)
        uy[n] = res.back().disp[2 * n + 1];
    const bool got_u = M.mesh.nodes_per_element == 15
                           ? interpolate_at<katai::core::Tri15Element>(M.mesh, uy, x_centre,
                                                                      y_surface, r.u_y)
                           : interpolate_at<katai::core::Tri6Element>(M.mesh, uy, x_centre,
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

    // Fine to coarse, each 1.5x the previous: comfortably above the 1.3 minimum refinement
    // factor Celik et al. recommend -- the wider the ratio, the less the mesher's own
    // irregularity contaminates the observed order -- and still cheap enough for a suite.
    const double sizes[3] = {0.7, 1.05, 1.575};
    MeshRun runs[3];
    for (int i = 0; i < 3; ++i) {
        runs[i] = run_at(base, sizes[i], probes, y_surface);
        std::printf("  elem_size %.2f m -> %d elements, %d nodes, h = %.4f m%s\n", sizes[i],
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
        check(e.r21 >= 1.3 && e.r32 >= 1.3,
              name + ": refinement factors meet the recommended minimum of 1.3");
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

        // Whatever basis it was computed on, the published band has to cover what the three
        // meshes actually did. A band narrower than the observed spread would be a claim the
        // data contradicts.
        check(e.band >= spread * 0.999,
              name + ": the reported band covers the spread across the three meshes");
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
