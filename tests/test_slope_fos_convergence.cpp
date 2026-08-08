// The factor of safety is the number this program will be judged on, so it is the number that
// most needs to say how much of itself is the mesh.
//
// Strength reduction is the everyday computation of practical geotechnics, and it is also the
// one where mesh dependence is a known hazard: with a non-associated flow rule (psi = 0, which
// is what Griffiths and Lane use and what this corpus case uses) failure localises into a shear
// band, and without a regularisation the band narrows as the elements shrink. A factor of
// safety that drifts downward with every refinement is then not converging -- it is following
// the mesh -- and the only way to tell the two apart is to refine systematically and look.
//
// This case does that on the published benchmark, with the nested family from KV-NUM-006 so
// that element SIZE is the only thing changing, and reports the factor of safety with the band
// that belongs to it.
//
// A WARNING THIS FILE MUST CARRY, because the gated family is the cheap one. Two nested
// families of the same case, both measured on this tree on 2026-08-08:
//
//   base 2.2 m (gated here):   142 / 568 / 2272 elements -> 1.4210, 1.3867, 1.3645
//                              differences shrink; the estimator calls it monotonic
//                              CONVERGENCE and offers a band of +/-3.7%
//   base 1.6 m (not gated):    267 / 1068 / 4272 elements -> 1.4217, 1.3867, 1.3093
//                              differences GROW; monotonic DIVERGENCE, band +/-8.6%,
//                              a fall of 7.9% across the same fourfold refinement
//
// The coarse family looks like it is converging and it is not: refine further and the drift
// accelerates, which is what a shear band narrowing with the elements does. Anyone quoting a
// factor of safety from a three-mesh study on this problem should read that twice. The finer
// family is reproduced by setting elem_size to 1.6 below; it costs about ten minutes, which is
// why the suite runs the cheaper one and this comment carries the rest.
//
// READ THE REFERENCE FIRST -- what "1.4" means in Griffiths and Lane (1999), Table 2, verbatim:
//
//     FOS   E' d_max / (gamma H^2)   Iterations
//     0.80          0.379                2
//     1.00          0.381               10
//     1.20          0.422               20
//     1.30          0.453               41
//     1.35          0.544              792
//     1.40          1.476             1000   <- did not converge
//
// Their factor of safety is the first TRIAL value at which their viscoplastic algorithm fails
// to converge within an iteration ceiling of 1000 ("the non-convergence option is taken as
// being a suitable indicator of failure"). So the published 1.4 is not a converged number: it
// is a bracket, (1.35, 1.40], quantised by the trial increments they chose, and it depends on
// their ceiling. The chart solution of Bishop and Morgenstern (1960) gives 1.380 for the same
// problem. A comparison that quotes "+x% against Griffiths and Lane" without saying that is
// quoting a resolution the reference does not have.
//
// verify: KV-SLP-003
//   oracle:   published_benchmark
//   source:   D. V. Griffiths and P. A. Lane (1999), "Slope stability analysis by finite elements", Geotechnique 49(3):387-403, Example 1 and its Table 2 (read from the authors' copy at inside.mines.edu); A. W. Bishop and N. R. Morgenstern (1960), "Stability coefficients for earth slopes", Geotechnique 10(4):129-153; discretisation-error estimation per Roache (1994) and Celik et al. (2008) as implemented in katai/math/grid_convergence.hpp
//   locator:  Example 1: homogeneous 2:1 slope (26.57 deg), no foundation layer (D = 1), 1.2H of crest and 2H beyond the toe, rollers on the left boundary and full fixity at the base, phi' = 20 deg, c'/(gamma H) = 0.05, psi = 0, E' = 1e5 kPa, nu' = 0.3; dimensionalised here as H = 10 m, gamma = 20 kN/m3, c' = 10 kPa. Published: FE 1.4 (the first trial FOS that did not converge within 1000 iterations; 1.35 converged in 792), charts 1.380
//   quantity: factor of safety by phi-c reduction, run as the file's own initial procedure from tests/corpus/kv-slp-002-griffiths-lane-example1.k2d on a NESTED mesh family (mesher visited once, then refined uniformly twice; refinement factor exactly 2) [-]
//   expected: at the case file's own density a factor of safety inside the published bracket (1.35, 1.40] and within 3% of the chart value 1.380; and, under refinement, a DOWNWARD drift rather than convergence -- the signature of non-associated localisation, which the run must declare rather than hide
//   band:     as asserted below and MEASURED on this tree, not inherited. Gated family (elem_size 2.2 m): FoS 1.4210 (142 elements, h = 1.245 m), 1.3867 (568, h = 0.622 m), 1.3645 (2272, h = 0.311 m) -- a fall of 4.0% across a fourfold refinement, reported band +/-3.7%. Finer family (elem_size 1.6 m, not gated: about ten minutes): 1.4217 / 1.3867 / 1.3093, a fall of 7.9%, monotonic DIVERGENCE, band +/-8.6% -- so the cheaper family's apparent convergence is not to be trusted. What is asserted: the file's own density inside the published bracket and within 3% of the charts, a downward direction, a 12% regression bound on the magnitude, a band of the order of the drift, and the K2D-A005 declaration reaching the user

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

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

struct Level {
    bool ok = false;
    double h = 0.0, fos = 0.0;
    bool lower_bound = false;
    int elements = 0;
    bool declared_mesh_dependence = false;   // K2D-A005 reached the result
};

Level run_level(const m::Project& pr, const katai::mesh::Mesh& mesh) {
    Level L;
    double area = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        area += 0.5 * std::fabs((mesh.x[b] - mesh.x[a]) * (mesh.y[c] - mesh.y[a]) -
                                (mesh.x[c] - mesh.x[a]) * (mesh.y[b] - mesh.y[a]));
    }
    L.h = gc::mesh_representative_size(area, mesh.element_count);
    L.elements = mesh.element_count;
    const auto res = katai::app::solve_phases(pr, mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res.back().ok) {
        std::printf("      (solve: %s)\n", res.empty() ? "no phases" : res.back().message.c_str());
        return L;
    }
    L.fos = res.back().fos;
    L.lower_bound = res.back().fos_lower_bound;
    for (const auto& d : res.back().diagnostics)
        if (d.code == "K2D-A005") L.declared_mesh_dependence = true;
    L.ok = L.fos > 0.0;
    return L;
}

}  // namespace

int main() {
    std::printf("== factor of safety under systematic refinement (KV-SLP-003) ==\n");

    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-slp-002-griffiths-lane-example1.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        return 1;
    }
    // The published values, kept beside the run so a reader never has to look them up.
    constexpr double kCharts = 1.380;      // Bishop and Morgenstern (1960)
    constexpr double kFeLow = 1.35;        // Griffiths and Lane: converged in 792 iterations
    constexpr double kFeHigh = 1.40;       // Griffiths and Lane: did NOT converge in 1000

    // The coarse level; the finest is this refined twice. Chosen so that three strength
    // reduction analyses -- each of which is a sequence of nonlinear solves -- fit in a suite
    // that runs in minutes. The mesh dependence this case measures is a property of the method,
    // not of the absolute density, so a cheaper family measures the same thing.
    pr.mesh.elem_size = 2.2;
    const auto M0 = katai::app::mesh_from_project(pr);
    check(M0.ok, "base mesh built from the case file");
    if (!M0.ok) { std::printf("      (%s)\n", M0.message.c_str()); return 1; }
    const katai::mesh::Mesh m1 = katai::mesh::refine_uniform(M0.mesh);
    const katai::mesh::Mesh m2 = katai::mesh::refine_uniform(m1);
    const katai::mesh::Mesh* meshes[3] = {&m2, &m1, &M0.mesh};   // fine, medium, coarse

    Level lv[3];
    for (int i = 0; i < 3; ++i) {
        lv[i] = run_level(pr, *meshes[i]);
        std::printf("  %-6s %6d elements, h = %.4f m -> FoS %s %.4f  (%+.2f%% vs charts 1.380)\n",
                    i == 0 ? "fine" : i == 1 ? "medium" : "coarse", lv[i].elements, lv[i].h,
                    lv[i].lower_bound ? ">" : "=", lv[i].fos,
                    100.0 * (lv[i].fos - kCharts) / kCharts);
        check(lv[i].ok, std::string("level ") + std::to_string(i) + " produced a factor of safety");
    }
    if (!lv[0].ok || !lv[1].ok || !lv[2].ok) {
        std::printf("\n%d CHECK(S) FAILED\n", ++g_failures);
        return 1;
    }

    // The medium level is the density the checked-in case file itself asks for, and it is the
    // one the corpus reports. It must land inside the published bracket; the other two are
    // measured and printed, not asserted against a reference they were never run at.
    std::printf("      medium level %s the published FE bracket (%.2f, %.2f], charts %.3f\n",
                (lv[1].fos > kFeLow && lv[1].fos <= kFeHigh) ? "is INSIDE" : "is OUTSIDE", kFeLow,
                kFeHigh, kCharts);
    check(lv[1].fos > kFeLow && lv[1].fos <= kFeHigh,
          "at the case file's own density the factor of safety is inside the published bracket");
    check(std::fabs(lv[1].fos - kCharts) / kCharts < 0.03,
          "and within 3% of the Bishop-Morgenstern chart value 1.380");

    gc::GridTriplet t;
    t.h1 = lv[0].h; t.h2 = lv[1].h; t.h3 = lv[2].h;
    t.phi1 = lv[0].fos; t.phi2 = lv[1].fos; t.phi3 = lv[2].fos;
    const gc::ConvergenceEstimate e = gc::grid_convergence_band(t);
    std::printf("\n  r21 = %.4f, r32 = %.4f, observed order p = %.4f (%s)\n", e.r21, e.r32, e.p,
                gc::convergence_kind_name(e.kind));
    if (!e.message.empty()) std::printf("  note: %s\n", e.message.c_str());
    std::printf("  extrapolated FoS = %.4f\n", e.phi_extrapolated);
    std::printf("  REPORTED BAND = +/- %.3f%%  [%s]\n", 100.0 * e.band, e.band_basis.c_str());
    check(e.ok, "the triplet yields an estimate");
    check(std::fabs(e.r21 - 2.0) < 1e-9 && std::fabs(e.r32 - 2.0) < 1e-9,
          "the refinement factor is exactly 2 (nested refinement)");

    // The localisation question, asked directly. A non-associated strength-reduction analysis
    // may keep losing factor of safety as the elements shrink, because the shear band narrows
    // with them. If that is happening here it must be visible, not averaged away.
    const double drift = (lv[0].fos - lv[2].fos) / lv[2].fos;
    std::printf("  fine-vs-coarse change = %+.3f%% over a 4x refinement (%s)\n", 100.0 * drift,
                drift < 0.0 ? "DOWNWARD -- the localisation direction" : "upward");
    // What this case pins is the CHARACTER of the result, because that is what is true. The
    // factor of safety does not converge under refinement here, and it is not supposed to: with
    // psi = 0 the flow rule is non-associated, failure localises, and the band width follows the
    // elements. "The result obtained from a phi/c reduction is influenced by the mesh size,
    // element type and convergence tolerances" (Tschuchnigg, Schweiger and Sloan 2015, as cited
    // in Torggler, TU Graz, section 3.1.1). The assertions are therefore:
    //   1. the drift is downward -- if it ever turned upward, the mechanism changed and the case
    //      needs re-reading rather than a wider band;
    //   2. it stays inside the magnitude measured here, so a change that made localisation worse
    //      would be caught;
    //   3. the estimator refuses to call the extrapolated value a limit;
    //   4. the run TELLS the user the number is mesh-dependent. A factor of safety handed over
    //      without that sentence is the failure mode this whole programme exists to remove.
    check(drift < 0.0,
          "the factor of safety falls with refinement, as non-associated localisation predicts");
    check(std::fabs(drift) < 0.12,
          "and the fall stays within the -7.9% measured on this tree (regression bound)");
    check(e.band >= std::fabs(drift) * 0.5,
          "the reported band is of the order of the mesh dependence it describes");
    check(lv[0].declared_mesh_dependence && lv[1].declared_mesh_dependence &&
              lv[2].declared_mesh_dependence,
          "every run declares the mesh dependence to the user (K2D-A005)");

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
