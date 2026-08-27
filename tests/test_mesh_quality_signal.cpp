// KV-DIA-002: a mesh that did not reach its own quality bound says so.
//
// Ruppert's refinement guarantees the requested minimum angle for the angles this tree asks for,
// so mesh quality here is STRUCTURAL rather than measured element by element -- and that is the
// right design. What was missing is the other half of any guard: the refinement loop carries a
// step cap as a final safety valve, and until this case existed the valve opened in silence. A
// mesh that met 20 degrees and a mesh that ran out of steps trying were the same object to
// everything downstream, and mesh quality is not something the answer is indifferent to.
//
// This is the same argument as K2D-A012 (the constitutive integrator running out of substeps) and
// the same shape: the guard was right, the silence was the defect.
//
// verify: KV-DIA-002
//   oracle:   independent_path
//   source:   Ruppert, J. (1995), "A Delaunay refinement algorithm for quality 2-dimensional mesh generation", Journal of Algorithms 18(3) 548-585 -- the termination guarantee holds for a minimum angle below about 20.7 deg, which is why this tree asks for 20 and why reaching the bound is the ordinary outcome rather than a hoped-for one. The implementation note is in kernel/jobs/include/katai/jobs/mesh_builder.hpp ("min_angle : Ruppert quality bound; <= ~20.7 deg guarantees termination")
//   locator:  the independent path is the mesh itself: the smallest angle over every produced triangle, computed here from the vertex coordinates with no help from the mesher, against the bound the mesher was given. A run that reports quality_met must satisfy that bound on EVERY element; the report and the geometry are checked against each other rather than the report being trusted
//   quantity: the minimum interior angle over all elements of a generated mesh [deg], and the quality_met flag the mesher returns beside it
//   expected: on an ordinary rectangular domain at a 20 deg bound, quality_met is true AND the measured minimum angle is >= 20 deg -- the two agree. The pair is the assertion: a flag that is always true proves nothing, so the same mesh is also measured
//   band:     the bound itself, with 1e-9 deg of arithmetic tolerance, as asserted below. This is not a discretisation band: Ruppert either met the bound on every element or it did not, and the flag says which
#include <katai/jobs/mesh_builder.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (ok) std::printf("ok:   %s\n", what.c_str());
    else { std::printf("FAIL: %s\n", what.c_str()); ++g_failures; }
}

// The smallest interior angle over every triangle, computed from the coordinates alone. This is
// the independent path: it asks the GEOMETRY what the mesher achieved rather than asking the
// mesher what it thinks it achieved.
double min_angle_deg(const katai::mesh::Triangulation& T) {
    const double pi = 3.14159265358979323846;
    double worst = 180.0;
    for (const auto& t : T.triangles) {
        const double ax = T.x[t[0]], ay = T.y[t[0]];
        const double bx = T.x[t[1]], by = T.y[t[1]];
        const double cx = T.x[t[2]], cy = T.y[t[2]];
        const double a = std::hypot(bx - cx, by - cy);
        const double b = std::hypot(ax - cx, ay - cy);
        const double c = std::hypot(ax - bx, ay - by);
        if (a <= 0.0 || b <= 0.0 || c <= 0.0) continue;
        const auto ang = [&](double o, double p, double q) {
            const double v = (p * p + q * q - o * o) / (2.0 * p * q);
            return std::acos(std::max(-1.0, std::min(1.0, v))) * 180.0 / pi;
        };
        worst = std::min({worst, ang(a, b, c), ang(b, a, c), ang(c, a, b)});
    }
    return worst;
}

void case_ordinary_mesh_meets_its_bound() {
    std::printf("\n== an ordinary domain: the flag and the geometry agree ==\n");
    // A plain rectangle: four corners, four segments. Ruppert's ordinary case.
    const std::vector<double> px{0.0, 20.0, 20.0, 0.0};
    const std::vector<double> py{0.0, 0.0, 10.0, 10.0};
    const std::vector<std::array<int, 2>> segs{{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    const double bound = 20.0;
    const katai::mesh::Triangulation T =
        katai::mesh::quality_mesh(px, py, segs, bound, 2.0, segs);

    check(!T.triangles.empty(), "the rectangle meshes");
    if (T.triangles.empty()) return;
    const double got = min_angle_deg(T);
    std::printf("  %zu elements, worst angle %.4f deg, %d refinement step(s), quality_met = %s\n",
                T.triangles.size(), got, T.refinement_steps, T.quality_met ? "true" : "false");
    check(T.quality_met, "the mesher reports it reached the bound");
    // The pair is the point: a flag that is always true proves nothing, so the SAME mesh is
    // measured independently and has to agree with the flag.
    check(got >= bound - 1e-9,
          "and the measured worst angle really is at or above the bound it reports meeting");
    check(T.refinement_steps > 0,
          "the refinement did work rather than returning on the first look");
}

void case_the_flag_is_carried_out_of_the_builder() {
    std::printf("\n== the signal survives the layer above, which is where a user could read it ==\n");
    m::Project pr;
    pr.x_min = 0.0; pr.x_max = 20.0; pr.y_min = 0.0; pr.y_max = 10.0;
    pr.materials.resize(1);
    m::SoilPolygon P;
    P.name = "Soil"; P.material = 0;
    P.x = {0.0, 20.0, 20.0, 0.0};
    P.y = {0.0, 0.0, 10.0, 10.0};
    P.edge_bc.assign(4, 0);
    pr.polygons = {P};
    pr.mesh.elem_size = 2.0;

    const auto R = katai::app::mesh_from_project(pr);
    check(R.ok, "the project meshes");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    std::printf("  %s\n", R.message.c_str());
    check(R.quality_met, "and the builder carries the quality answer out with it");
    check(R.min_angle_asked > 0.0, "...together with the bound that was asked for");
    // The message is what every front end already shows, so a met bound must NOT add noise to it.
    check(R.message.find("step cap") == std::string::npos,
          "a mesh that met its bound says nothing extra: the warning is for the case that did not");
}

}  // namespace

int main() {
    std::printf("== KV-DIA-002: a mesh reports whether it reached its own quality bound ==\n");
    case_ordinary_mesh_meets_its_bound();
    case_the_flag_is_carried_out_of_the_builder();
    if (g_failures == 0) {
        std::printf("\nOK: the quality bound is reported, and the report agrees with the geometry\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
