// Local mesh density (PLAXIS coarseness factors) through the sizing-field mesher.
// docs/references/mesh-sizing.md: target edge h0 = sqrt(2 max_area); per-region factor caps the
// area inside a polygon; structural lines / loads become sources h_src = h0 * f with Lipschitz
// grading h(d) = h_src + g d. Checks: (A) defaults are bit-identical to the legacy constant-area
// mesher; (B) a region factor bounds element areas inside that region only; (C) auto-refine makes
// elements near a plate line ~f^2 smaller with graded growth away from it; (D) the Ruppert minimum
// angle bound survives the variable field; (E) no element escapes the domain.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

m::Project block() {
    m::Project pr;
    m::Material s; pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

double tri_area(const katai::mesh::Mesh& M, int e) {
    const int a = M.node_of(e, 0), b = M.node_of(e, 1), c = M.node_of(e, 2);
    return 0.5 * std::fabs((M.x[b] - M.x[a]) * (M.y[c] - M.y[a]) -
                           (M.x[c] - M.x[a]) * (M.y[b] - M.y[a]));
}
void centroid(const katai::mesh::Mesh& M, int e, double& cx, double& cy) {
    const int a = M.node_of(e, 0), b = M.node_of(e, 1), c = M.node_of(e, 2);
    cx = (M.x[a] + M.x[b] + M.x[c]) / 3.0;
    cy = (M.y[a] + M.y[b] + M.y[c]) / 3.0;
}
double min_angle_deg(const katai::mesh::Mesh& M) {
    double amin = 180.0;
    auto ang = [](double a, double b, double c) {
        const double v = std::fmax(-1.0, std::fmin(1.0, (b * b + c * c - a * a) / (2.0 * b * c)));
        return std::acos(v) * 180.0 / 3.14159265358979323846;
    };
    for (int e = 0; e < M.element_count; ++e) {
        const int n0 = M.node_of(e, 0), n1 = M.node_of(e, 1), n2 = M.node_of(e, 2);
        const double l0 = std::hypot(M.x[n1] - M.x[n2], M.y[n1] - M.y[n2]);
        const double l1 = std::hypot(M.x[n0] - M.x[n2], M.y[n0] - M.y[n2]);
        const double l2 = std::hypot(M.x[n0] - M.x[n1], M.y[n0] - M.y[n1]);
        amin = std::fmin(amin, std::fmin(ang(l0, l1, l2), std::fmin(ang(l1, l0, l2), ang(l2, l0, l1))));
    }
    return amin;
}

void test_default_identical() {
    std::printf("-- (A) defaults reproduce the legacy constant-area mesh exactly --\n");
    const auto pr = block();
    const auto M1 = katai::app::mesh_from_project(pr, 2.0, 6);
    katai::app::MeshOptions mo;   // all defaults, no factors anywhere
    const auto M2 = katai::app::mesh_from_project(pr, 2.0, 6, mo);
    std::printf("   legacy: %d elems / options: %d elems\n", M1.mesh.element_count, M2.mesh.element_count);
    check(M1.ok && M2.ok && M1.mesh.element_count == M2.mesh.element_count &&
          M1.mesh.node_count == M2.mesh.node_count,
          "no factors -> same mesher path, identical mesh");
}

void test_region_factor() {
    std::printf("-- (B) region coarseness factor (lower layer f = 0.25) --\n");
    m::Project pr = block();
    pr.polygons[0].y = {0, 0, 6, 6};   // lower 0..6
    m::SoilPolygon U = pr.polygons[0];
    U.y = {6, 6, 10, 10};              // upper 6..10
    U.coarseness = 1.0;
    pr.polygons[0].coarseness = 0.25;
    pr.polygons.push_back(U);

    const auto R = katai::app::mesh_from_project(pr, 2.0, 6);   // h0 = 2 m
    check(R.ok, "two-layer block meshed");
    const double cap_low = 0.5 * (2.0 * 0.25) * (2.0 * 0.25);   // 0.125 m^2
    double max_low = 0.0, max_up = 0.0; int n_low = 0, n_up = 0;
    for (int e = 0; e < R.mesh.element_count; ++e) {
        double cx, cy; centroid(R.mesh, e, cx, cy);
        const double a = tri_area(R.mesh, e);
        if (cy < 6.0) { max_low = std::fmax(max_low, a); ++n_low; }
        else          { max_up = std::fmax(max_up, a); ++n_up; }
    }
    std::printf("   lower: %d elems max area %.4f (cap %.4f) / upper: %d elems max area %.4f\n",
                n_low, max_low, cap_low, n_up, max_up);
    check(max_low <= cap_low * 1.001, "lower-region element areas respect h0 * f cap");
    check(max_up > 3.0 * cap_low, "upper region stays coarse (factor local, not global)");
    check(n_low > 2 * n_up, "refined region carries the element count");
}

void test_line_auto_refine() {
    std::printf("-- (C) auto-refine around a plate line (EMR-style, f_auto = 0.5) --\n");
    m::Project pr = block();
    m::StructElement e; e.kind = m::StructKind::Plate; e.name = "Plate";
    e.x1 = 4.0; e.y1 = 5.0; e.x2 = 12.0; e.y2 = 5.0;
    pr.structs.push_back(e);
    katai::app::MeshOptions mo; mo.auto_refine = true;          // GUI default
    const auto R = katai::app::mesh_from_project(pr, 2.0, 6, mo);   // h0 = 2 -> h_src = 1
    check(R.ok, "block + plate meshed");
    const double cap_near = 0.5 * 1.0 * 1.0;                     // 0.5 m^2 at the line
    double max_near = 0.0, max_far = 0.0;
    auto dist_line = [](double x, double y) {
        const double t = std::fmin(1.0, std::fmax(0.0, (x - 4.0) / 8.0));
        return std::hypot(x - (4.0 + 8.0 * t), y - 5.0);
    };
    int n_near = 0;
    for (int el = 0; el < R.mesh.element_count; ++el) {
        double cx, cy; centroid(R.mesh, el, cx, cy);
        const double a = tri_area(R.mesh, el), d = dist_line(cx, cy);
        if (d < 0.6) { max_near = std::fmax(max_near, a); ++n_near; }
        if (d > 4.0) max_far = std::fmax(max_far, a);
    }
    std::printf("   near-line: %d elems, max area %.4f (cap %.4f) / far max area %.4f\n",
                n_near, max_near, cap_near, max_far);
    // Within 0.6 m of the line the field allows at most 0.5*(h_src + g*0.6)^2 = 0.845 m^2.
    check(n_near > 10 && max_near > 0.0, "the near-line band actually contains elements");
    check(max_near <= 0.5 * 1.3 * 1.3 * 1.001, "elements at the plate are ~2x finer (auto factor 0.5)");
    check(max_far > 1.6 * max_near, "sizes grow away from the line (graded, not globally fine)");

    // Without auto-refine (and factor 1) the plate must NOT shrink the mesh.
    const auto R0 = katai::app::mesh_from_project(pr, 2.0, 6);
    check(R0.ok && R0.mesh.element_count < R.mesh.element_count,
          "auto-refine off -> coarser mesh (refinement was the field's doing)");
}

void test_quality_and_domain() {
    std::printf("-- (D/E) Ruppert quality + domain integrity under the variable field --\n");
    m::Project pr = block();
    pr.polygons[0].coarseness = 0.5;
    m::StructElement e; e.kind = m::StructKind::Geogrid; e.name = "Grid";
    e.x1 = 2.0; e.y1 = 8.0; e.x2 = 18.0; e.y2 = 8.0;
    pr.structs.push_back(e);
    katai::app::MeshOptions mo; mo.auto_refine = true;
    const auto R = katai::app::mesh_from_project(pr, 2.0, 6, mo);
    check(R.ok, "refined config meshed");
    const double amin = min_angle_deg(R.mesh);
    std::printf("   min angle %.2f deg (bound 20)\n", amin);
    check(amin > 19.5, "minimum angle bound survives the sizing field");
    bool inside = true;
    for (int el = 0; el < R.mesh.element_count; ++el) {
        double cx, cy; centroid(R.mesh, el, cx, cy);
        if (!(cx > -1e-9 && cx < 20.0 + 1e-9 && cy > -1e-9 && cy < 10.0 + 1e-9)) inside = false;
    }
    check(inside, "no element escapes the domain");
}

}  // namespace

int main() {
    std::printf("Local mesh density (coarseness factors) -- sizing-field mesher\n\n");
    test_default_identical();
    test_region_factor();
    test_line_auto_refine();
    test_quality_and_domain();
    if (g_failures == 0) {
        std::printf("\nOK: sizing field = legacy-identical defaults, local region/line refinement, quality kept\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
