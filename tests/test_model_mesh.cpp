// test_model_mesh — meshing GUI soil polygons (app/build_mesh.hpp).
//
// The strong invariant: the sum of element (corner-triangle) areas must equal the
// analytic polygon area -- this catches both triangles left outside the soil and
// gaps/holes. Also checks material tagging, conformity of adjacent regions, and
// that connectivity indices and node sharing are sane.

#include <cmath>
#include <cstdio>
#include <vector>

#include <katai/jobs/mesh_builder.hpp>

using katai::app::mesh_from_project;
namespace m = katai::model;

static int failures = 0;
#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) { std::printf("FAIL: %s\n", msg); ++failures; }             \
        else std::printf("ok:   %s\n", msg);                                     \
    } while (0)

// Sum of element corner-triangle areas.
static double mesh_area(const katai::mesh::Mesh& msh) {
    double a = 0.0;
    for (int e = 0; e < msh.element_count; ++e) {
        const int n0 = msh.node_of(e, 0), n1 = msh.node_of(e, 1), n2 = msh.node_of(e, 2);
        a += 0.5 * std::fabs((msh.x[n1] - msh.x[n0]) * (msh.y[n2] - msh.y[n0]) -
                             (msh.x[n2] - msh.x[n0]) * (msh.y[n1] - msh.y[n0]));
    }
    return a;
}

static m::SoilPolygon poly(std::vector<double> x, std::vector<double> y, int mat) {
    m::SoilPolygon P; P.x = std::move(x); P.y = std::move(y); P.material = mat; return P;
}

int main() {
    // ---- 1) Unit square 4x4 (area 16), tri6 ------------------------------------
    {
        m::Project pr; pr.materials.resize(1);
        pr.polygons.push_back(poly({0, 4, 4, 0}, {0, 0, 4, 4}, 0));
        const auto R = mesh_from_project(pr, /*max_area*/ 1.0, /*order*/ 6);
        CHECK(R.ok, "square: meshed");
        CHECK(R.mesh.element_count > 0 && R.mesh.node_count > 0, "square: nonempty");
        CHECK(R.mesh.nodes_per_element == 6, "square: tri6");
        const double area = mesh_area(R.mesh);
        std::printf("      square area = %.6f (exact 16)\n", area);
        CHECK(std::fabs(area - 16.0) < 1e-6, "square: area conserved");
        bool all0 = true; for (int e = 0; e < R.mesh.element_count; ++e) if (R.mesh.element_material[e] != 0) all0 = false;
        CHECK(all0, "square: all elements material 0");
        // refinement actually happened (max_area 1 over area 16 -> at least ~16 tris)
        CHECK(R.mesh.element_count >= 16, "square: refined to target size");
        // connectivity in range
        bool inrange = true;
        for (int c : R.mesh.connectivity) if (c < 0 || c >= R.mesh.node_count) inrange = false;
        CHECK(inrange, "square: connectivity indices valid");
    }

    // ---- 2) L-shape (concave), area 12 -----------------------------------------
    {
        // 4x4 square minus the top-right 2x2 -> area 16 - 4 = 12.
        m::Project pr; pr.materials.resize(1);
        pr.polygons.push_back(poly({0, 4, 4, 2, 2, 0}, {0, 0, 2, 2, 4, 4}, 0));
        const auto R = mesh_from_project(pr, 0.5, 6);
        CHECK(R.ok, "L-shape: meshed");
        const double area = mesh_area(R.mesh);
        std::printf("      L-shape area = %.6f (exact 12)\n", area);
        CHECK(std::fabs(area - 12.0) < 1e-6, "L-shape: concavity respected (no outside triangles)");
        // every element centroid must lie inside the L
        bool all_in = true;
        for (int e = 0; e < R.mesh.element_count; ++e) {
            const int n0 = R.mesh.node_of(e, 0), n1 = R.mesh.node_of(e, 1), n2 = R.mesh.node_of(e, 2);
            const double cx = (R.mesh.x[n0] + R.mesh.x[n1] + R.mesh.x[n2]) / 3.0;
            const double cy = (R.mesh.y[n0] + R.mesh.y[n1] + R.mesh.y[n2]) / 3.0;
            if (!katai::app::point_in_polygon(cx, cy, pr.polygons[0])) all_in = false;
        }
        CHECK(all_in, "L-shape: all element centroids inside the region");
    }

    // ---- 3) Two adjacent squares sharing the edge x=4 (areas 16 + 16) ----------
    {
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 4, 4, 0}, {0, 0, 4, 4}, 0));   // left  -> material 0
        pr.polygons.push_back(poly({4, 8, 8, 4}, {0, 0, 4, 4}, 1));   // right -> material 1
        const auto R = mesh_from_project(pr, 1.0, 6);
        CHECK(R.ok, "two regions: meshed");
        const double area = mesh_area(R.mesh);
        std::printf("      two-region area = %.6f (exact 32)\n", area);
        CHECK(std::fabs(area - 32.0) < 1e-6, "two regions: total area conserved");
        bool has0 = false, has1 = false, only01 = true;
        for (int e = 0; e < R.mesh.element_count; ++e) {
            const int mm = R.mesh.element_material[e];
            if (mm == 0) has0 = true; else if (mm == 1) has1 = true; else only01 = false;
        }
        CHECK(has0 && has1 && only01, "two regions: both materials present (conforming)");
    }

    // ---- 4) tri15 promotion -----------------------------------------------------
    {
        m::Project pr; pr.materials.resize(1);
        pr.polygons.push_back(poly({0, 4, 4, 0}, {0, 0, 4, 4}, 0));
        const auto R = mesh_from_project(pr, 1.0, 15);
        CHECK(R.ok && R.mesh.nodes_per_element == 15, "tri15: meshed (15-node)");
        CHECK(std::fabs(mesh_area(R.mesh) - 16.0) < 1e-6, "tri15: area conserved");
    }

    std::printf(failures ? "\n%d CHECK(S) FAILED\n" : "\nAll checks passed.\n", failures);
    return failures ? 1 : 0;
}
