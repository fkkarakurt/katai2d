// test_mesh_multi_region -- soil regions that share, touch, overlap or nearly meet one another
// mesh as ONE conforming partition, in bounded time, with every region where it was drawn.
//
// Before the geometry was noded once to a model-scaled tolerance (katai/geometry/planar_graph.hpp)
// and pieces were kept by who owns their two sides:
//   - two regions overlapping in area meshed with the overlap as a HOLE (even-odd over every edge);
//   - a lens drawn inside a layer vanished and left a hole in its place;
//   - a corner typed 1e-5 off its neighbour, or a vertex 1e-4 above an edge, never finished
//     meshing (a 1e-5 sliver that Ruppert refinement tried to resolve; killed after 30 min);
//   - a single 3-degree corner, and a layer pinching out at 5.7 or 1.1 degrees, never finished
//     either (skinny corner triangles no insertion can improve, refined until the step cap).
// Oracles: total and per-material area against the drawn geometry (last polygon wins), and the
// BOUNDARY LENGTH of the mesh against the perimeter of the union -- an interior crack left by a
// non-conforming junction shows up as extra boundary edges, which the area alone would not see.

// verify: KV-GEO-002
//   oracle:   closed_form
//   source:   KATAI 2D input contract (docs/k2d-format.md, soil region object -- how regions combine: the last region containing a point owns it; edges are noded to 1e-5 of the model diagonal)
//   locator:  the shoelace area of each drawn polygon and the perimeter of their union, evaluated by hand from the coordinates in this file (stated in full: area = 1/2 |sum x_k y_k+1 - x_k+1 y_k|; overlap and lens areas by rectangle subtraction; the boundary of a conforming mesh is the set of element edges used once, whose length is the perimeter of the union)
//   quantity: total mesh area and mesh area per material [m2]; total length of the element edges used by one element only [m]
//   expected: share 80 = 50 + 30 / boundary 36; T-junction 68 = 50 + 18 / 36; overlap 88 = 38 + 50 / 44; lens 50 = 44 + 6 / 30; 2 mm gap 79.98 = 50 + 29.98 / 55.996; pinch-outs 50 = 45 + 5 and 49 + 1 / 30; 3-degree wedge 5 h / 10 + h + sqrt(100 + h^2); a corner 1e-5 off: identical to the exact model; a vertex 1e-4 above an edge: 68 / 36 within 1e-3
//   band:     1e-9 absolute, as asserted below -- round-off: the mesh follows the drawn edges exactly, so any larger deviation is a gap, an overlap, a lost region or an interior crack; 1e-6 for the 2 mm gap (its coordinates are not binary-exact); 1e-3 where a vertex 1e-4 off an edge is merged onto it by design (the edge bends through it, moving 5e-4 m2)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <katai/geometry/planar_graph.hpp>
#include <katai/jobs/mesh_builder.hpp>

using katai::app::mesh_from_project;
namespace m = katai::model;
namespace g = katai::geometry;

static int failures = 0;
#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) { std::printf("FAIL: %s\n", msg); ++failures; }             \
        else std::printf("ok:   %s\n", msg);                                     \
    } while (0)

static double mesh_area(const katai::mesh::Mesh& msh, int mat = -1) {
    double a = 0.0;
    for (int e = 0; e < msh.element_count; ++e) {
        if (mat >= 0 && msh.element_material[e] != mat) continue;
        const int n0 = msh.node_of(e, 0), n1 = msh.node_of(e, 1), n2 = msh.node_of(e, 2);
        a += 0.5 * std::fabs((msh.x[n1] - msh.x[n0]) * (msh.y[n2] - msh.y[n0]) -
                             (msh.x[n2] - msh.x[n0]) * (msh.y[n1] - msh.y[n0]));
    }
    return a;
}

// Total length of corner-triangle edges used by exactly one element.
static double boundary_length(const katai::mesh::Mesh& msh) {
    std::map<std::pair<int, int>, int> uses;
    for (int e = 0; e < msh.element_count; ++e)
        for (int k = 0; k < 3; ++k) {
            int a = msh.node_of(e, k), b = msh.node_of(e, (k + 1) % 3);
            if (a > b) std::swap(a, b);
            ++uses[{a, b}];
        }
    double L = 0.0;
    for (const auto& [ab, n] : uses)
        if (n == 1) L += std::hypot(msh.x[ab.second] - msh.x[ab.first], msh.y[ab.second] - msh.y[ab.first]);
    return L;
}

static m::SoilPolygon poly(std::vector<double> x, std::vector<double> y, int mat) {
    m::SoilPolygon P; P.x = std::move(x); P.y = std::move(y); P.material = mat; return P;
}

struct Run { katai::app::MeshResult R; double seconds; };
static Run run(const m::Project& pr) {
    const auto t0 = std::chrono::steady_clock::now();
    Run r{mesh_from_project(pr, 0.5, 6), 0.0};
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return r;
}

static void expect(const char* name, const m::Project& pr, double area, double mat0, double mat1,
                   double perimeter, double tol = 1e-9) {
    const Run r = run(pr);
    std::printf("      %-9s el=%d area=%.6f mat0=%.6f mat1=%.6f boundary=%.6f  %.3f s\n", name,
                r.R.mesh.element_count, mesh_area(r.R.mesh), mesh_area(r.R.mesh, 0),
                mesh_area(r.R.mesh, 1), boundary_length(r.R.mesh), r.seconds);
    char msg[160];
    std::snprintf(msg, sizeof msg, "%s: meshed, quality bound met", name);
    CHECK(r.R.ok && r.R.quality_met, msg);
    // Generous: every case here takes milliseconds; the defects this guards took forever.
    std::snprintf(msg, sizeof msg, "%s: meshed in bounded time (< 10 s)", name);
    CHECK(r.seconds < 10.0, msg);
    std::snprintf(msg, sizeof msg, "%s: total area = drawn union", name);
    CHECK(std::fabs(mesh_area(r.R.mesh) - area) < tol, msg);
    std::snprintf(msg, sizeof msg, "%s: each region has its drawn area (last polygon wins)", name);
    CHECK(std::fabs(mesh_area(r.R.mesh, 0) - mat0) < tol && std::fabs(mesh_area(r.R.mesh, 1) - mat1) < tol, msg);
    std::snprintf(msg, sizeof msg, "%s: boundary = perimeter of the union (no interior crack)", name);
    CHECK(std::fabs(boundary_length(r.R.mesh) - perimeter) < tol, msg);
}

static void planar_graph_checks() {
    using S = std::array<g::V2, 2>;
    // T-junction + a shared piece traced twice: A's top (0,5)-(10,5), B's bottom (0,5)-(6,5).
    {
        const std::vector<S> segs = {S{g::V2{0, 5}, g::V2{10, 5}}, S{g::V2{0, 5}, g::V2{6, 5}}};
        const auto G = g::planar_graph(segs, 1e-4);
        CHECK(G.nodes.size() == 3 && G.edges.size() == 2, "graph: T-junction splits the long segment, shared piece is one edge");
        bool shared = false;
        for (const auto& e : G.edges) if (e.src.size() == 2) shared = true;
        CHECK(shared, "graph: the shared piece remembers both sources");
        CHECK(G.chain[0].size() == 2 && G.chain[1].size() == 1, "graph: chains follow each segment");
        const auto c0 = G.chain_nodes(0);
        CHECK(c0.size() == 3 && G.nodes[c0[1]].x == 6.0, "graph: segment 0 walks through the junction in order");
    }
    // Near-misses below the tolerance merge; above it they stay.
    {
        const std::vector<S> segs = {S{g::V2{0, 0}, g::V2{10, 0}}, S{g::V2{10.00001, 0.00001}, g::V2{10, 5}},
                                     S{g::V2{4, 0.00005}, g::V2{4, 3}}};
        const auto G = g::planar_graph(segs, 1e-4);
        CHECK(G.nodes.size() == 5 && G.edges.size() == 4, "graph: a corner 1e-5 off merges; a vertex 5e-5 above an edge splits it");
        const auto G2 = g::planar_graph(segs, 1e-6);
        CHECK(G2.nodes.size() == 6, "graph: at a finer tolerance the same points stay apart");
    }
    // A crossing becomes a node.
    {
        const std::vector<S> segs = {S{g::V2{0, 0}, g::V2{2, 2}}, S{g::V2{0, 2}, g::V2{2, 0}}};
        const auto G = g::planar_graph(segs, 1e-6);
        CHECK(G.nodes.size() == 5 && G.edges.size() == 4, "graph: a crossing splits both segments");
    }
    // Faces: a box with a lens inside -> the lens face and the box face with one hole.
    {
        std::vector<S> segs;
        auto ring = [&](std::vector<g::V2> p) {
            for (size_t k = 0; k < p.size(); ++k) segs.push_back(S{p[k], p[(k + 1) % p.size()]});
        };
        ring({{0, 0}, {10, 0}, {10, 5}, {0, 5}});
        ring({{3, 1}, {6, 1}, {6, 3}, {3, 3}});
        const auto G = g::planar_graph(segs, 1e-6);
        const auto F = g::faces(G);
        bool box = false, lens = false;
        for (const auto& f : F) {
            if (std::fabs(f.area - 44.0) < 1e-9 && f.holes.size() == 1) box = true;
            if (std::fabs(f.area - 6.0) < 1e-9 && f.holes.empty()) lens = true;
        }
        CHECK(F.size() == 2 && box && lens, "faces: the box keeps the lens as a hole, the lens is its own face");
    }
    // simple_loops: a pinched walk comes apart into its two lobes.
    {
        const auto L = g::simple_loops({0, 1, 2, 3, 1, 4, 0});
        CHECK(L.size() == 2 && L[0].size() == 3 && L[1].size() == 3, "loops: a walk that revisits a node splits into simple loops");
    }
}

int main() {
    planar_graph_checks();

    {   // full shared edge, drawn by both
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({0, 10, 10, 0}, {5, 5, 8, 8}, 1));
        expect("share", pr, 80, 50, 30, 36);
    }
    {   // T-junction: B sits on part of A's top edge
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({0, 6, 6, 0}, {5, 5, 8, 8}, 1));
        expect("tjunc", pr, 68, 50, 18, 36);
    }
    {   // overlap in area: the later polygon wins it
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({4, 14, 14, 4}, {3, 3, 8, 8}, 1));
        expect("overlap", pr, 88, 38, 50, 44);
    }
    {   // a lens drawn inside a layer
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({3, 6, 6, 3}, {1, 1, 3, 3}, 1));
        expect("lens", pr, 50, 44, 6, 30);
    }
    {   // a 2 mm empty strip between two layers is real geometry and stays
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({0, 10, 10, 0}, {5.002, 5.002, 8, 8}, 1));
        expect("gap", pr, 79.98, 50, 29.98, 30 + 2 * 12.998, 1e-6);
    }
    {   // a layer pinching out at 5.7 and at 1.1 degrees
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 4}, 0));
        pr.polygons.push_back(poly({0, 10, 0}, {4, 5, 5}, 1));
        expect("pinch", pr, 50, 45, 5, 30);
        const Run r = run(pr);
        CHECK(r.R.corner_elements > 0 && r.R.message.find("narrower") != std::string::npos,
              "pinch: the elements kept at the pinch angle are counted and reported");
        m::Project p1; p1.materials.resize(2);
        p1.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 4.8}, 0));
        p1.polygons.push_back(poly({0, 10, 0}, {4.8, 5, 5}, 1));
        expect("pinch1", p1, 50, 49, 1, 30);
    }
    {   // one region with a 3-degree corner
        m::Project pr; pr.materials.resize(1);
        const double h = 10 * std::tan(3 * 3.14159265358979323846 / 180);
        pr.polygons.push_back(poly({0, 10, 10}, {0, 0, h}, 0));
        expect("sharp", pr, 5 * h, 5 * h, 0, 10 + h + std::hypot(10, h));
    }
    {   // a corner typed 1e-5 off its neighbour meshes exactly as the clean model
        m::Project clean; clean.materials.resize(2);
        clean.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        clean.polygons.push_back(poly({0, 10, 10, 0}, {5, 5, 8, 8}, 1));
        m::Project off = clean;
        off.polygons[1] = poly({0, 10.00001, 10, 0}, {5, 5.00001, 8, 8}, 1);
        expect("near", off, 80, 50, 30, 36);
        const Run a = run(clean), b = run(off);
        CHECK(a.R.mesh.x == b.R.mesh.x && a.R.mesh.y == b.R.mesh.y &&
              a.R.mesh.connectivity == b.R.mesh.connectivity,
              "near: the mesh is bit-identical to the one of the corner typed exactly");
    }
    {   // a vertex 1e-4 above an edge joins it (the edge bends through it by 1e-4)
        m::Project pr; pr.materials.resize(2);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        pr.polygons.push_back(poly({0, 6, 6, 0}, {5, 5.0001, 8, 8}, 1));
        expect("nearedge", pr, 68, 50, 18, 36, 1e-3);
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
