// test_region_edit -- the region edits of the drawing tools (katai/model/region_edit.hpp).
//
//   close_region: a region drawn against its neighbours closes along THEIR boundary, so a shared
//     boundary is never traced twice. Before it, a region above a layer whose top has several
//     vertices had to be drawn through every one of them; drawn with its two end points only, the
//     straight closing edge left slivers of gap and overlap against the real boundary.
//   split_polygons: a line cuts every region it crosses into pieces that keep the parent's
//     material, coarseness, per-edge conditions and -- in every phase -- activation.
//   erase_polygon / insert_polygons: the phase activation flags move with the polygons. Deleting
//     a polygon in the Studio left them in place, which handed its excavation to the next one.
// The meshing oracle is the one of test_mesh_multi_region: total and per-region area, and the
// boundary length of the mesh against the perimeter of the union (an interior crack would add).

// verify: KV-GEO-003
//   oracle:   closed_form
//   source:   KATAI 2D input contract (docs/k2d-format.md, soil region object; per-edge conditions run vertex j -> j+1; a phase activates regions by index)
//   locator:  the shoelace area of the regions before and after each edit, evaluated by hand from the coordinates in this file; the space a line closes against a layer is the rectangle under the line minus the layer's area (stated in full: 10 x 8 - area(A))
//   quantity: area of each produced region [m2]; vertex count of the closed region [-]; mesh area and boundary length of the edited model [m2; m]; phase activation and per-edge conditions of every piece
//   expected: closed region 80 - 50.4 = 29.6 with 6 vertices, meshed area 80, boundary 36; roof 10; split pieces 20 + 30, meshed area 65, boundary 46; lens cut 44 + 6, meshed area 50, lens 6, boundary 30; pieces inactive wherever the parent was, the other region untouched
//   band:     1e-9 absolute, as asserted below -- round-off: the regions are built from the drawn coordinates and the edits add no approximation

#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/region_edit.hpp>

namespace m = katai::model;
using V = katai::geometry::V2;

static int failures = 0;
#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) { std::printf("FAIL: %s\n", msg); ++failures; }             \
        else std::printf("ok:   %s\n", msg);                                     \
    } while (0)

static double poly_area(const m::SoilPolygon& P) {   // even-odd area via the shoelace (keyholes cancel)
    double a2 = 0.0;
    for (size_t k = 0; k < P.x.size(); ++k) {
        const size_t j = (k + 1) % P.x.size();
        a2 += P.x[k] * P.y[j] - P.x[j] * P.y[k];
    }
    return 0.5 * a2;
}
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
    m::SoilPolygon P; P.x = std::move(x); P.y = std::move(y); P.material = mat;
    P.edge_bc.assign(P.x.size(), 0); return P;
}

int main() {
    // ---- 1) close against a layer whose top has four vertices ---------------------------------
    {
        m::Project pr; pr.materials.resize(2);
        // A: top runs (10,5) -> (6,5.2) -> (3,4.9) -> (0,5).
        pr.polygons.push_back(poly({0, 10, 10, 6, 3, 0}, {0, 0, 5, 5.2, 4.9, 5}, 0));
        const double aA = poly_area(pr.polygons[0]);
        // The user starts on A's top-left corner, clicks two free vertices, ends on A's top-right.
        const int n = m::close_region(pr, {V{0, 5}, V{0, 8}, V{10, 8}, V{10, 5}});
        CHECK(n == 1 && pr.polygons.size() == 2, "close: one region closed by the line and A's top");
        const auto& B = pr.polygons[1];
        std::printf("      B: %zu vertices, area %.9f (A %.9f)\n", B.x.size(), poly_area(B), aA);
        CHECK(B.x.size() == 6, "close: B follows all four vertices of A's top (2 drawn + 4 shared)");
        CHECK(std::fabs(poly_area(B) - (80.0 - aA)) < 1e-9, "close: B is exactly the space between the line and A");
        CHECK(B.edge_bc.size() == B.x.size() && B.material == -1, "close: a new region with no material and free edges");
        pr.polygons[1].material = 1;
        const auto R = katai::app::mesh_from_project(pr, 0.5, 6);
        CHECK(R.ok && R.quality_met, "close: the two regions mesh");
        CHECK(std::fabs(mesh_area(R.mesh) - 80.0) < 1e-9 && std::fabs(mesh_area(R.mesh, 1) - (80.0 - aA)) < 1e-9,
              "close: no gap and no overlap left between them (area 80, B exact)");
        CHECK(std::fabs(boundary_length(R.mesh) - 36.0) < 1e-9, "close: no interior crack (boundary = 36)");
    }
    // ---- 2) lines that close nothing -----------------------------------------------------------
    {
        m::Project pr; pr.materials.resize(1);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        CHECK(m::close_region(pr, {V{0, 5}, V{0, 8}, V{10, 8}}) == 0 && pr.polygons.size() == 1,
              "close: a line ending in the air closes nothing");
        CHECK(m::close_region(pr, {V{2, 1}, V{4, 4}, V{8, 1}}) == 0,
              "close: a line whose enclosed space is already soil closes nothing");
        // Across a gap: a line from A's top-left corner over to A's top-right corner closes the
        // space above A only.
        CHECK(m::close_region(pr, {V{0, 5}, V{5, 7}, V{10, 5}}) == 1 && std::fabs(poly_area(pr.polygons[1]) - 10.0) < 1e-9,
              "close: a roof over A's top closes the triangle between them (area 10)");
    }
    // ---- 3) split with phase flags and per-edge conditions -------------------------------------
    {
        m::Project pr; pr.materials.resize(2);
        m::SoilPolygon A = poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0);
        A.edge_bc[0] = 3; A.edge_bc[1] = 1;   // bottom and right edges carry conditions
        A.edge_flow.assign(4, 0); A.edge_flow[2] = 1; A.edge_head.assign(4, 0.0); A.edge_head[2] = 4.5;
        A.coarseness = 0.5; A.name = "Clay";
        pr.polygons.push_back(A);
        pr.polygons.push_back(poly({12, 15, 15, 12}, {0, 0, 5, 5}, 1));
        m::Phase dig; dig.poly_active = {0, 1};   // A excavated, B stays
        pr.phases.push_back(dig);
        const int cut = m::split_polygons(pr, {V{4, -1}, V{4, 6}});
        CHECK(cut == 1 && pr.polygons.size() == 3, "split: a line across A cuts it in two, B untouched");
        const double a0 = poly_area(pr.polygons[0]), a1 = poly_area(pr.polygons[1]);
        std::printf("      pieces: %.6f + %.6f\n", a0, a1);
        CHECK((std::fabs(a0 - 20) < 1e-9 && std::fabs(a1 - 30) < 1e-9) || (std::fabs(a0 - 30) < 1e-9 && std::fabs(a1 - 20) < 1e-9),
              "split: pieces of 20 and 30");
        CHECK(pr.polygons[0].material == 0 && pr.polygons[1].material == 0 && pr.polygons[1].coarseness == 0.5 &&
              pr.polygons[2].x[0] == 12.0, "split: pieces keep the material and coarseness; B moved behind them");
        CHECK(pr.phases[0].poly_active.size() == 3 && !pr.phases[0].active_poly(0) && !pr.phases[0].active_poly(1) &&
              pr.phases[0].active_poly(2), "split: both pieces excavated where A was, B still active");
        // Conditions follow the edges: every piece edge on y = 0 keeps 3, the right edge keeps 1,
        // the top keeps its head, and the cut edges (x = 4) start free.
        bool bc_ok = true;
        for (int pi = 0; pi < 2; ++pi) {
            const auto& P = pr.polygons[pi];
            for (size_t k = 0; k < P.x.size(); ++k) {
                const size_t j = (k + 1) % P.x.size();
                const bool bottom = P.y[k] == 0 && P.y[j] == 0, right = P.x[k] == 10 && P.x[j] == 10;
                const bool top = P.y[k] == 5 && P.y[j] == 5, cutedge = P.x[k] == 4 && P.x[j] == 4;
                if (bottom && P.edge_bc[k] != 3) bc_ok = false;
                if (right && P.edge_bc[k] != 1) bc_ok = false;
                if (top && (P.edge_flow[k] != 1 || P.edge_head[k] != 4.5)) bc_ok = false;
                if (cutedge && (P.edge_bc[k] != 0 || P.edge_flow[k] != 0)) bc_ok = false;
            }
        }
        CHECK(bc_ok, "split: per-edge conditions follow their edges; the cut edge starts free");
        const auto R = katai::app::mesh_from_project(pr, 0.5, 6);
        CHECK(R.ok && std::fabs(mesh_area(R.mesh) - 65.0) < 1e-9 && std::fabs(boundary_length(R.mesh) - 46.0) < 1e-9,
              "split: the pieces mesh as the original did (area 65, no crack)");
        CHECK(m::split_polygons(pr, {V{4, 1}, V{6, 1}}) == 0, "split: a line that does not cross a region cuts nothing");
    }
    // ---- 4) a closed line inside a region cuts out a lens --------------------------------------
    {
        m::Project pr; pr.materials.resize(1);
        pr.polygons.push_back(poly({0, 10, 10, 0}, {0, 0, 5, 5}, 0));
        const int cut = m::split_polygons(pr, {V{3, 1}, V{6, 1}, V{6, 3}, V{3, 3}, V{3, 1}});
        CHECK(cut == 1 && pr.polygons.size() == 2, "lens: a closed line inside cuts the region in two");
        const double a0 = std::fabs(poly_area(pr.polygons[0])), a1 = std::fabs(poly_area(pr.polygons[1]));
        std::printf("      pieces: %.6f + %.6f\n", a0, a1);
        CHECK((std::fabs(a0 - 44) < 1e-9 && std::fabs(a1 - 6) < 1e-9) || (std::fabs(a0 - 6) < 1e-9 && std::fabs(a1 - 44) < 1e-9),
              "lens: the ring (44, keyholed) and the lens (6)");
        pr.materials.resize(2);
        pr.polygons[a0 < a1 ? 0 : 1].material = 1;
        const auto R = katai::app::mesh_from_project(pr, 0.5, 6);
        CHECK(R.ok && std::fabs(mesh_area(R.mesh) - 50.0) < 1e-9 && std::fabs(mesh_area(R.mesh, 1) - 6.0) < 1e-9 &&
              std::fabs(boundary_length(R.mesh) - 30.0) < 1e-9, "lens: meshes as ring + lens, the keyhole bridge leaves no crack");
    }
    // ---- 5) erase keeps the flags on the right polygons ----------------------------------------
    {
        m::Project pr;
        for (int i = 0; i < 3; ++i) pr.polygons.push_back(poly({0.0 + i, 1.0 + i, 1.0 + i}, {0, 0, 1}, 0));
        m::Phase ph; ph.poly_active = {1, 0, 1}; pr.phases.push_back(ph);
        pr.initial.poly_active = {1, 1, 0};
        m::erase_polygon(pr, 0);
        CHECK(pr.polygons.size() == 2 && !pr.phases[0].active_poly(0) && pr.phases[0].active_poly(1) &&
              pr.initial.active_poly(0) && !pr.initial.active_poly(1),
              "erase: the excavated polygon is still the excavated one, in every phase");
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
