// The mesh is opened along the drawn walls and nowhere else.
//
// A plate with an interface is a wall: the mesh is split along its line, each node on it doubled,
// the elements on one side rewired to the twins. Where the line ends -- at a corner with another
// wall, at a junction, at a tip inside the ground -- the elements round the end node are
// classified, and a classification by the side of the line's infinite extension opens the ground
// beyond the end: two elements that shared an edge in the mesh no longer share it, along an edge
// no wall was drawn on -- a slit nobody drew.
//
// Two rules close it: round the end of a split line the elements take their side by the sector
// the walls leaving that point cut (an L corner, a vertex of a ring, a junction), and the two
// copies of a top that ends in the ground, where no other wall leaves it, are tied -- one point of
// the ground, which the wall's top meets through its interfaces on both sides.
//
// THE INVARIANT. Every edge two elements share in the mesh before the split is shared after it,
// unless it lies on a drawn wall's line. No closed form is needed: the topology is the oracle.
//
// verify: KV-STR-027
//   oracle:   closed_form
//   source:   KATAI 2D input contract (docs/k2d-format.md, structs[i].iface_pos / iface_neg): a plate with an interface separates the ground along its own line only
//   locator:  a linear elastic block 30 x 15 m under its own weight (gravity loading), tri6, element size 1 m, with walls with interfaces on both sides: an L corner (a wall down to a slab), a T junction (a slab ending on a wall), a slab inside the ground, a wall whose top is buried, a raft on the surface, a ring of 16 chords
//   quantity: the number of mesh edges two elements share before the split and not after it -- by a node of each that is not one point (two copies the run ties move as one) -- off every drawn wall line
//   expected: 0 in every configuration
//   band:     exact (a count; measured before 2026-10-26: 2 at the L corner either way round -- one under the corner, one past the slab's free end --, 1 at the T junction, the slab inside the ground, the buried top and the ring of 16 chords; the raft 0)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kPi = 3.14159265358979323846;

struct Seg { double x1, y1, x2, y2; };

m::Project block(const std::vector<Seg>& walls) {
    m::Project pr;
    pr.x_min = 0; pr.x_max = 30; pr.y_min = 0; pr.y_max = 15;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
    s.E = 50000.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 20.0;
    pr.materials.push_back(s);
    m::SoilPolygon g;
    g.name = "Ground"; g.material = 0; g.x = {0, 30, 30, 0}; g.y = {0, 0, 15, 15};
    g.edge_bc.assign(4, 0);
    pr.polygons.push_back(g);
    m::PlateMaterial pl;
    pl.name = "Wall"; pl.EA = 1.0e6; pl.EI = 1.0e4;
    pr.plates.push_back(pl);
    for (size_t i = 0; i < walls.size(); ++i) {
        m::StructElement e;
        e.kind = m::StructKind::Plate; e.name = "W" + std::to_string(i + 1);
        e.x1 = walls[i].x1; e.y1 = walls[i].y1; e.x2 = walls[i].x2; e.y2 = walls[i].y2;
        e.material = 0;
        e.iface_pos = e.iface_neg = true;
        pr.structs.push_back(e);
    }
    pr.initial.struct_active.assign(walls.size(), 1);
    return pr;
}

bool on_a_wall(double x, double y, const std::vector<Seg>& walls) {
    for (const auto& w : walls) {
        const double dx = w.x2 - w.x1, dy = w.y2 - w.y1, l2 = dx * dx + dy * dy;
        const double t = ((x - w.x1) * dx + (y - w.y1) * dy) / l2;
        if (t < -1e-9 || t > 1.0 + 1e-9) continue;
        if (std::hypot(w.x1 + t * dx - x, w.y1 + t * dy - y) < 1e-6) return true;
    }
    return false;
}

// The edges shared before the split and not after it, off every wall line; the first few printed.
int spurious_slits(const char* what, const std::vector<Seg>& walls) {
    const m::Project pr = block(walls);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { check(false, std::string(what) + ": meshes (" + M.message + ")"); return -1; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res[0].ok) { check(false, std::string(what) + ": solves"); return -1; }
    const auto& a = M.mesh;
    const auto& b = res[0].mesh;
    if (a.element_count != b.element_count || a.nodes_per_element != b.nodes_per_element) {
        check(false, std::string(what) + ": the split keeps the elements");
        return -1;
    }
    // Edge (two corner nodes, before the split) -> the elements that hold it, and where.
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edges;   // -> (element, local corner of the lower id)
    for (int e = 0; e < a.element_count; ++e)
        for (int k = 0; k < 3; ++k) {
            const int p = a.node_of(e, k), q = a.node_of(e, (k + 1) % 3);
            edges[{std::min(p, q), std::max(p, q)}].push_back({e, k});
        }
    const auto copy_of = [&](int e, int orig) {   // element e's node, after the split, where it had `orig`
        for (int k = 0; k < 3; ++k)
            if (a.node_of(e, k) == orig) return b.node_of(e, k);
        return -1;
    };
    // Two copies of a point that the run ties are one point of the ground: they move as one.
    const auto& u = res[0].disp;
    double umax = 0.0;
    for (double v : u) umax = std::max(umax, std::fabs(v));
    const auto one_point = [&](int p, int q) {
        return p == q || (std::fabs(u[2 * p] - u[2 * q]) <= 1e-12 * umax && std::fabs(u[2 * p + 1] - u[2 * q + 1]) <= 1e-12 * umax);
    };
    int count = 0;
    for (const auto& [ed, els] : edges) {
        if (els.size() != 2) continue;
        const int e = els[0].first, f = els[1].first;
        const bool apart = !one_point(copy_of(e, ed.first), copy_of(f, ed.first)) ||
                           !one_point(copy_of(e, ed.second), copy_of(f, ed.second));
        if (!apart) continue;
        const double mx = 0.5 * (a.x[ed.first] + a.x[ed.second]), my = 0.5 * (a.y[ed.first] + a.y[ed.second]);
        if (on_a_wall(mx, my, walls)) continue;
        if (count < 4)
            std::printf("      %s: opened off the walls between (%.3f, %.3f) and (%.3f, %.3f)\n", what,
                        a.x[ed.first], a.y[ed.first], a.x[ed.second], a.y[ed.second]);
        ++count;
    }
    std::printf("      %s: %d edge(s) opened off the walls\n", what, count);
    return count;
}

}  // namespace

int main() {
    std::printf("== the mesh is opened along the drawn walls and nowhere else (KV-STR-027) ==\n");
    struct Case { const char* what; std::vector<Seg> walls; };
    std::vector<Seg> ring;
    for (int k = 0; k < 16; ++k) {
        const double a0 = -0.5 * kPi + 2.0 * kPi * k / 16, a1 = -0.5 * kPi + 2.0 * kPi * (k + 1) / 16;
        ring.push_back({15 + 3 * std::cos(a0), 7 + 3 * std::sin(a0), 15 + 3 * std::cos(a1), 7 + 3 * std::sin(a1)});
    }
    const std::vector<Case> cases = {
        {"an L corner: a wall down to a slab", {{10, 15, 10, 8}, {10, 8, 22, 8}}},
        {"an L corner drawn the other way", {{22, 8, 10, 8}, {10, 8, 10, 15}}},
        {"a T junction: a slab ending on a wall", {{15, 15, 15, 5}, {15, 10, 25, 10}}},
        {"a slab inside the ground", {{5, 7, 12, 7}}},
        {"a wall whose top is buried", {{20, 4, 20, 11}}},
        {"a raft on the surface", {{6, 15, 14, 15}}},
        {"a ring of 16 chords", ring},
    };
    for (const auto& c : cases) {
        const int n = spurious_slits(c.what, c.walls);
        check(n == 0, std::string(c.what) + ": the mesh is opened along the walls only");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
