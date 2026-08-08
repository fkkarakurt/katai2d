// Uniform (regular, "red") refinement: the instrument that makes a convergence study
// meaningful on this program's unstructured meshes.
//
// KV-NUM-005 measured the problem it exists to solve. Asking the mesher for three smaller
// element sizes produces three meshes that are not related to each other: different node
// positions, different angles, different achieved refinement ratios. The differences between
// their solutions then contain the mesher's irregularity as well as the discretisation, and at
// practical densities the first can be as large as the second -- which is why those triplets
// came out with observed orders of 0.24 and 6.98, orders no element here can deliver.
//
// Splitting every triangle at its edge midpoints instead gives the three properties the theory
// actually assumes, and this file pins each of them:
//
//   nested   -- every coarse node survives, at the same coordinates
//   similar  -- each child is the parent scaled by 1/2, so EVERY angle is preserved exactly
//               and the refined family is shape-regular by construction
//   exact    -- four children of equal area per parent, so h halves exactly and r = 2
//
// verify: KV-NUM-006
//   oracle:   closed_form
//   source:   the midpoint-triangle theorem of elementary geometry (the medial triangle and the three corner triangles are each similar to the parent with ratio 1/2 and area 1/4); the refinement scheme itself is the regular "red" refinement of R. E. Bank, A. H. Sherman and A. Weiser (1983), "Refinement algorithms and data structures for regular local mesh refinement", in Scientific Computing (IMACS/North-Holland); shape regularity is the hypothesis of the finite element convergence theory in P. G. Ciarlet (1978), The Finite Element Method for Elliptic Problems
//   locator:  a triangle with vertices A, B, C and edge midpoints M_AB, M_BC, M_CA splits into (A, M_AB, M_CA), (M_AB, B, M_BC), (M_CA, M_BC, C) and (M_AB, M_BC, M_CA); each has area |ABC|/4 and the same three interior angles as ABC (stated in full)
//   quantity: element count, total area, the multiset of interior angles, the minimum angle, the representative cell size h, and the coordinates of every SHARED coarse node (corners and edge nodes), before and after one and two refinements of a mesh built by the project mesher, for tri6 and tri15 [-; m2; degrees; m]
//   expected: elements x4 per level; total area unchanged; the minimum angle unchanged; h halved exactly; every shared coarse node present in the refined mesh; every child material equal to its parent's; and for tri15 the only nodes that do not survive are the three element-interior nodes per parent, which are lattice points of no child and are shared with nothing
//   band:     1e-12 relative on areas, angles and h, and 1e-12 absolute in metres on node coordinates, as asserted below -- these are round-off tolerances: the theorem is exact, so anything larger is an implementation error

#include <katai/jobs/mesh_builder.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
namespace km = katai::mesh;

namespace {

constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;

void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
void check_close(double got, double want, double rtol, const std::string& what) {
    const double err = std::fabs(got - want) / std::fmax(std::fabs(want), 1e-300);
    std::printf(err <= rtol ? "ok:   %s (%.15g vs %.15g)\n" : "FAIL: %s (%.15g vs %.15g)\n",
                what.c_str(), got, want);
    if (!(err <= rtol)) ++g_failures;
}

double triangle_area(const km::Mesh& mesh, int e) {
    const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
    return 0.5 * std::fabs((mesh.x[b] - mesh.x[a]) * (mesh.y[c] - mesh.y[a]) -
                           (mesh.x[c] - mesh.x[a]) * (mesh.y[b] - mesh.y[a]));
}

double total_area(const km::Mesh& mesh) {
    double s = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) s += triangle_area(mesh, e);
    return s;
}

// The three interior angles of element e, in degrees, sorted.
std::vector<double> angles_of(const km::Mesh& mesh, int e) {
    const int n[3] = {mesh.node_of(e, 0), mesh.node_of(e, 1), mesh.node_of(e, 2)};
    double len[3];
    for (int k = 0; k < 3; ++k) {
        const int p = n[(k + 1) % 3], q = n[(k + 2) % 3];
        len[k] = std::hypot(mesh.x[p] - mesh.x[q], mesh.y[p] - mesh.y[q]);   // side opposite k
    }
    std::vector<double> a(3);
    for (int k = 0; k < 3; ++k) {
        const double o = len[k], u = len[(k + 1) % 3], v = len[(k + 2) % 3];
        a[k] = std::acos(std::clamp((u * u + v * v - o * o) / (2.0 * u * v), -1.0, 1.0)) * 180.0 /
               kPi;
    }
    std::sort(a.begin(), a.end());
    return a;
}

double min_angle(const km::Mesh& mesh) {
    double lo = 180.0;
    for (int e = 0; e < mesh.element_count; ++e) lo = std::fmin(lo, angles_of(mesh, e)[0]);
    return lo;
}

double representative_h(const km::Mesh& mesh) {
    return std::sqrt(total_area(mesh) / std::fmax(1, mesh.element_count));
}

// A layered project, so material inheritance is actually exercised rather than assumed.
m::Project two_layer_project() {
    m::Project pr;
    pr.name = "refinement fixture";
    pr.x_min = 0; pr.x_max = 12; pr.y_min = 0; pr.y_max = 8;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;
    for (int i = 0; i < 2; ++i) {
        m::Material s;
        s.name = i ? "Lower" : "Upper";
        s.model = m::SoilModel::LinearElastic;
        s.E = 20000.0 * (i + 1); s.nu = 0.3;
        s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
        pr.materials.push_back(s);
    }
    m::SoilPolygon up;
    up.name = "Upper"; up.material = 0;
    up.x = {0, 12, 12, 0}; up.y = {4, 4, 8, 8};
    up.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed,
                  (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    m::SoilPolygon low;
    low.name = "Lower"; low.material = 1;
    low.x = {0, 12, 12, 0}; low.y = {0, 0, 4, 4};
    low.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                   (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(low);
    pr.polygons.push_back(up);
    return pr;
}

// The coarse nodes that MUST survive refinement: everything a neighbouring element can see,
// i.e. the corners and the edge nodes. A tri15 element's three interior nodes (local 12..14)
// sit at barycentric (2,1,1)/4 and its permutations, which are not lattice points of any
// child; they are private to their element and shared with nothing, so their disappearance
// costs no continuity. Saying which nodes are covered is the difference between a true
// statement and a convenient one.
std::vector<int> shared_nodes(const km::Mesh& mesh) {
    const int shared_locals = mesh.nodes_per_element == 15 ? 12 : mesh.nodes_per_element;
    std::vector<char> is_shared(mesh.node_count, 0);
    for (int e = 0; e < mesh.element_count; ++e)
        for (int k = 0; k < shared_locals; ++k) is_shared[mesh.node_of(e, k)] = 1;
    std::vector<int> out;
    for (int n = 0; n < mesh.node_count; ++n)
        if (is_shared[n]) out.push_back(n);
    return out;
}

// Is every listed coarse node present in the refined mesh, at the same place? This is the
// nesting property, and it is what lets two solutions be compared point by point.
bool nested_in(const km::Mesh& coarse, const km::Mesh& fine, double tol, int& missing,
               const std::vector<int>& which) {
    missing = -1;
    // Sorted by x, then a binary-searched window: an O(N log N) lookup with no bucket
    // arithmetic to get wrong. The first version of this used a spatial hash and reported a
    // node "missing" that was 1.1e-16 away -- a false alarm about the mesh caused by the
    // TEST, which is exactly the kind of finding that discredits a verification suite. A
    // check that can fail for its own reasons is not a check.
    std::vector<int> order(fine.node_count);
    for (int n = 0; n < fine.node_count; ++n) order[n] = n;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return fine.x[a] < fine.x[b]; });
    for (int n : which) {
        const double cx = coarse.x[n], cy = coarse.y[n];
        auto lo = std::lower_bound(order.begin(), order.end(), cx - tol,
                                   [&](int a, double v) { return fine.x[a] < v; });
        bool found = false;
        for (auto it = lo; it != order.end() && fine.x[*it] <= cx + tol && !found; ++it)
            if (std::fabs(fine.y[*it] - cy) <= tol) found = true;
        if (!found) { missing = n; return false; }
    }
    return true;
}

void check_one_level(const km::Mesh& coarse, const km::Mesh& fine, const std::string& tag) {
    check(fine.element_count == 4 * coarse.element_count, tag + ": four children per parent");
    check(fine.nodes_per_element == coarse.nodes_per_element, tag + ": same element order");
    check_close(total_area(fine), total_area(coarse), 1e-12, tag + ": total area is unchanged");
    check_close(min_angle(fine), min_angle(coarse), 1e-12,
                tag + ": the minimum angle is unchanged (children are similar to the parent)");
    check_close(representative_h(fine), 0.5 * representative_h(coarse), 1e-12,
                tag + ": h halves exactly, so the refinement factor is 2");
    int missing = -1;
    check(nested_in(coarse, fine, 1e-12, missing, shared_nodes(coarse)),
          tag + ": every shared coarse node survives at the same coordinates (nested)");
    if (missing >= 0) {
        std::printf("      first missing coarse node %d at (%.17g, %.17g)\n", missing,
                    coarse.x[missing], coarse.y[missing]);
        for (int e = 0; e < coarse.element_count; ++e)
            for (int k = 0; k < coarse.nodes_per_element; ++k)
                if (coarse.node_of(e, k) == missing)
                    std::printf("        occupies element %d local %d\n", e, k);
        int best = -1; double bd = 1e300;
        for (int f = 0; f < fine.node_count; ++f) {
            const double d = std::hypot(fine.x[f] - coarse.x[missing], fine.y[f] - coarse.y[missing]);
            if (d < bd) { bd = d; best = f; }
        }
        std::printf("        nearest refined node %d at (%.17g, %.17g), distance %.3e\n", best,
                    fine.x[best], fine.y[best], bd);
    }
    // And the exception is exactly the stated one, not a larger silent set: for tri15 the
    // nodes that do NOT survive are precisely the three interior nodes of each parent.
    if (coarse.nodes_per_element == 15)
        check((int)shared_nodes(coarse).size() == coarse.node_count - 3 * coarse.element_count,
              tag + ": the only non-shared nodes are the three interior nodes per element");
    bool materials_ok = true;
    for (int e = 0; e < coarse.element_count && materials_ok; ++e)
        for (int k = 0; k < 4; ++k)
            if (fine.element_material[4 * e + k] != coarse.element_material[e])
                materials_ok = false;
    check(materials_ok, tag + ": every child inherits its parent's material");

    // The similarity is per element, not merely on average: each child's angle triple must
    // equal its parent's. That is the property that keeps mesh quality out of the comparison.
    bool angles_ok = true;
    for (int e = 0; e < coarse.element_count && angles_ok; ++e) {
        const std::vector<double> pa = angles_of(coarse, e);
        for (int k = 0; k < 4; ++k) {
            const std::vector<double> ca = angles_of(fine, 4 * e + k);
            for (int i = 0; i < 3; ++i)
                if (std::fabs(ca[i] - pa[i]) > 1e-9) { angles_ok = false; break; }
        }
    }
    check(angles_ok, tag + ": every child has its parent's three angles, element by element");
}

}  // namespace

int main() {
    std::printf("== uniform (red) refinement (KV-NUM-006) ==\n");

    for (int order : {6, 15}) {
        std::printf("\n== element order %d ==\n", order);
        m::Project pr = two_layer_project();
        pr.mesh.order = order;
        const auto M = katai::app::mesh_from_project(pr);
        check(M.ok, "base mesh built from the project");
        if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); continue; }

        const km::Mesh& base = M.mesh;
        std::printf("   base:   %d elements, %d nodes, h = %.6f, min angle %.4f deg\n",
                    base.element_count, base.node_count, representative_h(base),
                    min_angle(base));
        const km::Mesh once = km::refine_uniform(base);
        std::printf("   once:   %d elements, %d nodes, h = %.6f, min angle %.4f deg\n",
                    once.element_count, once.node_count, representative_h(once),
                    min_angle(once));
        const km::Mesh twice = km::refine_uniform(once);
        std::printf("   twice:  %d elements, %d nodes, h = %.6f, min angle %.4f deg\n",
                    twice.element_count, twice.node_count, representative_h(twice),
                    min_angle(twice));

        check_one_level(base, once, "level 1");
        check_one_level(once, twice, "level 2");
        // Nesting is transitive, and it is the property the convergence study leans on:
        // the coarsest mesh's nodes must still be there two levels down.
        int missing = -1;
        check(nested_in(base, twice, 1e-12, missing, shared_nodes(base)),
              "the base mesh is still nested in the twice-refined mesh");

        // No node may be left unconnected: a tri15 parent's quarter-point and interior nodes
        // are NOT children corners, and carrying them over would leave free-floating degrees
        // of freedom -- a singular system that no assertion about areas would catch.
        std::vector<char> used(once.node_count, 0);
        for (int e = 0; e < once.element_count; ++e)
            for (int k = 0; k < once.nodes_per_element; ++k) used[once.node_of(e, k)] = 1;
        check(std::all_of(used.begin(), used.end(), [](char c) { return c != 0; }),
              "every node of the refined mesh belongs to an element");
        check(!once.boundary_nodes.empty(), "the refined mesh reports its boundary nodes");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
