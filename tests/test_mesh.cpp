// Structured tri6 mesh generator unit test (P0.4 geometry + P0.5 mesh).
// Verified: node/element counts, connectivity validity, element areas summing to
// the domain area (via tri6 detJ, positive orientation), boundary node sets.
#include <katai/fem/elements/tri6.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>

using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;
namespace tri6 = katai::core::tri6;

namespace {

int g_failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_failures;
    }
}

// Gather the element node coordinates into tri6::NodeCoords.
tri6::NodeCoords element_coords(const Mesh& m, int e) {
    tri6::NodeCoords nodes;
    for (int k = 0; k < tri6::kNodeCount; ++k) {
        const int n = m.node_of(e, k);
        nodes(k, 0) = m.x[n];
        nodes(k, 1) = m.y[n];
    }
    return nodes;
}

double element_area(const tri6::NodeCoords& nodes) {
    double area = 0.0;
    for (const auto& gp : tri6::gauss_points()) {
        const auto g = tri6::strain_displacement(nodes, gp.xi, gp.eta);
        area += gp.weight * g.det_jacobian;
    }
    return area;
}

void test_counts_and_connectivity() {
    const RectangularDomain domain{/*x0=*/1.0, /*y0=*/2.0, /*w=*/6.0,
                                   /*h=*/4.0, /*mat=*/7};
    const int nx = 3, ny = 2;
    const Mesh m = katai::mesh::generate_structured_tri6(domain, nx, ny);

    check(m.node_count == (2 * nx + 1) * (2 * ny + 1), "node count");
    check(m.element_count == 2 * nx * ny, "element count");
    check(static_cast<int>(m.connectivity.size()) == m.element_count * 6,
          "connectivity size");

    // All connectivity indices in the valid range.
    bool all_valid = true;
    for (int idx : m.connectivity)
        if (idx < 0 || idx >= m.node_count) all_valid = false;
    check(all_valid, "connectivity indices valid");

    // The material is right on every element.
    bool mat_ok = true;
    for (int mat : m.element_material)
        if (mat != domain.material_id) mat_ok = false;
    check(mat_ok, "element material id");

    // Boundary node counts (mid-edge included).
    check(static_cast<int>(m.bottom_nodes.size()) == 2 * nx + 1, "bottom nodes");
    check(static_cast<int>(m.top_nodes.size()) == 2 * nx + 1, "top nodes");
    check(static_cast<int>(m.left_nodes.size()) == 2 * ny + 1, "left nodes");
    check(static_cast<int>(m.right_nodes.size()) == 2 * ny + 1, "right nodes");
}

// Element areas sum to the domain area; every detJ positive (correct orientation).
void test_area_conservation() {
    const RectangularDomain domain{0.0, 0.0, 5.0, 3.0, 0};
    const Mesh m = katai::mesh::generate_structured_tri6(domain, 4, 3);

    double total = 0.0;
    bool positive = true;
    for (int e = 0; e < m.element_count; ++e) {
        const auto nodes = element_coords(m, e);
        const double a = element_area(nodes);
        if (a <= 0.0) positive = false;
        total += a;
    }
    check(positive, "all element areas positive (CCW)");
    check(std::fabs(total - domain.area()) < 1e-9, "area sum = domain area");
}

// nx=1, ny=1 → 2 elements, 9 nodes; the diagonal mid node (local 4) is shared.
void test_single_cell_shared_diagonal() {
    const RectangularDomain domain{0.0, 0.0, 1.0, 1.0, 0};
    const Mesh m = katai::mesh::generate_structured_tri6(domain, 1, 1);

    check(m.node_count == 9 && m.element_count == 2, "single-cell counts");
    // Lower element local-4 (mid 2-3) = the diagonal; upper element local-5 (mid 3-1).
    check(m.node_of(0, 4) == m.node_of(1, 5), "diagonal mid node shared");
}

} // namespace

int main() {
    test_counts_and_connectivity();
    test_area_conservation();
    test_single_cell_shared_diagonal();

    if (g_failures == 0) {
        std::printf("OK: mesh passed all checks\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
