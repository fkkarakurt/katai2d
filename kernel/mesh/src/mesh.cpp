#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <unordered_map>

namespace katai::mesh {
namespace {
std::uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
}
} // namespace

Mesh tri6_from_triangulation(const Triangulation& tg, int material_id) {
    Mesh mesh;
    mesh.nodes_per_element = 6;
    mesh.x = tg.x;  // triangulation vertices become the corner nodes
    mesh.y = tg.y;

    std::unordered_map<std::uint64_t, int> midnode;  // edge -> mid-node index
    std::unordered_map<std::uint64_t, int> edge_use;  // edge -> incident triangles
    auto mid_of = [&](int a, int b) {
        const std::uint64_t k = edge_key(a, b);
        const auto it = midnode.find(k);
        if (it != midnode.end()) return it->second;
        const int idx = static_cast<int>(mesh.x.size());
        mesh.x.push_back(0.5 * (tg.x[a] + tg.x[b]));
        mesh.y.push_back(0.5 * (tg.y[a] + tg.y[b]));
        midnode.emplace(k, idx);
        return idx;
    };

    mesh.connectivity.reserve(tg.triangles.size() * 6);
    for (const auto& t : tg.triangles) {
        const int c0 = t[0], c1 = t[1], c2 = t[2];
        const int m01 = mid_of(c0, c1), m12 = mid_of(c1, c2), m20 = mid_of(c2, c0);
        // tri6 local ordering: corners 0,1,2; mid-edges on 0-1, 1-2, 2-0.
        mesh.connectivity.insert(mesh.connectivity.end(),
                                 {c0, c1, c2, m01, m12, m20});
        mesh.element_material.push_back(material_id);
        ++edge_use[edge_key(c0, c1)];
        ++edge_use[edge_key(c1, c2)];
        ++edge_use[edge_key(c2, c0)];
    }
    mesh.node_count = static_cast<int>(mesh.x.size());
    mesh.element_count = static_cast<int>(tg.triangles.size());

    // Boundary = nodes on edges used by exactly one triangle (corner + mid-edge).
    std::vector<char> on_boundary(mesh.node_count, 0);
    for (const auto& t : tg.triangles) {
        const int c[3] = {t[0], t[1], t[2]};
        for (int e = 0; e < 3; ++e) {
            const int a = c[e], b = c[(e + 1) % 3];
            if (edge_use[edge_key(a, b)] == 1) {
                on_boundary[a] = on_boundary[b] = 1;
                on_boundary[midnode[edge_key(a, b)]] = 1;
            }
        }
    }
    for (int n = 0; n < mesh.node_count; ++n)
        if (on_boundary[n]) mesh.boundary_nodes.push_back(n);

    return mesh;
}

Mesh tri15_from_triangulation(const Triangulation& tg, int material_id) {
    Mesh mesh;
    mesh.nodes_per_element = 15;
    mesh.x = tg.x;  // triangulation vertices become the corner nodes
    mesh.y = tg.y;

    // Per (undirected) edge: three nodes at 1/4, 2/4, 3/4 in canonical (lo->hi)
    // orientation, shared by adjacent triangles. Returned in the requested u->v
    // direction so the tri15 local edge ordering is consistent.
    std::unordered_map<std::uint64_t, std::array<int, 3>> edge_nodes;
    std::unordered_map<std::uint64_t, int> edge_use;
    auto nodes_on_edge = [&](int u, int v) -> std::array<int, 3> {
        const std::uint64_t k = edge_key(u, v);
        const int lo = std::min(u, v), hi = std::max(u, v);
        auto it = edge_nodes.find(k);
        std::array<int, 3> canon;
        if (it == edge_nodes.end()) {
            for (int s = 0; s < 3; ++s) {
                const double t = (s + 1) / 4.0;
                canon[s] = static_cast<int>(mesh.x.size());
                mesh.x.push_back(tg.x[lo] + t * (tg.x[hi] - tg.x[lo]));
                mesh.y.push_back(tg.y[lo] + t * (tg.y[hi] - tg.y[lo]));
            }
            edge_nodes.emplace(k, canon);
        } else {
            canon = it->second;
        }
        if (u == lo) return canon;                      // u->v is lo->hi
        return {canon[2], canon[1], canon[0]};          // reversed
    };

    mesh.connectivity.reserve(tg.triangles.size() * 15);
    for (const auto& t : tg.triangles) {
        const int c0 = t[0], c1 = t[1], c2 = t[2];
        const auto e01 = nodes_on_edge(c0, c1);
        const auto e12 = nodes_on_edge(c1, c2);
        const auto e20 = nodes_on_edge(c2, c0);
        // Three private interior nodes at barycentric (2,1,1),(1,2,1),(1,1,2)/4.
        int in[3];
        const int wa[3][3] = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
        for (int j = 0; j < 3; ++j) {
            in[j] = static_cast<int>(mesh.x.size());
            mesh.x.push_back((wa[j][0] * tg.x[c0] + wa[j][1] * tg.x[c1] +
                              wa[j][2] * tg.x[c2]) / 4.0);
            mesh.y.push_back((wa[j][0] * tg.y[c0] + wa[j][1] * tg.y[c1] +
                              wa[j][2] * tg.y[c2]) / 4.0);
        }
        // tri15 local ordering: 3 corners, edge 0-1, edge 1-2, edge 2-0, interior.
        mesh.connectivity.insert(
            mesh.connectivity.end(),
            {c0, c1, c2, e01[0], e01[1], e01[2], e12[0], e12[1], e12[2],
             e20[0], e20[1], e20[2], in[0], in[1], in[2]});
        mesh.element_material.push_back(material_id);
        ++edge_use[edge_key(c0, c1)];
        ++edge_use[edge_key(c1, c2)];
        ++edge_use[edge_key(c2, c0)];
    }
    mesh.node_count = static_cast<int>(mesh.x.size());
    mesh.element_count = static_cast<int>(tg.triangles.size());

    // Boundary = nodes on edges used by exactly one triangle (corners + 3 edge nodes).
    std::vector<char> on_boundary(mesh.node_count, 0);
    for (const auto& t : tg.triangles) {
        const int c[3] = {t[0], t[1], t[2]};
        for (int e = 0; e < 3; ++e) {
            const int a = c[e], b = c[(e + 1) % 3];
            if (edge_use[edge_key(a, b)] == 1) {
                on_boundary[a] = on_boundary[b] = 1;
                for (int idx : edge_nodes[edge_key(a, b)]) on_boundary[idx] = 1;
            }
        }
    }
    for (int n = 0; n < mesh.node_count; ++n)
        if (on_boundary[n]) mesh.boundary_nodes.push_back(n);

    return mesh;
}

Mesh generate_structured_tri6(const geometry::RectangularDomain& domain, int nx,
                              int ny) {
    assert(nx >= 1 && ny >= 1);

    // Fine node grid: quadratic element → intermediate nodes per cell.
    const int fine_nx = 2 * nx + 1;
    const int fine_ny = 2 * ny + 1;
    const auto fine_node = [fine_nx](int i, int j) { return j * fine_nx + i; };

    Mesh mesh;
    mesh.nodes_per_element = 6;
    mesh.node_count = fine_nx * fine_ny;
    mesh.element_count = 2 * nx * ny;

    // --- Node coordinates (uniformly spaced) ---
    mesh.x.resize(mesh.node_count);
    mesh.y.resize(mesh.node_count);
    const double dx = domain.width / (fine_nx - 1);
    const double dy = domain.height / (fine_ny - 1);
    for (int j = 0; j < fine_ny; ++j) {
        for (int i = 0; i < fine_nx; ++i) {
            const int n = fine_node(i, j);
            mesh.x[n] = domain.x0 + i * dx;
            mesh.y[n] = domain.y0 + j * dy;
        }
    }

    // --- Connectivity: each macro cell → 2 tri6 (corners CCW) ---
    // tri6 node order: c1, c2, c3, mid12, mid23, mid31.
    mesh.connectivity.reserve(mesh.element_count * 6);
    mesh.element_material.assign(mesh.element_count, domain.material_id);
    for (int cy = 0; cy < ny; ++cy) {
        for (int cx = 0; cx < nx; ++cx) {
            const int bi = 2 * cx;  // the cell's fine-grid base indices
            const int bj = 2 * cy;
            const auto g = [&](int a, int b) { return fine_node(bi + a, bj + b); };

            // Lower triangle: corners (0,0),(2,0),(0,2); diagonal (2,0)-(0,2).
            mesh.connectivity.push_back(g(0, 0));
            mesh.connectivity.push_back(g(2, 0));
            mesh.connectivity.push_back(g(0, 2));
            mesh.connectivity.push_back(g(1, 0));  // mid 1-2
            mesh.connectivity.push_back(g(1, 1));  // mid 2-3 (diagonal)
            mesh.connectivity.push_back(g(0, 1));  // mid 3-1

            // Upper triangle: corners (2,0),(2,2),(0,2); shares the diagonal.
            mesh.connectivity.push_back(g(2, 0));
            mesh.connectivity.push_back(g(2, 2));
            mesh.connectivity.push_back(g(0, 2));
            mesh.connectivity.push_back(g(2, 1));  // mid 1-2
            mesh.connectivity.push_back(g(1, 2));  // mid 2-3
            mesh.connectivity.push_back(g(1, 1));  // mid 3-1 (diagonal)
        }
    }

    // --- Boundary node sets (mid-edge nodes included) ---
    for (int i = 0; i < fine_nx; ++i) {
        mesh.bottom_nodes.push_back(fine_node(i, 0));
        mesh.top_nodes.push_back(fine_node(i, fine_ny - 1));
    }
    for (int j = 0; j < fine_ny; ++j) {
        mesh.left_nodes.push_back(fine_node(0, j));
        mesh.right_nodes.push_back(fine_node(fine_nx - 1, j));
    }

    return mesh;
}

Mesh generate_structured_tri15(const geometry::RectangularDomain& domain, int nx,
                               int ny) {
    assert(nx >= 1 && ny >= 1);

    // Quarter-node fine grid: quartic element → nodes at 1/4 spacing per cell.
    const int fine_nx = 4 * nx + 1;
    const int fine_ny = 4 * ny + 1;
    const auto fine_node = [fine_nx](int i, int j) { return j * fine_nx + i; };

    Mesh mesh;
    mesh.nodes_per_element = 15;
    mesh.node_count = fine_nx * fine_ny;
    mesh.element_count = 2 * nx * ny;

    mesh.x.resize(mesh.node_count);
    mesh.y.resize(mesh.node_count);
    const double dx = domain.width / (fine_nx - 1);
    const double dy = domain.height / (fine_ny - 1);
    for (int j = 0; j < fine_ny; ++j)
        for (int i = 0; i < fine_nx; ++i) {
            const int n = fine_node(i, j);
            mesh.x[n] = domain.x0 + i * dx;
            mesh.y[n] = domain.y0 + j * dy;
        }

    // Connectivity: each macro cell → 2 tri15. Node order as in tri15.hpp
    // (3 corners, edge 0-1, edge 1-2, edge 2-0, 3 interior). Shared edge/diagonal
    // nodes are shared automatically through the same grid index.
    mesh.connectivity.reserve(mesh.element_count * 15);
    mesh.element_material.assign(mesh.element_count, domain.material_id);
    for (int cy = 0; cy < ny; ++cy)
        for (int cx = 0; cx < nx; ++cx) {
            const int bi = 4 * cx, bj = 4 * cy;
            const auto g = [&](int a, int b) { return fine_node(bi + a, bj + b); };
            // Lower triangle: corners (0,0),(4,0),(0,4); diagonal (4,0)-(0,4).
            const int lower[15] = {
                g(0, 0), g(4, 0), g(0, 4),                  // corners
                g(1, 0), g(2, 0), g(3, 0),                  // edge 0-1
                g(3, 1), g(2, 2), g(1, 3),                  // edge 1-2 (diagonal)
                g(0, 3), g(0, 2), g(0, 1),                  // edge 2-0
                g(1, 1), g(2, 1), g(1, 2)};                 // interior
            // Upper triangle: corners (4,0),(4,4),(0,4); shares the diagonal.
            const int upper[15] = {
                g(4, 0), g(4, 4), g(0, 4),                  // corners
                g(4, 1), g(4, 2), g(4, 3),                  // edge 0-1
                g(3, 4), g(2, 4), g(1, 4),                  // edge 1-2
                g(1, 3), g(2, 2), g(3, 1),                  // edge 2-0 (diagonal)
                g(3, 2), g(3, 3), g(2, 3)};                 // interior
            for (int k = 0; k < 15; ++k) mesh.connectivity.push_back(lower[k]);
            for (int k = 0; k < 15; ++k) mesh.connectivity.push_back(upper[k]);
        }

    // Ordered boundary node sets (edge nodes included).
    for (int i = 0; i < fine_nx; ++i) {
        mesh.bottom_nodes.push_back(fine_node(i, 0));
        mesh.top_nodes.push_back(fine_node(i, fine_ny - 1));
    }
    for (int j = 0; j < fine_ny; ++j) {
        mesh.left_nodes.push_back(fine_node(0, j));
        mesh.right_nodes.push_back(fine_node(fine_nx - 1, j));
    }
    return mesh;
}

} // namespace katai::mesh
