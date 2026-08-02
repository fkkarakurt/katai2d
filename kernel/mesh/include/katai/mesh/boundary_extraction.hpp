#pragma once
// Boundary extraction over a triangle mesh (Stage B5: extracted from the
// application driver -- both services are pure mesh geometry/topology and
// belong to the mesh module; nothing here knows a project, a material or a
// solver).

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

#include <katai/mesh/mesh.hpp>

namespace katai::mesh {

// Mesh nodes lying on the segment (x1,y1)-(x2,y2), ordered along it (corner, mid, corner, ...).
// The mesh was made to conform to structural lines, so these nodes exist exactly on the segment.
inline std::vector<int> collect_chain(const Mesh& mesh,
                                      double x1, double y1, double x2, double y2) {
    const double dx = x2 - x1, dy = y2 - y1, L2 = dx * dx + dy * dy;
    if (L2 < 1e-18) return {};
    const double L = std::sqrt(L2), tol = 1e-6 * L + 1e-9;
    std::vector<std::pair<double, int>> hits;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double t = ((mesh.x[n] - x1) * dx + (mesh.y[n] - y1) * dy) / L2;
        if (t < -1e-9 || t > 1.0 + 1e-9) continue;
        const double px = x1 + t * dx, py = y1 + t * dy;
        if (std::hypot(mesh.x[n] - px, mesh.y[n] - py) < tol) hits.push_back({t, n});
    }
    std::sort(hits.begin(), hits.end());
    std::vector<int> chain;
    double prev = -1e30;
    for (const auto& h : hits) if (h.first - prev > 1e-9) { chain.push_back(h.second); prev = h.first; }
    return chain;
}

// One domain-boundary element edge, for the dynamic absorbing / free-field boundaries. A boundary
// edge belongs to exactly ONE element; `node` holds its npe ordered nodes (corner, mid(s), corner --
// 3 for tri6, 5 for tri15), `material` the adjacent element's material, and (nx,ny) the outward unit
// normal (so lateral edges can be told from the base/top). See assemble_boundary_dashpot (D3).
struct BoundaryEdgeChain {
    std::array<int, 5> node{};
    int npe = 3;
    int material = 0;
    double nx = 0.0, ny = 0.0;
};

inline std::vector<BoundaryEdgeChain>
extract_boundary_edges(const Mesh& mesh, const std::vector<char>& active) {
    const int per = mesh.nodes_per_element;               // 6 or 15
    const int npe = per == 15 ? 5 : 3;                    // nodes along one edge
    // Local edge node ordering (corner -> corner), per element type (tri6.hpp / tri15.hpp).
    static const int e6[3][3]  = {{0, 3, 1}, {1, 4, 2}, {2, 5, 0}};
    static const int e15[3][5] = {{0, 3, 4, 5, 1}, {1, 6, 7, 8, 2}, {2, 9, 10, 11, 0}};
    static const int ec[3][2]  = {{0, 1}, {1, 2}, {2, 0}};   // corner pair of each local edge
    const auto on = [&](int e) { return active.empty() || active[e]; };
    // Count corner-edge occurrences: an edge shared by two elements is interior; boundary edges once.
    std::unordered_map<long long, int> count;
    const auto key = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        return static_cast<long long>(a) * mesh.node_count + b;
    };
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!on(e)) continue;
        for (int k = 0; k < 3; ++k)
            count[key(mesh.connectivity[e * per + ec[k][0]], mesh.connectivity[e * per + ec[k][1]])]++;
    }
    std::vector<BoundaryEdgeChain> out;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!on(e)) continue;
        const int c0 = mesh.connectivity[e * per + 0], c1 = mesh.connectivity[e * per + 1],
                  c2 = mesh.connectivity[e * per + 2];
        const double cx = (mesh.x[c0] + mesh.x[c1] + mesh.x[c2]) / 3.0;
        const double cy = (mesh.y[c0] + mesh.y[c1] + mesh.y[c2]) / 3.0;
        for (int k = 0; k < 3; ++k) {
            const int a = mesh.connectivity[e * per + ec[k][0]], b = mesh.connectivity[e * per + ec[k][1]];
            if (count[key(a, b)] != 1) continue;             // interior edge
            BoundaryEdgeChain bc; bc.npe = npe; bc.material = mesh.element_material[e];
            for (int i = 0; i < npe; ++i)
                bc.node[i] = mesh.connectivity[e * per + (npe == 5 ? e15[k][i] : e6[k][i])];
            double tx = mesh.x[b] - mesh.x[a], ty = mesh.y[b] - mesh.y[a];
            const double tl = std::hypot(tx, ty); if (tl < 1e-30) continue;
            double nx = ty / tl, ny = -tx / tl;              // unit normal (edge tangent rotated -90)
            const double mx = 0.5 * (mesh.x[a] + mesh.x[b]), my = 0.5 * (mesh.y[a] + mesh.y[b]);
            if (nx * (mx - cx) + ny * (my - cy) < 0) { nx = -nx; ny = -ny; }   // orient outward
            bc.nx = nx; bc.ny = ny;
            out.push_back(bc);
        }
    }
    return out;
}

} // namespace katai::mesh
