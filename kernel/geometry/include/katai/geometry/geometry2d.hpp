#pragma once
// 2D planar arrangement: given a soup of line segments, compute the bounded faces (closed regions).
// This is how a modern FE pre-processor defines soil clusters: the user draws geometry lines and the
// closed areas they enclose become regions. Robust enough for engineering input (proper crossings +
// T-junctions; collinear overlaps are ignored). Output faces are CCW polygons (world y-up).
//
// Algorithm: split all segments at pairwise intersections -> planar graph (nodes + undirected edges)
// -> half-edge face traversal (next = clockwise neighbour of the twin at the shared node) -> keep
// the CCW (positive-area) cycles as bounded faces; the CW cycle is the unbounded outer face.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace katai::geometry {

struct V2 { double x = 0, y = 0; };

namespace detail {

// Proper/touching intersection of segments p->p2 and q->q2. Returns params t,u in [0,1].
inline bool seg_intersect(const V2& p, const V2& p2, const V2& q, const V2& q2,
                          double& t, double& u, double tol) {
    const double rx = p2.x - p.x, ry = p2.y - p.y;
    const double sx = q2.x - q.x, sy = q2.y - q.y;
    const double rxs = rx * sy - ry * sx;
    if (std::fabs(rxs) < 1e-12) return false;   // parallel / collinear (overlaps ignored)
    const double qpx = q.x - p.x, qpy = q.y - p.y;
    t = (qpx * sy - qpy * sx) / rxs;
    u = (qpx * ry - qpy * rx) / rxs;
    return t >= -tol && t <= 1.0 + tol && u >= -tol && u <= 1.0 + tol;
}

}  // namespace detail

inline std::vector<std::vector<V2>> arrangement_faces(const std::vector<std::array<V2, 2>>& segs,
                                                      double tol = 1e-7) {
    const int n = (int)segs.size();
    // 1) Split each segment at every intersection with the others.
    std::vector<std::vector<double>> params(n);
    for (int i = 0; i < n; ++i) { params[i].push_back(0.0); params[i].push_back(1.0); }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double t, u;
            if (detail::seg_intersect(segs[i][0], segs[i][1], segs[j][0], segs[j][1], t, u, 1e-9)) {
                params[i].push_back(std::clamp(t, 0.0, 1.0));
                params[j].push_back(std::clamp(u, 0.0, 1.0));
            }
        }

    // 2) Unique nodes (snap by tolerance) + undirected edges from consecutive split points.
    const double snap = 1e-4;
    std::map<std::pair<long long, long long>, int> node_of;
    std::vector<V2> nodes;
    auto node_id = [&](const V2& p) {
        const long long kx = (long long)std::llround(p.x / snap), ky = (long long)std::llround(p.y / snap);
        auto it = node_of.find({kx, ky});
        if (it != node_of.end()) return it->second;
        const int id = (int)nodes.size(); nodes.push_back(p); node_of[{kx, ky}] = id; return id;
    };
    std::vector<std::vector<int>> adj;  // neighbours per node (deduped)
    auto add_edge = [&](int a, int b) {
        if (a == b) return;
        if ((int)adj.size() <= std::max(a, b)) adj.resize(std::max(a, b) + 1);
        if (std::find(adj[a].begin(), adj[a].end(), b) == adj[a].end()) adj[a].push_back(b);
        if (std::find(adj[b].begin(), adj[b].end(), a) == adj[b].end()) adj[b].push_back(a);
    };
    for (int i = 0; i < n; ++i) {
        auto& pr = params[i];
        std::sort(pr.begin(), pr.end());
        pr.erase(std::unique(pr.begin(), pr.end(), [](double a, double b) { return std::fabs(a - b) < 1e-7; }), pr.end());
        const V2 a = segs[i][0], b = segs[i][1];
        int prev = -1;
        for (double tparam : pr) {
            const V2 p{a.x + tparam * (b.x - a.x), a.y + tparam * (b.y - a.y)};
            const int id = node_id(p);
            if (prev >= 0) add_edge(prev, id);
            prev = id;
        }
    }
    if (adj.empty()) return {};
    adj.resize(nodes.size());

    // 3) Sort each node's neighbours CCW by angle.
    auto ang = [&](int from, int to) { return std::atan2(nodes[to].y - nodes[from].y, nodes[to].x - nodes[from].x); };
    for (int v = 0; v < (int)nodes.size(); ++v)
        std::sort(adj[v].begin(), adj[v].end(), [&](int a, int b) { return ang(v, a) < ang(v, b); });

    // 4) Trace faces via the half-edge "next" rule (clockwise neighbour of the twin at the node).
    auto key = [](int u, int v) { return ((std::uint64_t)(std::uint32_t)u << 32) | (std::uint32_t)v; };
    std::unordered_map<std::uint64_t, bool> visited;
    std::vector<std::vector<V2>> faces;
    for (int u = 0; u < (int)nodes.size(); ++u)
        for (int v : adj[u]) {
            if (visited[key(u, v)]) continue;
            std::vector<int> cycle;
            int cu = u, cv = v;
            for (int guard = 0; guard < 100000; ++guard) {
                visited[key(cu, cv)] = true;
                cycle.push_back(cu);
                // at cv, next = clockwise neighbour of the twin (cv -> cu)
                const auto& out = adj[cv];
                int idx = (int)(std::find(out.begin(), out.end(), cu) - out.begin());
                const int nx = out[(idx - 1 + (int)out.size()) % out.size()];
                cu = cv; cv = nx;
                if (cu == u && cv == v) break;
            }
            // signed area (shoelace); keep CCW (positive) bounded faces
            double area2 = 0.0;
            for (size_t k = 0; k < cycle.size(); ++k) {
                const V2& a = nodes[cycle[k]]; const V2& b = nodes[cycle[(k + 1) % cycle.size()]];
                area2 += a.x * b.y - b.x * a.y;
            }
            if (area2 > 1e-7 && cycle.size() >= 3) {
                std::vector<V2> poly; poly.reserve(cycle.size());
                for (int id : cycle) poly.push_back(nodes[id]);
                faces.push_back(std::move(poly));
            }
        }
    return faces;
}

}  // namespace katai::geometry
