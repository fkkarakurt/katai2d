#pragma once
// A planar graph from a soup of segments, NODED to one geometric tolerance.
//
// The drawn geometry of a model is a set of independent polygons and lines. Where two of them
// are meant to meet they rarely meet to the last bit: a corner is typed 1e-5 off its neighbour,
// a vertex lands a hair above an edge, a shared boundary is traced twice. Left as it is, such a
// near-miss is a real feature of the input -- a sliver 0.01 mm thick and 10 m long -- and a
// quality mesher is obliged to resolve it, which in practice means it never finishes. This is
// the one place those near-misses are decided, so every consumer (the mesher, the drawing
// tools) sees the same topology.
//
// Construction, with eps the tolerance:
//   1) every segment end point becomes a node; a point within eps of an existing node IS that
//      node (first come keeps its coordinates), so nodes are pairwise at least eps apart;
//   2) every proper crossing of two segments becomes a node (same snapping);
//   3) every node lying within eps of a segment's interior splits that segment -- this is
//      what turns a T-junction, a collinear overlap and a near-touch into shared topology;
//   4) the pieces between consecutive nodes along a segment are the graph edges, and a piece
//      produced by more than one segment is ONE edge that remembers all its sources.
// `chain[i]` lists, in order from segment i's start to its end, the edges segment i became,
// with the direction each one is walked in -- so a polygon can be re-read on the graph.
//
// Faces: `faces()` traces the bounded faces of the graph (each an outer ring plus its holes),
// which is what the drawing tools use to close a region against existing boundaries.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <katai/geometry/geometry2d.hpp>   // V2

namespace katai::geometry {

// The tolerance a model's geometry is noded to: a fixed fraction of the extent of the drawing.
// 1e-5 of the diagonal is 1 mm on a 100 m model -- far below anything drawn on purpose in a
// geotechnical section, far above the round-off of any coordinate typed or snapped.
inline double geometry_tolerance(double minx, double miny, double maxx, double maxy) {
    const double d = std::hypot(maxx - minx, maxy - miny);
    return d > 0.0 && std::isfinite(d) ? 1e-5 * d : 1e-9;
}

struct PlanarGraph {
    struct Edge {
        int a = -1, b = -1;        // node indices, a != b
        std::vector<int> src;      // input segments that produced this edge (ascending)
    };
    struct Step { int edge; bool forward; };   // forward: walked a -> b
    std::vector<V2> nodes;
    std::vector<Edge> edges;
    std::vector<std::vector<Step>> chain;      // per input segment
    double eps = 0.0;

    // Node sequence walked by input segment i (its start node first). Empty when the segment
    // collapsed to a single node.
    std::vector<int> chain_nodes(int i) const {
        std::vector<int> out;
        for (const Step& s : chain[i]) {
            const int u = s.forward ? edges[s.edge].a : edges[s.edge].b;
            const int w = s.forward ? edges[s.edge].b : edges[s.edge].a;
            if (out.empty()) out.push_back(u);
            out.push_back(w);
        }
        return out;
    }
};

namespace detail {

inline double dist_point_segment(const V2& p, const V2& a, const V2& b, double& t) {
    const double dx = b.x - a.x, dy = b.y - a.y, l2 = dx * dx + dy * dy;
    t = l2 > 0.0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / l2 : 0.0;
    const double tc = std::clamp(t, 0.0, 1.0);
    return std::hypot(p.x - (a.x + tc * dx), p.y - (a.y + tc * dy));
}

}  // namespace detail

inline PlanarGraph planar_graph(const std::vector<std::array<V2, 2>>& segs, double eps) {
    PlanarGraph G;
    G.eps = eps > 0.0 ? eps : 1e-9;
    const double e = G.eps;
    const int n = (int)segs.size();

    // Node snapping through a hash grid of cell eps: a point snaps to the lowest-index node
    // within eps (deterministic, independent of hash order).
    std::unordered_map<std::uint64_t, std::vector<int>> grid;
    auto cell = [&](double v) { return (std::int64_t)std::floor(v / e); };
    auto ckey = [](std::int64_t cx, std::int64_t cy) {
        return ((std::uint64_t)(std::uint32_t)cx << 32) | (std::uint32_t)cy;
    };
    auto node_id = [&](const V2& p) {
        const std::int64_t cx = cell(p.x), cy = cell(p.y);
        int best = -1;
        for (std::int64_t i = cx - 1; i <= cx + 1; ++i)
            for (std::int64_t j = cy - 1; j <= cy + 1; ++j) {
                auto it = grid.find(ckey(i, j));
                if (it == grid.end()) continue;
                for (int id : it->second)
                    if (std::hypot(G.nodes[id].x - p.x, G.nodes[id].y - p.y) < e &&
                        (best < 0 || id < best))
                        best = id;
            }
        if (best >= 0) return best;
        G.nodes.push_back(p);
        const int id = (int)G.nodes.size() - 1;
        grid[ckey(cx, cy)].push_back(id);
        return id;
    };

    // 1) End points. The end is formed as a + (b - a), the way a split parameter of 1 would
    // place it, so a model with nothing to snap reproduces the coordinates it always had.
    std::vector<int> na(n), nb(n);
    for (int i = 0; i < n; ++i) {
        const V2 a = segs[i][0];
        const V2 b{a.x + (segs[i][1].x - a.x), a.y + (segs[i][1].y - a.y)};
        na[i] = node_id(a);
        nb[i] = node_id(b);
    }
    auto P = [&](int i, int end) { return G.nodes[end ? nb[i] : na[i]]; };

    // Per-segment bounding boxes, inflated by eps, for cheap rejection.
    struct Box { double x0, y0, x1, y1; };
    std::vector<Box> box(n);
    for (int i = 0; i < n; ++i) {
        const V2 a = P(i, 0), b = P(i, 1);
        box[i] = {std::min(a.x, b.x) - e, std::min(a.y, b.y) - e,
                  std::max(a.x, b.x) + e, std::max(a.y, b.y) + e};
    }
    auto overlap = [&](const Box& p, const Box& q) {
        return p.x0 <= q.x1 && q.x0 <= p.x1 && p.y0 <= q.y1 && q.y0 <= p.y1;
    };

    // 2) Proper crossings.
    for (int i = 0; i < n; ++i) {
        if (na[i] == nb[i]) continue;
        for (int j = i + 1; j < n; ++j) {
            if (na[j] == nb[j] || !overlap(box[i], box[j])) continue;
            if (na[i] == na[j] || na[i] == nb[j] || nb[i] == na[j] || nb[i] == nb[j]) continue;
            const V2 p = P(i, 0), p2 = P(i, 1), q = P(j, 0), q2 = P(j, 1);
            const double rx = p2.x - p.x, ry = p2.y - p.y, sx = q2.x - q.x, sy = q2.y - q.y;
            const double rxs = rx * sy - ry * sx;
            if (std::fabs(rxs) <= 1e-12 * std::hypot(rx, ry) * std::hypot(sx, sy)) continue;
            const double qpx = q.x - p.x, qpy = q.y - p.y;
            const double t = (qpx * sy - qpy * sx) / rxs, u = (qpx * ry - qpy * rx) / rxs;
            if (t < 0.0 || t > 1.0 || u < 0.0 || u > 1.0) continue;
            node_id({p.x + t * rx, p.y + t * ry});   // step 3 attaches it to both segments
        }
    }

    // 3) Every node within eps of a segment's interior splits it.
    std::vector<std::vector<std::pair<double, int>>> on(n);
    for (int i = 0; i < n; ++i) {
        if (na[i] == nb[i]) continue;
        const V2 a = P(i, 0), b = P(i, 1);
        for (int v = 0; v < (int)G.nodes.size(); ++v) {
            if (v == na[i] || v == nb[i]) continue;
            const V2& q = G.nodes[v];
            if (q.x < box[i].x0 || q.x > box[i].x1 || q.y < box[i].y0 || q.y > box[i].y1) continue;
            double t;
            if (detail::dist_point_segment(q, a, b, t) < e && t > 0.0 && t < 1.0)
                on[i].push_back({t, v});
        }
    }

    // 4) Edges between consecutive nodes along each segment, shared pieces merged.
    std::unordered_map<std::uint64_t, int> edge_of;
    G.chain.assign(n, {});
    for (int i = 0; i < n; ++i) {
        if (na[i] == nb[i]) continue;
        auto& l = on[i];
        std::sort(l.begin(), l.end());
        std::vector<int> seq{na[i]};
        for (const auto& [t, v] : l) if (v != seq.back()) seq.push_back(v);
        if (nb[i] != seq.back()) seq.push_back(nb[i]);
        for (size_t k = 0; k + 1 < seq.size(); ++k) {
            const int u = seq[k], w = seq[k + 1];
            if (u == w) continue;
            const std::uint64_t key = ((std::uint64_t)(std::uint32_t)std::min(u, w) << 32) |
                                      (std::uint32_t)std::max(u, w);
            auto it = edge_of.find(key);
            int id;
            if (it == edge_of.end()) {
                id = (int)G.edges.size();
                G.edges.push_back({u, w, {}});
                edge_of.emplace(key, id);
            } else {
                id = it->second;
            }
            auto& src = G.edges[id].src;
            if (src.empty() || src.back() != i) src.push_back(i);
            G.chain[i].push_back({id, G.edges[id].a == u});
        }
    }
    return G;
}

// Signed area of a node ring (CCW positive).
inline double ring_area(const std::vector<V2>& nodes, const std::vector<int>& ring) {
    double a2 = 0.0;
    for (size_t k = 0; k < ring.size(); ++k) {
        const V2& p = nodes[ring[k]];
        const V2& q = nodes[ring[(k + 1) % ring.size()]];
        a2 += p.x * q.y - q.x * p.y;
    }
    return 0.5 * a2;
}

// Even-odd containment of (x, y) in a node ring.
inline bool ring_contains(const std::vector<V2>& nodes, const std::vector<int>& ring,
                          double x, double y) {
    bool in = false;
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
        const V2& a = nodes[ring[i]];
        const V2& b = nodes[ring[j]];
        if (((a.y > y) != (b.y > y)) && (x < (b.x - a.x) * (y - a.y) / (b.y - a.y) + a.x))
            in = !in;
    }
    return in;
}

// Split a closed node walk into simple loops at every node it revisits. On a planar graph a
// walk with no repeated node cannot cross itself, so each loop is a simple polygon and has a
// well-defined orientation -- a pinched layer (two lobes touching at a point) comes apart into
// its two lobes, a bow tie into its two opposite-handed halves. Loops of fewer than three
// nodes (a spike walked out and back) enclose nothing and are dropped.
inline std::vector<std::vector<int>> simple_loops(const std::vector<int>& walk) {
    std::vector<std::vector<int>> loops;
    std::vector<int> stack;
    std::unordered_map<int, size_t> pos;
    auto close_at = [&](int v) {
        const size_t p = pos[v];
        std::vector<int> loop(stack.begin() + (long)p, stack.end());
        for (size_t k = p + 1; k < stack.size(); ++k) pos.erase(stack[k]);
        stack.resize(p + 1);
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    };
    for (int v : walk) {
        if (!stack.empty() && stack.back() == v) continue;
        if (pos.count(v)) { close_at(v); continue; }
        pos[v] = stack.size();
        stack.push_back(v);
    }
    if (stack.size() >= 3) loops.push_back(stack);
    return loops;
}

// A bounded face of the graph: its outer boundary (CCW) and the boundaries of the holes in it
// (CW), as node rings.
struct Face {
    std::vector<int> outer;
    std::vector<std::vector<int>> holes;
    double area = 0.0;   // outer area minus the holes
};

// Bounded faces of the graph restricted to the edges with keep[edge] (all edges when empty).
// Half-edge traversal (the next half-edge is the clockwise neighbour of the twin): CCW cycles
// are outer boundaries of bounded faces, CW cycles are holes (or the unbounded face). Each hole
// belongs to the smallest outer boundary that contains it; a hole contained in none bounds the
// unbounded face and is not returned.
inline std::vector<Face> faces(const PlanarGraph& G, const std::vector<char>& keep = {}) {
    const int nn = (int)G.nodes.size();
    std::vector<std::vector<int>> adj(nn);
    for (int k = 0; k < (int)G.edges.size(); ++k) {
        if (!keep.empty() && !keep[k]) continue;
        adj[G.edges[k].a].push_back(G.edges[k].b);
        adj[G.edges[k].b].push_back(G.edges[k].a);
    }
    for (int v = 0; v < nn; ++v) {
        auto ang = [&](int w) { return std::atan2(G.nodes[w].y - G.nodes[v].y, G.nodes[w].x - G.nodes[v].x); };
        std::sort(adj[v].begin(), adj[v].end(), [&](int a, int b) { return ang(a) < ang(b); });
    }
    auto hkey = [](int u, int v) { return ((std::uint64_t)(std::uint32_t)u << 32) | (std::uint32_t)v; };
    std::unordered_map<std::uint64_t, char> seen;
    std::vector<std::vector<int>> ccw, cw;
    for (int u = 0; u < nn; ++u)
        for (int v : adj[u]) {
            if (seen[hkey(u, v)]) continue;
            std::vector<int> cyc;
            int cu = u, cv = v;
            for (size_t guard = 0; guard <= 2 * G.edges.size() + 2; ++guard) {
                seen[hkey(cu, cv)] = 1;
                cyc.push_back(cu);
                const auto& out = adj[cv];
                const int idx = (int)(std::find(out.begin(), out.end(), cu) - out.begin());
                const int nx = out[(idx - 1 + (int)out.size()) % (int)out.size()];
                cu = cv; cv = nx;
                if (cu == u && cv == v) break;
            }
            const double a = ring_area(G.nodes, cyc);
            if (a > 0.0) ccw.push_back(std::move(cyc));
            else cw.push_back(std::move(cyc));
        }
    std::vector<Face> out(ccw.size());
    std::vector<double> outer_area(ccw.size());
    for (size_t f = 0; f < ccw.size(); ++f) {
        out[f].outer = ccw[f];
        out[f].area = outer_area[f] = ring_area(G.nodes, ccw[f]);
    }
    for (auto& h : cw) {
        // A CW cycle never shares a node with the face that encloses it (touching would have
        // joined the two into one cycle), but it shares every node with the faces it wraps --
        // so an outer ring that shares a node is not a candidate, and any node of the hole
        // then tests containment in the others away from their boundaries.
        int best = -1; double best_area = std::numeric_limits<double>::infinity();
        for (size_t f = 0; f < ccw.size(); ++f) {
            bool shares = false;
            for (int v : h) if (std::find(ccw[f].begin(), ccw[f].end(), v) != ccw[f].end()) { shares = true; break; }
            if (shares) continue;
            const V2& p = G.nodes[h[0]];
            if (ring_contains(G.nodes, ccw[f], p.x, p.y) && outer_area[f] < best_area) {
                best = (int)f; best_area = outer_area[f];
            }
        }
        if (best >= 0) { out[best].area += ring_area(G.nodes, h); out[best].holes.push_back(std::move(h)); }
    }
    return out;
}

}  // namespace katai::geometry
