#pragma once
// Region editing on drawn rings: close a new region against the boundaries that already exist,
// and cut a region in pieces with a line. Both are faces of one planar graph
// (katai/geometry/planar_graph.hpp), so they meet existing geometry exactly where the mesher
// will see it meet, at the same tolerance.
//
//   close_against(rings, line)  -- `line` is an open polyline whose ends lie on existing
//       boundaries. The regions it encloses together with those boundaries, and that no
//       existing ring already covers, are returned. This is how a region is drawn against its
//       neighbours without tracing the shared boundary a second time: start on a neighbour,
//       click the free vertices, finish on a neighbour.
//   split_ring(ring, line)      -- the pieces `line` cuts `ring` into (fewer than two: no cut).
//       A closed line drawn inside the ring cuts out the piece it encloses.
//
// A returned ring is CCW. A face with holes comes back as ONE ring joined to each hole by a
// zero-width bridge (a keyhole): under the even-odd rule every consumer here applies
// (point_in_polygon, the mesher's loop split) the bridge encloses nothing and the hole stays out.

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <katai/geometry/planar_graph.hpp>

namespace katai::geometry {

using Ring = std::vector<V2>;

namespace detail {

inline bool ring_contains_pt(const Ring& r, double x, double y) {
    bool in = false;
    for (size_t i = 0, j = r.size() - 1; i < r.size(); j = i++)
        if (((r[i].y > y) != (r[j].y > y)) &&
            (x < (r[j].x - r[i].x) * (y - r[i].y) / (r[j].y - r[i].y) + r[i].x))
            in = !in;
    return in;
}

// Edges that bound something: drop, repeatedly, every edge with an end of degree one (a line
// that stops short of a boundary, or runs past it).
inline std::vector<char> bounding_edges(const PlanarGraph& G) {
    std::vector<char> keep(G.edges.size(), 1);
    std::vector<int> deg(G.nodes.size(), 0);
    for (const auto& e : G.edges) { ++deg[e.a]; ++deg[e.b]; }
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t k = 0; k < G.edges.size(); ++k) {
            if (!keep[k]) continue;
            if (deg[G.edges[k].a] <= 1 || deg[G.edges[k].b] <= 1) {
                keep[k] = 0; --deg[G.edges[k].a]; --deg[G.edges[k].b]; changed = true;
            }
        }
    }
    return keep;
}

inline Ring ring_points(const PlanarGraph& G, const std::vector<int>& ids) {
    Ring r; r.reserve(ids.size());
    for (int id : ids) r.push_back(G.nodes[id]);
    return r;
}

// A point strictly inside the face: just inside some edge of the outer ring, checked against
// the outer ring and every hole.
inline bool face_point(const PlanarGraph& G, const Face& f, double& px, double& py) {
    const Ring outer = ring_points(G, f.outer);
    std::vector<Ring> holes;
    for (const auto& h : f.holes) holes.push_back(ring_points(G, h));
    for (double frac : {0.25, 0.05, 0.01, 1e-3}) {
        for (size_t k = 0; k < outer.size(); ++k) {
            const V2& a = outer[k];
            const V2& b = outer[(k + 1) % outer.size()];
            const double dx = b.x - a.x, dy = b.y - a.y, L = std::hypot(dx, dy);
            if (L <= 0.0) continue;
            const double d = std::max(frac * L, 4.0 * G.eps);
            const double x = 0.5 * (a.x + b.x) - dy / L * d, y = 0.5 * (a.y + b.y) + dx / L * d;   // left of a->b
            if (!ring_contains_pt(outer, x, y)) continue;
            bool in_hole = false;
            for (const auto& h : holes) if (ring_contains_pt(h, x, y)) { in_hole = true; break; }
            if (in_hole) continue;
            px = x; py = y;
            return true;
        }
    }
    return false;
}

inline bool proper_cross(const V2& p, const V2& p2, const V2& q, const V2& q2) {
    auto orient = [](const V2& a, const V2& b, const V2& c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    };
    const double d1 = orient(q, q2, p), d2 = orient(q, q2, p2);
    const double d3 = orient(p, p2, q), d4 = orient(p, p2, q2);
    return ((d1 > 0) != (d2 > 0)) && d1 != 0 && d2 != 0 && ((d3 > 0) != (d4 > 0)) && d3 != 0 && d4 != 0;
}

// Join each hole to the outer ring by the shortest bridge that crosses nothing.
inline Ring keyhole(const PlanarGraph& G, const Face& f) {
    Ring outer = ring_points(G, f.outer);
    std::vector<Ring> holes;
    for (const auto& h : f.holes) holes.push_back(ring_points(G, h));
    for (size_t hi = 0; hi < holes.size(); ++hi) {
        const Ring& h = holes[hi];
        auto blocked = [&](const V2& a, const V2& b) {
            auto test = [&](const Ring& r) {
                for (size_t k = 0; k < r.size(); ++k)
                    if (proper_cross(a, b, r[k], r[(k + 1) % r.size()])) return true;
                return false;
            };
            if (test(outer)) return true;
            for (size_t j = hi; j < holes.size(); ++j) if (test(holes[j])) return true;
            return false;
        };
        size_t bo = 0, bh = 0; double best = -1.0;
        for (size_t i = 0; i < outer.size(); ++i)
            for (size_t j = 0; j < h.size(); ++j) {
                const double d = std::hypot(outer[i].x - h[j].x, outer[i].y - h[j].y);
                if (best >= 0.0 && d >= best) continue;
                if (blocked(outer[i], h[j])) continue;
                best = d; bo = i; bh = j;
            }
        if (best < 0.0) continue;   // no clear bridge (cannot happen on a planar face); leave the hole
        Ring joined(outer.begin(), outer.begin() + (long)bo + 1);
        for (size_t k = 0; k <= h.size(); ++k) joined.push_back(h[(bh + k) % h.size()]);
        joined.insert(joined.end(), outer.begin() + (long)bo, outer.end());
        outer = std::move(joined);
    }
    return outer;
}

inline void add_ring(std::vector<std::array<V2, 2>>& segs, const Ring& r) {
    for (size_t k = 0; k < r.size(); ++k) segs.push_back({r[k], r[(k + 1) % r.size()]});
}

inline double extent_tolerance(const std::vector<std::array<V2, 2>>& segs) {
    if (segs.empty()) return 1e-9;
    double x0 = segs[0][0].x, x1 = x0, y0 = segs[0][0].y, y1 = y0;
    for (const auto& s : segs)
        for (const V2& p : s) {
            x0 = std::min(x0, p.x); x1 = std::max(x1, p.x);
            y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
        }
    return geometry_tolerance(x0, y0, x1, y1);
}

}  // namespace detail

// eps <= 0: the model tolerance of the input's extent (geometry_tolerance).
inline std::vector<Ring> close_against(const std::vector<Ring>& rings, const std::vector<V2>& line,
                                       double eps = 0.0) {
    if (line.size() < 2) return {};
    std::vector<std::array<V2, 2>> segs;
    for (const auto& r : rings) if (r.size() >= 3) detail::add_ring(segs, r);
    const int first_line = (int)segs.size();
    for (size_t k = 0; k + 1 < line.size(); ++k) segs.push_back({line[k], line[k + 1]});
    const PlanarGraph G = planar_graph(segs, eps > 0.0 ? eps : detail::extent_tolerance(segs));
    const std::vector<char> keep = detail::bounding_edges(G);
    std::vector<char> from_line(G.edges.size(), 0);
    for (size_t k = 0; k < G.edges.size(); ++k)
        for (int s : G.edges[k].src) if (s >= first_line) from_line[k] = 1;
    // Node pair -> edge, to ask whether a face boundary runs along the new line.
    auto touches_line = [&](const std::vector<int>& cyc) {
        for (size_t k = 0; k < cyc.size(); ++k) {
            const int u = cyc[k], w = cyc[(k + 1) % cyc.size()];
            for (size_t e = 0; e < G.edges.size(); ++e)
                if (from_line[e] && keep[e] &&
                    ((G.edges[e].a == u && G.edges[e].b == w) || (G.edges[e].a == w && G.edges[e].b == u)))
                    return true;
        }
        return false;
    };
    std::vector<Ring> out;
    for (const Face& f : faces(G, keep)) {
        bool adj = touches_line(f.outer);
        for (const auto& h : f.holes) adj = adj || touches_line(h);
        if (!adj) continue;
        double px, py;
        if (!detail::face_point(G, f, px, py)) continue;
        bool covered = false;
        for (const auto& r : rings) if (r.size() >= 3 && detail::ring_contains_pt(r, px, py)) { covered = true; break; }
        if (covered) continue;
        out.push_back(detail::keyhole(G, f));
    }
    return out;
}

inline std::vector<Ring> split_ring(const Ring& ring, const std::vector<V2>& line, double eps = 0.0) {
    if (ring.size() < 3 || line.size() < 2) return {};
    std::vector<std::array<V2, 2>> segs;
    detail::add_ring(segs, ring);
    for (size_t k = 0; k + 1 < line.size(); ++k) segs.push_back({line[k], line[k + 1]});
    const PlanarGraph G = planar_graph(segs, eps > 0.0 ? eps : detail::extent_tolerance(segs));
    std::vector<Ring> out;
    for (const Face& f : faces(G, detail::bounding_edges(G))) {
        double px, py;
        if (!detail::face_point(G, f, px, py)) continue;
        if (!detail::ring_contains_pt(ring, px, py)) continue;
        out.push_back(detail::keyhole(G, f));
    }
    return out;
}

}  // namespace katai::geometry
