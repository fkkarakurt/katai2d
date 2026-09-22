#pragma once
// Soil-region edits the drawing tools make, on the project model (header-only).
//
//   close_region(project, line)    -- the new regions an open line closes against the regions
//                                     that exist (katai::geometry::close_against), appended.
//   split_polygons(project, line)  -- every region the line cuts is replaced by its pieces.
//
// A piece is its parent in everything but shape: material, coarseness and, in every phase,
// activation (katai::model::insert_polygons keeps the phase flags on the right polygons). The
// per-edge conditions (deformation BC, flow BC, head, flux) follow the edges: a piece edge lying
// on an edge of the parent keeps that edge's conditions, and a cut edge -- which is interior to
// the old region -- starts free, as the interior it was.

#include <cmath>
#include <string>
#include <vector>

#include <katai/geometry/region_ops.hpp>
#include <katai/model/project.hpp>

namespace katai::model {

namespace detail {

inline geometry::Ring ring_of(const SoilPolygon& P) {
    geometry::Ring r;
    for (size_t k = 0; k < P.x.size() && k < P.y.size(); ++k) r.push_back({P.x[k], P.y[k]});
    return r;
}

inline double model_tolerance(const Project& p, const std::vector<geometry::V2>& line) {
    double x0 = 0, x1 = 0, y0 = 0, y1 = 0; bool any = false;
    auto grow = [&](double x, double y) {
        if (!any) { x0 = x1 = x; y0 = y1 = y; any = true; return; }
        x0 = std::min(x0, x); x1 = std::max(x1, x); y0 = std::min(y0, y); y1 = std::max(y1, y);
    };
    for (const auto& P : p.polygons) for (size_t k = 0; k < P.x.size(); ++k) grow(P.x[k], P.y[k]);
    for (const auto& v : line) grow(v.x, v.y);
    return geometry::geometry_tolerance(x0, y0, x1, y1);
}

// A polygon with shape `r` and every other property of `parent` (none when parent is null).
inline SoilPolygon piece(const geometry::Ring& r, const SoilPolygon* parent, double eps) {
    SoilPolygon P;
    if (parent) { P.name = parent->name; P.material = parent->material; P.coarseness = parent->coarseness; }
    const size_t n = r.size();
    for (const auto& v : r) { P.x.push_back(v.x); P.y.push_back(v.y); }
    P.edge_bc.assign(n, 0);
    if (!parent) return P;
    const bool flow = parent->edge_flow.size() == parent->x.size();
    const bool head = parent->edge_head.size() == parent->x.size();
    const bool flux = parent->edge_flux.size() == parent->x.size();
    if (flow) P.edge_flow.assign(n, 0);
    if (head) P.edge_head.assign(n, 0.0);
    if (flux) P.edge_flux.assign(n, 0.0);
    const size_t m = parent->x.size();
    auto on_edge = [&](const geometry::V2& q, size_t k) {
        const geometry::V2 a{parent->x[k], parent->y[k]}, b{parent->x[(k + 1) % m], parent->y[(k + 1) % m]};
        double t;
        return geometry::detail::dist_point_segment(q, a, b, t) < eps;
    };
    for (size_t i = 0; i < n; ++i) {
        const geometry::V2 &a = r[i], &b = r[(i + 1) % n];
        for (size_t k = 0; k < m; ++k) {
            if (!on_edge(a, k) || !on_edge(b, k)) continue;
            if (k < parent->edge_bc.size()) P.edge_bc[i] = parent->edge_bc[k];
            if (flow) P.edge_flow[i] = parent->edge_flow[k];
            if (head) P.edge_head[i] = parent->edge_head[k];
            if (flux) P.edge_flux[i] = parent->edge_flux[k];
            break;
        }
    }
    return P;
}

}  // namespace detail

// Returns how many regions were added (0: the line closes nothing new -- its ends do not both
// reach existing boundaries, or everything it encloses is already soil).
inline int close_region(Project& p, const std::vector<geometry::V2>& line) {
    std::vector<geometry::Ring> rings;
    for (const auto& P : p.polygons) rings.push_back(detail::ring_of(P));
    const double eps = detail::model_tolerance(p, line);
    const auto made = geometry::close_against(rings, line, eps);
    int added = 0;
    for (const auto& r : made) {
        SoilPolygon P = detail::piece(r, nullptr, eps);
        P.name = "Soil " + std::to_string(p.polygons.size() + 1);
        p.polygons.push_back(std::move(P));
        ++added;
    }
    return added;
}

// Returns how many regions were cut (each into two or more pieces).
inline int split_polygons(Project& p, const std::vector<geometry::V2>& line) {
    const double eps = detail::model_tolerance(p, line);
    int cut = 0;
    for (size_t i = p.polygons.size(); i-- > 0;) {   // back to front: insertions never move a polygon still to visit
        const auto pieces = geometry::split_ring(detail::ring_of(p.polygons[i]), line, eps);
        if (pieces.size() < 2) continue;
        const SoilPolygon parent = p.polygons[i];
        std::vector<SoilPolygon> rest;
        for (size_t k = 0; k < pieces.size(); ++k) {
            SoilPolygon P = detail::piece(pieces[k], &parent, eps);
            if (k > 0) P.name = parent.name + " (" + std::to_string(k + 1) + ")";
            if (k == 0) p.polygons[i] = std::move(P);
            else rest.push_back(std::move(P));
        }
        insert_polygons(p, i + 1, rest, (int)i);
        ++cut;
    }
    return cut;
}

}  // namespace katai::model
