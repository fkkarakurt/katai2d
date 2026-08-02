#pragma once
// Incremental Delaunay triangulation of a 2D point set (P1.4b).
//
// Bowyer-Watson insertion on a triangle-based data structure with O(1) neighbour
// adjacency. Robustness rests entirely on the adaptive orient2d / incircle
// predicates (P1.4a), so the construction never produces inverted elements or
// loops on near-degenerate (e.g. cocircular grid) input. Points are inserted in
// Morton (Z-order) sequence so that point location -- a walk from the previous
// triangle -- stays short, giving near-linear behaviour in practice. The same
// incremental kernel will host constrained-segment recovery (P1.4c) and Ruppert
// refinement (P1.4d).

#include <array>
#include <functional>
#include <vector>

namespace katai::mesh {

// Position-dependent sizing field: the MAXIMUM triangle area [m^2] allowed at point
// (x, y). Returning <= 0 means no area bound at that point (angle quality only). A
// constant max_area is its special case. Local refinement (the PLAXIS coarseness factor /
// refine) is done through this field; the field must have a POSITIVE lower bound
// (Ruppert's termination guarantee).
using SizeField = std::function<double(double x, double y)>;

// A triangulation. Triangle vertices are CCW indices into (x, y). Input points
// keep their original indices 0..n_input-1; any Steiner points added by
// refinement are appended after them. Triangles incident to the auxiliary
// super-triangle are removed, so the result covers the (constrained) domain.
struct Triangulation {
    std::vector<double> x, y;                    // vertex coordinates
    std::vector<std::array<int, 3>> triangles;   // CCW vertex-index triples
    int point_count = 0;                         // number of input points
};

// Delaunay-triangulate the points (px, py). Inputs must be distinct; the result
// is the (unique, in general position) Delaunay triangulation of their convex hull.
Triangulation delaunay_triangulate(const std::vector<double>& px,
                                   const std::vector<double>& py);

// Constrained Delaunay triangulation of a planar straight-line graph: every
// segment (a pair of indices into the points) is forced to appear as an edge via
// flip-based recovery (Anglada), after which the Delaunay property is restored on
// all non-constrained edges. Segments must not cross one another and no vertex
// may lie in the interior of a segment.
Triangulation constrained_delaunay(
    const std::vector<double>& px, const std::vector<double>& py,
    const std::vector<std::array<int, 2>>& segments);

// Quality mesh by Ruppert's Delaunay refinement of the constrained triangulation:
// encroached subsegments are split at their midpoints and skinny triangles are
// eliminated by inserting circumcentres, until no triangle has a minimum angle
// below min_angle_deg (use <= ~20.7 deg for guaranteed termination) and, when
// max_area > 0, no triangle exceeds that area. Steiner points are appended to the
// returned (x, y); input point indices are preserved.
//
// `segments` are ALL the constraints the mesh must conform to (recovered as edges and
// subdivided for quality). `outline` is the subset that forms the domain BOUNDARY, used
// for the inside/outside (in_domain) test that drops exterior triangles -- it must contain
// ONLY the closed domain polygon(s), never internal lines (plates/geogrids); an internal
// segment treated as a boundary would flip the ray-casting parity and delete the mesh on
// one side of it. When `outline` is empty, every segment is treated as boundary (the
// original behaviour, correct only when there are no internal constraint lines).
Triangulation quality_mesh(
    const std::vector<double>& px, const std::vector<double>& py,
    const std::vector<std::array<int, 2>>& segments, double min_angle_deg,
    double max_area, const std::vector<std::array<int, 2>>& outline = {});

// Quality mesh with a sizing field: instead of a constant max_area, a
// position-dependent area bound evaluated at each triangle's CENTROID (local refinement —
// region/line/point coarseness). Everything else (Ruppert min-angle, encroachment,
// boundary) is identical to the constant version.
Triangulation quality_mesh(
    const std::vector<double>& px, const std::vector<double>& py,
    const std::vector<std::array<int, 2>>& segments, double min_angle_deg,
    const SizeField& max_area_at, const std::vector<std::array<int, 2>>& outline = {});

} // namespace katai::mesh
