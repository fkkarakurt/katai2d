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
// constant max_area is its special case. Local refinement (a per-region coarseness factor /
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
    // DID THE REFINEMENT REACH ITS OWN QUALITY BOUND? Ruppert's algorithm terminates for the
    // angles this tree asks for, and when it does every in-domain triangle satisfies the bound --
    // which is why quality here is structural rather than measured per element. But the loop
    // carries a step cap as a final safety valve, and a valve that opens in silence is not a
    // safety valve: until this field existed, a mesh that met 20 degrees and a mesh that gave up
    // trying were the same object to everything downstream, and mesh quality is not a detail the
    // answer is indifferent to.
    //
    // Same argument, and the same shape, as the constitutive integrator running out of substeps
    // (K2D-A012): the guard is right, and the silence was the defect.
    //
    // One exception is not a failure of the refinement and is counted apart: a triangle in the
    // corner of an INPUT angle narrower than the bound. No mesh of that geometry can put an
    // angle wider than the corner's own there, and trying (Miller, Pav & Walkington; see
    // in_small_angle_corner) only makes smaller copies of the same triangle. quality_met stays
    // true for them; corner_elements says how many there are, so a corner drawn narrow by
    // accident (a vertex a hair off a straight line) is still visible.
    bool quality_met = true;
    int refinement_steps = 0;                    // how many the refinement actually took
    int corner_elements = 0;                     // below the bound only because of a narrow input corner
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
