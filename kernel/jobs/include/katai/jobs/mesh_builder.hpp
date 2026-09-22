#pragma once
// Build an FE mesh from a project's soil polygons (layer 2, katai/jobs; historical
// name build_mesh.hpp). The model-to-mesh half of the driver seam.
//
// Pipeline: the polygon edges and the lines the mesh must follow are noded into one
// planar graph at a tolerance scaled to the model (katai/geometry/planar_graph.hpp:
// near-coincident points merge, T-junctions and overlaps become shared vertices, a
// shared edge is one segment) -> each piece is kept when the owners on its two sides
// differ (owner = the last polygon containing the point, as everywhere else), so
// overlapping polygons and a lens inside a layer are regions, not holes -> Ruppert
// quality mesh (constrained Delaunay + refinement) -> triangles outside every soil
// polygon are dropped -> used vertices are compacted (no orphan/singular nodes) ->
// the linear triangulation is promoted to a tri6 / tri15 FE mesh, and every
// element takes the material of the last polygon containing its centroid.
//
// Pure geometry on katai::model data; no GUI dependency, so it is unit-testable.

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

namespace katai::app {

// Even-odd point-in-polygon test (world coordinates).
inline bool point_in_polygon(double x, double y, const model::SoilPolygon& P) {
    const size_t n = P.x.size(); if (n < 3) return false;
    bool in = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++)
        if (((P.y[i] > y) != (P.y[j] > y)) &&
            (x < (P.x[j] - P.x[i]) * (y - P.y[i]) / (P.y[j] - P.y[i]) + P.x[i]))
            in = !in;
    return in;
}

struct MeshResult {
    katai::mesh::Mesh mesh;
    bool ok = false;
    std::string message;
    // Did the refinement reach the quality bound it was given, or did its step cap stop it first?
    // Ruppert terminates for the angles asked here, so this is true on every ordinary model; it
    // exists because the alternative was a mesh that had given up being indistinguishable from one
    // that had not (see katai::mesh::Triangulation::quality_met for the argument).
    bool quality_met = true;
    double min_angle_asked = 0.0;   // the bound the refinement was given [deg]
    // Elements below that bound only because they sit in a corner of the geometry narrower than
    // it (Triangulation::corner_elements) -- a pinched-out layer, a sharp toe. Not a defect of
    // the mesh, but reported, because a corner can also be narrow by accident.
    int corner_elements = 0;
};

// Local mesh-density options (coarseness factors; docs/references/mesh-sizing.md).
// The target edge length h0 = sqrt(2 * max_area) is modulated by per-object coarseness factors:
//   - inside a polygon with factor f      -> h_region = h0 * f (last polygon wins, like material);
//   - near a structural line / load with factor f -> a line/point SOURCE of size h0 * f
//     (times auto_factor when auto_refine is on -- automatic refinement around
//     structures and loads), growing away from the source with the Lipschitz grading slope
//     h(d) = h_src + grading * d (bounded size transition; Shewchuk / Persson practice).
// h(x) = min(h_region(x), min over sources). Factors are clamped to [1/16, 4], so the field has
// a positive lower bound (Ruppert termination stays guaranteed).
struct MeshOptions {
    double min_angle_deg = 20.0;
    bool auto_refine = false;    // refine around structures + loads automatically (GUI default ON)
    double auto_factor = 0.5;    // automatic refinement factor for structures/loads
    double grading = 0.5;        // size growth per unit distance from a source
};

// The CONNECTION POINT of an embedded beam (docs/k2d-format.md, `conn`): for a pile it is the
// point on the embedded beam with the highest y-coordinate in the model, and for an exactly
// horizontal pile (a rare case) the end with the lowest x-coordinate. Returns true and writes
// (cx, cy) for an embedded beam.
//
// The mesher carries this ONE point as a vertex so a hinged connection is an exact degree-of-
// freedom identity rather than an interpolation; the shaft stays mesh-nonconforming. Both the
// mesher and the driver must agree on it to the last bit, which is why there is one definition.
inline bool embedded_connection_point(const model::StructElement& s, double& cx, double& cy) {
    if (s.kind != model::StructKind::EmbeddedBeam) return false;
    const bool second = (s.y2 != s.y1) ? (s.y2 > s.y1) : (s.x2 < s.x1);
    cx = second ? s.x2 : s.x1;
    cy = second ? s.y2 : s.y1;
    return true;
}

// max_area  : target maximum triangle area [m^2] (> 0).
// order     : 6 (tri6) or 15 (tri15).
// min_angle : Ruppert quality bound; <= ~20.7 deg guarantees termination.
// Definition in kernel/jobs/src/mesh_builder.cpp (section 5.2).
MeshResult mesh_from_project(const model::Project& pr, double max_area,
                             int order, const MeshOptions& opt);

inline MeshResult mesh_from_project(const model::Project& pr, double max_area,
                                    int order = 6, double min_angle_deg = 20.0) {
    MeshOptions opt; opt.min_angle_deg = min_angle_deg;
    return mesh_from_project(pr, max_area, order, opt);
}

// Mesh from the project's OWN MeshSettings -- the .k2d-determined discretization.
// Every front end that meshes "the file" must come through here, so a checked-in
// project reproduces the same mesh everywhere.
inline MeshResult mesh_from_project(const model::Project& pr) {
    MeshOptions opt;
    opt.auto_refine = pr.mesh.auto_refine;
    return mesh_from_project(pr, 0.5 * pr.mesh.elem_size * pr.mesh.elem_size, pr.mesh.order, opt);
}

} // namespace katai::app
