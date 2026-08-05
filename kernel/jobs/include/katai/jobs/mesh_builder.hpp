#pragma once
// Build an FE mesh from a project's soil polygons (layer 2, katai/jobs; historical
// name build_mesh.hpp). The model-to-mesh half of the driver seam.
//
// Pipeline: the polygon edges form a planar straight-line graph (PSLG, shared
// vertices deduplicated so adjacent regions stay conforming) -> Ruppert quality
// mesh (constrained Delaunay + refinement) -> triangles whose centroid falls
// outside every soil polygon are dropped (handles concavity, holes and the
// convex-hull fill) -> used vertices are compacted (no orphan/singular nodes) ->
// the linear triangulation is promoted to a tri6 / tri15 FE mesh, and every
// element takes the material of the polygon containing its centroid.
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
};

// Local mesh-density options (PLAXIS coarseness semantics; docs/references/mesh-sizing.md).
// The target edge length h0 = sqrt(2 * max_area) is modulated by per-object coarseness factors:
//   - inside a polygon with factor f      -> h_region = h0 * f (last polygon wins, like material);
//   - near a structural line / load with factor f -> a line/point SOURCE of size h0 * f
//     (times auto_factor when auto_refine is on -- PLAXIS-style automatic refinement around
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
