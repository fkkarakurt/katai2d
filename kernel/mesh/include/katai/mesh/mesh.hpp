#pragma once
// Mesh (minimal) — Phase 0. Splits a rectangle into structured tri6 triangles.
// The full Delaunay/Ruppert mesher was deferred to Phase 1; this structured mesh
// suffices to exercise the solver pipeline (assembly→solve→post).
//
// Data-oriented SoA (see docs/ARCHITECTURE): contiguous arrays, flat connectivity.

#include <vector>

#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/delaunay.hpp>

namespace katai::mesh {

// Triangle mesh, structure-of-arrays. Flat connectivity: element e's nodes are
// connectivity[e*nodes_per_element .. (e+1)*nodes_per_element).
struct Mesh {
    std::vector<double> x;                 // node x coordinates (node_count)
    std::vector<double> y;                 // node y coordinates (node_count)
    std::vector<int> connectivity;         // element_count * nodes_per_element
    std::vector<int> element_material;     // material id per element (element_count)
    int nodes_per_element = 6;
    int node_count = 0;
    int element_count = 0;

    // Boundary node sets (for Phase 0 BCs: base clamped, side rollers).
    std::vector<int> bottom_nodes;
    std::vector<int> top_nodes;
    std::vector<int> left_nodes;
    std::vector<int> right_nodes;

    // Generic domain-boundary nodes (corner + mid-edge), for meshes produced by
    // the unstructured mesher where the axis-aligned named sets do not apply.
    std::vector<int> boundary_nodes;

    // The global index of element e's n-th node.
    int node_of(int element, int local) const {
        return connectivity[element * nodes_per_element + local];
    }
};

// Splits the rectangular domain into nx*ny cells, each cell into 2 tri6 triangles.
// The generated node grid is (2*nx+1) x (2*ny+1) (including the quadratic mid-edge
// nodes). All elements are positive-area (CCW) with domain.material_id.
Mesh generate_structured_tri6(const geometry::RectangularDomain& domain,
                              int nx, int ny);

// Structured 15-node (quartic) mesh: nx*ny cells, 2 tri15 per cell. Uses a
// (4nx+1) x (4ny+1) fine grid -- every tri15 node (corners, the 3 per edge, and
// the 3 interior) lands on a quarter-grid point, so node sharing is automatic.
// Ordered bottom/top/left/right boundary node sets are provided (incl. edge nodes)
// for boundary conditions and surface tractions.
Mesh generate_structured_tri15(const geometry::RectangularDomain& domain,
                               int nx, int ny);

// Build a quadratic (6-node) FEM mesh from a linear triangulation produced by the
// mesher: the triangulation vertices become corner nodes and one shared mid-edge
// node is created per unique edge. Element connectivity follows the tri6 local
// ordering (corners 0,1,2; mid-edge nodes 3,4,5 on edges 0-1, 1-2, 2-0). All
// elements inherit material_id. boundary_nodes collects the nodes on edges that
// bound a single element (the domain boundary).
Mesh tri6_from_triangulation(const Triangulation& triangulation,
                             int material_id = 0);

// Build a quartic (15-node) FEM mesh from a linear triangulation. Each triangle
// gets its 3 corner nodes, 3 shared nodes per edge (at 1/4, 2/4, 3/4 along the
// edge), and 3 element-interior nodes, in the tri15 local ordering (see
// elements/tri15.hpp). Edge nodes are shared by adjacent triangles via a canonical
// orientation; interior nodes are private. boundary_nodes collects nodes on edges
// bounding a single element (corners + that edge's three nodes).
Mesh tri15_from_triangulation(const Triangulation& triangulation,
                              int material_id = 0);

// Uniform ("red", or regular) refinement: every triangle is split into four by the
// midpoints of its edges (Bank, Sherman and Weiser 1983). Works on tri6 and tri15,
// and returns a mesh of the same element order, with each child inheriting its
// parent's material.
//
// The three properties that make this the right instrument for a convergence study,
// and that a fresh mesh at half the element size does NOT have:
//
//   NESTED   -- every CORNER and EDGE node of the coarse mesh is a node of the refined
//               mesh, at the same coordinates. Successive solutions are then comparable
//               point by point, and their differences are discretisation rather than a
//               different mesh's idea of where to put nodes. (The exception, measured
//               and stated rather than discovered later: a tri15 element's three
//               INTERIOR nodes sit at barycentric (2,1,1)/4 and its permutations, which
//               are not lattice points of any child, so they do not survive. They are
//               private to their element and shared with nothing, so no continuity or
//               comparability is lost -- but "every node survives" would be false, and
//               tests/test_mesh_refine_uniform.cpp asserts the true statement.)
//   SIMILAR  -- the four children are similar to the parent (midpoint triangle), so
//               every angle is preserved EXACTLY. The refined family is shape-regular
//               by construction, which is the hypothesis the O(h^p) convergence
//               theory is proved under (Ciarlet 1978), and it removes mesh quality
//               as a variable between the runs.
//   EXACT r  -- h halves exactly, so the refinement factor is 2, comfortably above
//               the 1.3 that Celik et al. (2008) require of a grid-refinement study.
//
// Rebuilding a mesh from scratch at a smaller target size gives none of the three:
// the meshes are not nested, the angles differ, and the achieved ratio is whatever
// the mesher happened to produce. Measured on this tree, that difference is what put
// KV-NUM-005's first triplets outside the asymptotic range (docs/validation/
// numerical-uncertainty.md).
//
// Geometry is preserved exactly: the boundary of the refined mesh is the boundary of
// the coarse one, because a midpoint of a straight edge lies on that edge. A CURVED
// domain would be a different matter -- refinement here does not know the true
// geometry and cannot pull new boundary nodes onto it -- and this program's domains
// are straight-sided polygons, so the question does not arise.
Mesh refine_uniform(const Mesh& coarse);

} // namespace katai::mesh
