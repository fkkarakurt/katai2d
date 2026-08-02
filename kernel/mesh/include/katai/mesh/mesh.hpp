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

} // namespace katai::mesh
