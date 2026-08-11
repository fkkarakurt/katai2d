#pragma once
// Staged construction — element active/passive (P2.2). The core workflow of geotechnical
// FEM: element sets change between phases (soil goes passive on excavation, active on
// fill); equilibrium is re-established each phase from the carried initial stress + active
// set + active weight. The excavation mechanics is AUTOMATIC: when the solver is given
// (i) the previous phase's stress as committed pre-stress (solve_nonlinear initial_state),
// (ii) the active mask, (iii) active gravity, the first residual = f_active − ∫_active Bᵀσ0
// contains exactly the support loss of the removed soil (the released surface load) → the
// remaining mass comes to equilibrium (base heave). The phase displacement starts from zero
// (incremental).
//
// Math + sources: docs/references/initial-stress-k0.md §5; PLAXIS Reference/
// Scientific Manual (staged construction).

#include <algorithm>
#include <cmath>
#include <vector>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Mesh node SPLITTING for an embedded wall (barrier). A wall is a barrier separating the
// soil into two sides: in a continuous mesh the nodes along the wall line are shared, so
// the two sides stay tied (the wall does not separate). This helper PAIRS the nodes on the
// vertical line x≈x_wall in the range y_toe < y ≤ y_top: elements with x < x_wall (left)
// are re-wired to twin (new) nodes → the left and right soils remain connected only BELOW
// the TOE (y_toe) (continuous soil). The wall (plate + interfaces on both sides) attaches
// to these (right, left-twin) pairs. The mesh changes in place: node_count grows, the
// left-element connectivity is updated. x_wall must be a mesh column line (nodes aligned in
// a structured mesh). Returns: (right, left, y) seam pairs (not sorted top-down — the
// caller may sort by y). Must be called BEFORE the DofMap is built.
struct SeamPair { int right; int left; double y; };
inline std::vector<SeamPair> split_mesh_at_wall(mesh::Mesh& mesh, double x_wall,
                                                double y_toe, double y_top,
                                                double tol = 1e-6) {
    std::vector<int> dup(mesh.node_count, -1);  // original node → twin (left) node
    std::vector<SeamPair> pairs;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - x_wall) < tol && mesh.y[n] > y_toe + tol &&
            mesh.y[n] < y_top + tol) {
            dup[n] = mesh.node_count + static_cast<int>(pairs.size());
            pairs.push_back({n, dup[n], mesh.y[n]});
        }
    for (const auto& p : pairs) {  // add the twin nodes (same location)
        mesh.x.push_back(x_wall);
        mesh.y.push_back(p.y);
    }
    mesh.node_count += static_cast<int>(pairs.size());
    for (int e = 0; e < mesh.element_count; ++e) {  // re-wire the left elements to the twins
        double xc = 0.0;
        for (int k = 0; k < mesh.nodes_per_element; ++k) xc += mesh.x[mesh.node_of(e, k)];
        if (xc / mesh.nodes_per_element >= x_wall) continue;  // right side → unchanged
        for (int k = 0; k < mesh.nodes_per_element; ++k) {
            const int idx = e * mesh.nodes_per_element + k;
            if (dup[mesh.connectivity[idx]] >= 0)
                mesh.connectivity[idx] = dup[mesh.connectivity[idx]];
        }
    }
    return pairs;
}

// GENERAL (arbitrarily oriented) mesh node SPLITTING — the generalization of
// split_mesh_at_wall to any direction. Along the segment (ax,ay)→(bx,by), PAIRS the nodes
// sitting on mesh edges: elements on the NEGATIVE side of the segment (signed distance d<0)
// are re-wired to twin (new) nodes → the two sides remain connected only through the
// unsplit nodes, and an interface (Coulomb slip/separation) can be placed between them.
// `s` = the along-segment parameter (∈[0,L]); only nodes in (s_lo, s_hi) are split → the
// end nodes stay SHARED (instead of a free in-soil tip; like a wall toe). The segment must
// be a mesh edge chain (a structural/interface line is embedded in the mesh as a
// constraint → nodes aligned). The mesh changes in place (node_count grows). Returns: seam
// pairs SORTED by s (corner,mid,corner… tri6 / corner,q,m,q,corner… tri15). Must be called
// BEFORE the DofMap is built. (PLAXIS interface: a structural/in-soil slip surface.)
struct SegSeam { int orig; int dup; double x, y, s; };
inline std::vector<SegSeam> split_mesh_at_segment(mesh::Mesh& mesh, double ax, double ay,
                                                  double bx, double by, double s_lo, double s_hi,
                                                  double tol = 1e-6) {
    const double L = std::hypot(bx - ax, by - ay);
    if (L < 1e-12) return {};
    const double tx = (bx - ax) / L, ty = (by - ay) / L;  // unit tangent
    const double nx = ty, ny = -tx;                        // unit normal (right-handed)
    std::vector<int> dup(mesh.node_count, -1);
    std::vector<SegSeam> seam;
    for (int i = 0; i < mesh.node_count; ++i) {
        const double rx = mesh.x[i] - ax, ry = mesh.y[i] - ay;
        const double s = rx * tx + ry * ty;       // along-segment parameter
        const double d = rx * nx + ry * ny;       // signed perpendicular distance
        if (std::fabs(d) < tol && s > s_lo - tol && s < s_hi + tol) {
            dup[i] = mesh.node_count + static_cast<int>(seam.size());
            seam.push_back({i, dup[i], mesh.x[i], mesh.y[i], s});
        }
    }
    for (const auto& p : seam) { mesh.x.push_back(p.x); mesh.y.push_back(p.y); }
    mesh.node_count += static_cast<int>(seam.size());
    for (int e = 0; e < mesh.element_count; ++e) {  // re-wire the negative-side elements to the twins
        double cx = 0.0, cy = 0.0;
        for (int k = 0; k < mesh.nodes_per_element; ++k) { cx += mesh.x[mesh.node_of(e, k)]; cy += mesh.y[mesh.node_of(e, k)]; }
        cx /= mesh.nodes_per_element; cy /= mesh.nodes_per_element;
        if ((cx - ax) * nx + (cy - ay) * ny >= 0.0) continue;   // positive side → unchanged
        for (int k = 0; k < mesh.nodes_per_element; ++k) {
            const int idx = e * mesh.nodes_per_element + k;
            if (dup[mesh.connectivity[idx]] >= 0)
                mesh.connectivity[idx] = dup[mesh.connectivity[idx]];
        }
    }
    std::sort(seam.begin(), seam.end(), [](const SegSeam& a, const SegSeam& b) { return a.s < b.s; });
    return seam;
}

// Marks the nodes touched by at least one active element (active_element empty = all active).
inline std::vector<char> active_nodes(const mesh::Mesh& mesh,
                                      const std::vector<char>& active_element) {
    std::vector<char> na(mesh.node_count, 0);
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        for (int k = 0; k < mesh.nodes_per_element; ++k)
            na[mesh.node_of(e, k)] = 1;
    }
    return na;
}

// Fixes both components of nodes touching no active element (orphaned by excavation) →
// the reduced system stays non-singular. Must be called BEFORE DofMap::finalize().
//
// `carried` (optional, node_count) exempts nodes that a STRUCTURAL element already holds:
// such a node is not orphaned, and fixing it welds the structure to the outside world --
// silently, because the run then converges on a model in which that structure carries
// nothing. The exemption is for PLATES only. A plate is a Timoshenko beam: its axial,
// bending and shear stiffness together hold BOTH in-plane translations of every node of
// its chain, so releasing them leaves a solvable system. A geogrid or an anchor is
// axial-only and an embedded beam has its own DOFs; exempting those would leave the
// transverse direction singular, so they keep the fixity. Empty = nothing exempted (the
// old behaviour, bit-for-bit).
//
// The case this exists for is the manual's own: PLAXIS 2D Validation Manual §2.3 builds a
// beam by deactivating the soil cluster so that only the beams remain.
inline void fix_inactive_nodes(const mesh::Mesh& mesh,
                               const std::vector<char>& active_element,
                               DofMap& dofs,
                               const std::vector<char>& carried = {}) {
    const std::vector<char> na = active_nodes(mesh, active_element);
    for (int n = 0; n < mesh.node_count; ++n) {
        if (!na[n] && (carried.empty() || !carried[(size_t)n])) {
            dofs.fix_node_component(n, 0);
            dofs.fix_node_component(n, 1);
        }
    }
}

}  // namespace katai::core
