#pragma once
// Displacement boundary conditions on the domain boundary (Stage B3: extracted
// from the application driver). The engine's contract is a neutral list of
// boundary edges -- segment, unit outward normal, fixity -- and it owns every
// load-bearing convention: node-on-edge matching with its tie tolerance, the
// union semantics at corners, and the compliant-base horizontal release. The
// caller that owns a project schema walks its polygons once and hands the edge
// list over; the engine never sees the schema.

#include <algorithm>
#include <cmath>
#include <vector>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Fixity of one boundary edge. The engine's own vocabulary; a schema enum is
// mapped once, at the caller's seam.
enum class EdgeFixity { Free, NormallyFixed, HorizontallyFixed, VerticallyFixed, FullyFixed };

// One boundary edge: the segment (ax,ay)-(bx,by), its UNIT outward normal, and
// its fixity. Free edges belong in the list too -- they contribute to the
// nearest-edge tie tolerance exactly as they always have, and only the apply
// step skips them.
struct BcEdge {
    double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
    double nx = 0.0, ny = 0.0;
    EdgeFixity fixity = EdgeFixity::Free;
};

// Fix the DOFs of one mesh node for a given edge fixity. (nx, ny) is the edge
// outward normal; "normally fixed" constrains the dominant normal component.
inline void fix_node(DofMap& dofs, int node, EdgeFixity bc, double nx, double ny) {
    using BC = EdgeFixity;
    switch (bc) {
        case BC::FullyFixed:        dofs.fix_node_component(node, 0); dofs.fix_node_component(node, 1); break;
        case BC::HorizontallyFixed: dofs.fix_node_component(node, 0); break;
        case BC::VerticallyFixed:   dofs.fix_node_component(node, 1); break;
        case BC::NormallyFixed:     if (std::fabs(nx) >= std::fabs(ny)) dofs.fix_node_component(node, 0);
                                    else                                dofs.fix_node_component(node, 1); break;
        case BC::Free: default: break;
    }
}

// Apply per-edge BCs: each domain-boundary node takes the UNION of the fixities of every edge it
// lies on (within tolerance). A node at the junction of two edges -- a polygon corner, or a layer
// interface meeting the side boundary -- belongs to BOTH edges; picking only the "nearest" one is
// an arbitrary tie-break that can let a Free internal edge override a fixed side edge, leaving
// e.g. the side node at a layer boundary unconstrained (a real support hole: the K0 lateral
// traction there had nothing to push against). Union semantics match PLAXIS corner practice:
// constraints accumulate, Free contributes nothing.
// `free_base_ux` (compliant-base Dynamic phase only): a node on the model's BOTTOM extreme plane
// must keep its horizontal DOF free -- the base moves there (total-motion formulation, the Lysmer
// dashpot + upward-wave traction take over the restraint role). Only u_x is released; u_y stays
// fixed (v1 horizontal-SH scope, disclosed in the phase message). Static phases never pass it, so
// the static chain (and the rigid-base dynamic default) is bit-for-bit unchanged.
// `released` (optional, one entry per mesh node): a node marked here takes NO edge fixity at all.
// It exists for one measured reason. An interface drawn ALONG a fixed boundary splits the mesh into
// two node sets at the SAME coordinates -- one carrying the soil, one carrying nothing -- and this
// function matches by coordinate, so it fixed both and welded the joint shut. The support belongs
// to the empty side, which plays the rigid outside world; the soil side has to be free to slide,
// and only the caller knows which is which. Nodes are released wholesale rather than per component
// because the seam's two sides are geometrically indistinguishable here: the caller has already
// decided that this node is not the one the boundary is talking about.
inline void apply_boundary_conditions(const std::vector<BcEdge>& edges, const katai::mesh::Mesh& mesh,
                                      DofMap& dofs, bool free_base_ux = false,
                                      const std::vector<char>* released = nullptr) {
    double ymin = 1e300, yspan_max = -1e300;
    if (free_base_ux) {
        for (int node : mesh.boundary_nodes) {
            ymin = std::fmin(ymin, mesh.y[node]);
            yspan_max = std::fmax(yspan_max, mesh.y[node]);
        }
    }
    const double ytol = free_base_ux ? 1e-6 * std::fmax(1.0, yspan_max - ymin) : 0.0;
    for (int node : mesh.boundary_nodes) {
        if (released && node < static_cast<int>(released->size()) && (*released)[node]) continue;
        const double xn = mesh.x[node], yn = mesh.y[node];
        struct EdgeHit { EdgeFixity bc; double nx, ny; double d; };
        std::vector<EdgeHit> hits;
        hits.reserve(edges.size());
        double best_d = 1e300;
        for (const BcEdge& e : edges) {
            const double ex = e.bx - e.ax, ey = e.by - e.ay, l2 = ex * ex + ey * ey;
            if (l2 < 1e-18) continue;
            double t = ((xn - e.ax) * ex + (yn - e.ay) * ey) / l2; t = std::clamp(t, 0.0, 1.0);
            const double qx = e.ax + t * ex, qy = e.ay + t * ey;
            const double d = std::hypot(xn - qx, yn - qy);
            best_d = std::fmin(best_d, d);
            hits.push_back({e.fixity, e.nx, e.ny, d});
        }
        // "On the edge" tolerance: tie window around the nearest edge (node-edge distances are
        // either ~round-off for edges through the node or >= a fraction of the element size).
        const double tol = best_d + 1e-9 + 1e-7 * std::fmax(std::fabs(xn), std::fabs(yn));
        const bool on_base = free_base_ux && yn <= ymin + ytol;
        for (const auto& h : hits) {
            if (h.d > tol || h.bc == EdgeFixity::Free) continue;
            EdgeFixity bc = h.bc;
            if (on_base) {   // release u_x on the moving base; keep every u_y fixity
                if (bc == EdgeFixity::FullyFixed) bc = EdgeFixity::VerticallyFixed;
                else if (bc == EdgeFixity::HorizontallyFixed) continue;
                else if (bc == EdgeFixity::NormallyFixed &&
                         std::fabs(h.nx) >= std::fabs(h.ny)) continue;
            }
            fix_node(dofs, node, bc, h.nx, h.ny);
        }
    }
}

} // namespace katai::core
