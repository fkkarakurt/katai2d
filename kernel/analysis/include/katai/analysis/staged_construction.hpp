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
// Math + sources: docs/references/initial-stress-k0.md §5 (staged construction).

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

// THE SIDE OF AN ELEMENT AT A LINE'S END IS THE SECTOR IT LIES IN. Round the node where a split
// line ends, the elements were classified by the side of the line's infinite extension -- and
// where another wall leaves that node (an L corner, a vertex of a lining), the extension cuts
// across the ground the other wall bounds: an element beyond it took the copy of the other side,
// and the mesh opened along an edge no wall was drawn on (measured on an L corner of a wall and a
// slab: one edge opened under the slab). The rays of the walls that leave the node cut its
// neighbourhood into sectors, and an element takes the side of the sector it lies in.
// `own` is this line's ray from the node, `others` the other walls' rays, `ccw_positive` whether
// the sector counter-clockwise next to `own` is the positive (unchanged) side. Returns +1 positive,
// -1 negative, 0 for a sector the line does not bound (the caller keeps its own rule there).
inline int sector_side(double px, double py, double cx, double cy, double own,
                       const std::vector<double>& others, bool ccw_positive) {
    constexpr double kTwoPi = 6.28318530717958647692;
    const auto rel = [&](double a) {
        double d = std::fmod(a - own, kTwoPi);
        if (d < 0.0) d += kTwoPi;
        return d;
    };
    double lo = kTwoPi, hi = 0.0;   // the nearest other ray counter-clockwise, and clockwise
    for (const double o : others) {
        const double r = rel(o);
        if (r < 1e-9 || r > kTwoPi - 1e-9) continue;   // a ray along this line is this line
        lo = std::min(lo, r);
        hi = std::max(hi, r);
    }
    if (lo >= kTwoPi) return 0;
    const double t = rel(std::atan2(cy - py, cx - px));
    if (t < lo) return ccw_positive ? 1 : -1;
    if (t > hi) return ccw_positive ? -1 : 1;
    return 0;
}

inline std::vector<SeamPair> split_mesh_at_wall(mesh::Mesh& mesh, double x_wall,
                                                double y_toe, double y_top,
                                                double tol = 1e-6,
                                                const std::vector<double>& rays_top = {}) {
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
        double xc = 0.0, yc = 0.0;
        bool at_top = false;
        for (int k = 0; k < mesh.nodes_per_element; ++k) {
            const int nd = mesh.node_of(e, k);
            xc += mesh.x[nd]; yc += mesh.y[nd];
            at_top = at_top || (dup[nd] >= 0 && std::fabs(mesh.y[nd] - y_top) < tol);
        }
        xc /= mesh.nodes_per_element; yc /= mesh.nodes_per_element;
        // At the top, where other walls leave it: the sector (the wall's ray points down; the
        // right side is counter-clockwise of it).
        const int sec = at_top && !rays_top.empty()
                            ? sector_side(x_wall, y_top, xc, yc, -1.57079632679489661923, rays_top, true)
                            : 0;
        if (sec > 0 || (sec == 0 && xc >= x_wall)) continue;  // right side → unchanged
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
// BEFORE the DofMap is built. (Interface: a structural/in-soil slip surface.)
struct SegSeam { int orig; int dup; double x, y, s; };
// `rays_a` / `rays_b`: the directions [rad] of the other walls that leave (ax, ay) / (bx, by);
// the elements round a node split there take their side by sector (sector_side above).
inline std::vector<SegSeam> split_mesh_at_segment(mesh::Mesh& mesh, double ax, double ay,
                                                  double bx, double by, double s_lo, double s_hi,
                                                  double tol = 1e-6,
                                                  const std::vector<double>& rays_a = {},
                                                  const std::vector<double>& rays_b = {}) {
    const double L = std::hypot(bx - ax, by - ay);
    if (L < 1e-12) return {};
    const double tx = (bx - ax) / L, ty = (by - ay) / L;  // unit tangent
    const double nx = ty, ny = -tx;                        // unit normal (right-handed)
    std::vector<int> dup(mesh.node_count, -1);
    std::vector<SegSeam> seam;
    std::vector<int> on;
    std::vector<double> on_s;
    for (int i = 0; i < mesh.node_count; ++i) {
        const double rx = mesh.x[i] - ax, ry = mesh.y[i] - ay;
        const double s = rx * tx + ry * ty;       // along-segment parameter
        const double d = rx * nx + ry * ny;       // signed perpendicular distance
        if (std::fabs(d) < tol && s > s_lo - tol && s < s_hi + tol) { on.push_back(i); on_s.push_back(s); }
    }
    // A POINT ANOTHER LINE HAS ALREADY SPLIT, WHERE THIS ONE ENDS, IS NOT SPLIT AGAIN. Where two
    // split lines end at one point -- the crown of a lining drawn as chords -- the first line left
    // two nodes there, one for each side, and splitting both again made four: the second line's
    // seam held the point twice, its wall had a position that belonged to nothing, and the ring
    // was cut at the crown. Measured on a deep tunnel lined with an interface: the lining carried
    // 23 kN/m of the 2690 the closed form gives, and its two halves moved apart by 0.55 m. The
    // two nodes already there ARE the pair: the one whose elements lie on this line's positive
    // side keeps that side, the other the negative one, and no element is rewired. A point where
    // two lines CROSS is split as before -- both of its nodes again, four quadrants for four sides.
    std::vector<int> new_nodes;
    std::vector<char> paired(on.size(), 0);
    // Which side a node's elements lie on: by sector where other walls leave the node (an
    // element on the far side of another wall says nothing of this line's sides), else by the
    // signed distance from the line.
    const auto side_of = [&](int node) {
        const bool at_a = std::hypot(mesh.x[node] - ax, mesh.y[node] - ay) < tol && !rays_a.empty();
        const bool at_b = std::hypot(mesh.x[node] - bx, mesh.y[node] - by) < tol && !rays_b.empty();
        double sum = 0.0, by_sector = 0.0;
        for (int e = 0; e < mesh.element_count; ++e) {
            bool has = false;
            double cx = 0.0, cy = 0.0;
            for (int k = 0; k < mesh.nodes_per_element; ++k) {
                const int nd = mesh.node_of(e, k);
                has = has || nd == node;
                cx += mesh.x[nd]; cy += mesh.y[nd];
            }
            if (!has) continue;
            cx /= mesh.nodes_per_element; cy /= mesh.nodes_per_element;
            sum += (cx - ax) * nx + (cy - ay) * ny;
            if (at_a) by_sector += sector_side(ax, ay, cx, cy, std::atan2(ty, tx), rays_a, false);
            else if (at_b) by_sector += sector_side(bx, by, cx, cy, std::atan2(-ty, -tx), rays_b, true);
        }
        return by_sector != 0.0 ? by_sector : sum;
    };
    for (size_t a = 0; a < on.size(); ++a) {
        if (paired[a]) continue;
        size_t twin = on.size();
        int count = 0;
        for (size_t b = 0; b < on.size(); ++b)
            if (b != a && std::fabs(mesh.x[on[b]] - mesh.x[on[a]]) < tol &&
                std::fabs(mesh.y[on[b]] - mesh.y[on[a]]) < tol) { twin = b; ++count; }
        const int i = on[a];
        const bool at_end = on_s[a] < tol || on_s[a] > L - tol;
        if (count == 1 && at_end) {
            const int j = on[twin];
            const bool i_positive = side_of(i) >= side_of(j);
            seam.push_back({i_positive ? i : j, i_positive ? j : i, mesh.x[i], mesh.y[i], on_s[a]});
            paired[a] = paired[twin] = 1;
            continue;
        }
        dup[i] = mesh.node_count + static_cast<int>(new_nodes.size());
        new_nodes.push_back(i);
        seam.push_back({i, dup[i], mesh.x[i], mesh.y[i], on_s[a]});
    }
    for (const int i : new_nodes) { mesh.x.push_back(mesh.x[i]); mesh.y.push_back(mesh.y[i]); }
    mesh.node_count += static_cast<int>(new_nodes.size());
    const double own_a = std::atan2(ty, tx), own_b = std::atan2(-ty, -tx);
    for (int e = 0; e < mesh.element_count; ++e) {  // re-wire the negative-side elements to the twins
        double cx = 0.0, cy = 0.0;
        bool at_a = false, at_b = false;
        for (int k = 0; k < mesh.nodes_per_element; ++k) {
            const int nd = mesh.node_of(e, k);
            cx += mesh.x[nd]; cy += mesh.y[nd];
            if (dup[nd] >= 0) {
                at_a = at_a || std::hypot(mesh.x[nd] - ax, mesh.y[nd] - ay) < tol;
                at_b = at_b || std::hypot(mesh.x[nd] - bx, mesh.y[nd] - by) < tol;
            }
        }
        cx /= mesh.nodes_per_element; cy /= mesh.nodes_per_element;
        // At an end other walls leave: the sector. The positive side is clockwise of the line's
        // ray at its start and counter-clockwise of it at its end (the normal is the tangent
        // turned clockwise).
        int sec = 0;
        if (at_a && !rays_a.empty()) sec = sector_side(ax, ay, cx, cy, own_a, rays_a, false);
        if (sec == 0 && at_b && !rays_b.empty()) sec = sector_side(bx, by, cx, cy, own_b, rays_b, true);
        if (sec > 0) continue;
        if (sec == 0 && (cx - ax) * nx + (cy - ay) * ny >= 0.0) continue;   // positive side → unchanged
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
// The case this exists for is a structural benchmark built inside a soil model: a beam is
// modelled by deactivating the soil cluster so that only the beams remain.
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
