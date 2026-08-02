#pragma once
// Hydraulic boundary description (Stage B4: extracted from the application
// driver). Three engine services over neutral inputs:
//
//   - the phreatic-surface elevation along a water polyline;
//   - the drainage mask of a time-dependent flow phase (consolidation /
//     fully-coupled);
//   - the prescribed-head node set of a transient flow phase.
//
// The flow contract mirrors Stage B3's displacement BCs: the engine takes a
// neutral list of flow edges -- segment, kind, head -- and owns the node
// matching with its tie tolerance and the corner rules. Closed edges belong in
// the list exactly like Free edges did in B3: every edge always contributed to
// the nearest-edge tie window, and only the apply step filters by kind.
//
// One asymmetry is deliberate and preserved from the shipped behaviour:
// flow_drained_nodes takes `have_flow_bcs` as a separate statement from the
// caller rather than deriving it from the edge list. The original scan counted
// flow-BC declarations even on polygons too degenerate to produce a matchable
// edge; on such a model the top-drain fallback must NOT silently activate, and
// an edge-list-derived flag would activate it.

#include <algorithm>
#include <cmath>
#include <vector>

#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Kind of one flow boundary edge. The engine's own vocabulary; the schema enum
// is mapped once, at the caller's seam.
enum class FlowEdgeKind { Closed, Head, Seepage };

// One flow boundary edge: the segment (ax,ay)-(bx,by), its kind, and -- for a
// Head edge -- the prescribed hydraulic head [m elevation].
struct FlowEdge {
    double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
    FlowEdgeKind kind = FlowEdgeKind::Closed;
    double head = 0.0;
};

// Phreatic-surface elevation at horizontal position x, interpolated along the
// (possibly sloped) water polyline. Returns a very low value when the polyline
// cannot define a surface (fewer than two points -> everywhere "above table").
inline double phreatic_surface_at(const std::vector<double>& wx,
                                  const std::vector<double>& wy, double x) {
    if (wx.size() < 2) return -1e30;
    const int n = (int)wx.size();
    // Polyline may run either direction; scan segments and interpolate by x.
    double y_at = wy[0]; double best_d = 1e300;
    for (int i = 0; i + 1 < n; ++i) {
        const double xa = wx[i], ya = wy[i], xb = wx[i + 1], yb = wy[i + 1];
        const double lo = std::min(xa, xb), hi = std::max(xa, xb);
        if (x >= lo - 1e-9 && x <= hi + 1e-9) {
            const double dx = xb - xa;
            const double t = std::fabs(dx) < 1e-12 ? 0.0 : (x - xa) / dx;
            return ya + t * (yb - ya);
        }
        // Outside the polyline span: remember the nearest endpoint to clamp to.
        const double da = std::fabs(x - xa); if (da < best_d) { best_d = da; y_at = ya; }
        const double db = std::fabs(x - xb); if (db < best_d) { best_d = db; y_at = yb; }
    }
    return y_at;   // clamp to nearest end elevation beyond the polyline span
}

// Drainage-boundary mask for a time-dependent flow phase (consolidation / fully-coupled): a boundary
// node on a prescribed-head or seepage flow edge drains (excess pore = 0). With NO flow BCs declared
// (`have_flow_bcs` false) the model top is taken as draining (sensible foundation/embankment
// default); everything else is impermeable (natural). Nodes touched only by inactive (excavated)
// elements carry no pore DOF -> drained.
inline std::vector<char> flow_drained_nodes(const std::vector<FlowEdge>& edges, bool have_flow_bcs,
                                            const katai::mesh::Mesh& mesh,
                                            const std::vector<char>& act, double yscale) {
    std::vector<char> drained(mesh.node_count, 0);
    if (have_flow_bcs) {
        for (int node : mesh.boundary_nodes) {
            const double xn = mesh.x[node], yn = mesh.y[node];
            double best_d = 1e300;
            struct Hit { FlowEdgeKind kind; double d; };
            std::vector<Hit> hits;
            hits.reserve(edges.size());
            for (const FlowEdge& e : edges) {
                const double ex = e.bx - e.ax, ey = e.by - e.ay, l2 = ex * ex + ey * ey;
                if (l2 < 1e-18) continue;
                double t = ((xn - e.ax) * ex + (yn - e.ay) * ey) / l2; t = std::clamp(t, 0.0, 1.0);
                const double d = std::hypot(xn - (e.ax + t * ex), yn - (e.ay + t * ey));
                best_d = std::fmin(best_d, d); hits.push_back({e.kind, d});
            }
            const double tol = best_d + 1e-9 + 1e-7 * std::fmax(std::fabs(xn), std::fabs(yn));
            for (const auto& h : hits)
                if (h.d <= tol && (h.kind == FlowEdgeKind::Head ||
                                   h.kind == FlowEdgeKind::Seepage)) drained[node] = 1;
        }
    } else {
        double ytop = -1e30;
        for (int node : mesh.boundary_nodes) ytop = std::fmax(ytop, mesh.y[node]);
        for (int node : mesh.boundary_nodes)
            if (mesh.y[node] > ytop - 1e-6 * std::fmax(1.0, yscale)) drained[node] = 1;
    }
    if (!act.empty()) {
        std::vector<char> touched(mesh.node_count, 0);
        const int npe = mesh.nodes_per_element;
        for (int e = 0; e < mesh.element_count; ++e)
            if (act[e]) for (int k = 0; k < npe; ++k) touched[mesh.node_of(e, k)] = 1;
        for (int n = 0; n < mesh.node_count; ++n) if (!touched[n]) drained[n] = 1;
    }
    return drained;
}

// Prescribed-head boundary for a transient flow phase (PLAXIS GroundwaterFlow head BC): a boundary
// node on a Head flow edge is prescribed to that edge's head [m elevation]; the nearest Head edge
// wins at corners. Fills `is_presc` (node) + `head_val` (node).
inline void flow_head_nodes(const std::vector<FlowEdge>& edges, const katai::mesh::Mesh& mesh,
                            std::vector<char>& is_presc, std::vector<double>& head_val) {
    is_presc.assign(mesh.node_count, 0);
    head_val.assign(mesh.node_count, 0.0);
    for (int node : mesh.boundary_nodes) {
        const double xn = mesh.x[node], yn = mesh.y[node];
        double best_d = 1e300;
        struct Hit { FlowEdgeKind kind; double head, d; };
        std::vector<Hit> hits;
        hits.reserve(edges.size());
        for (const FlowEdge& e : edges) {
            const double ex = e.bx - e.ax, ey = e.by - e.ay, l2 = ex * ex + ey * ey;
            if (l2 < 1e-18) continue;
            double t = ((xn - e.ax) * ex + (yn - e.ay) * ey) / l2; t = std::clamp(t, 0.0, 1.0);
            const double d = std::hypot(xn - (e.ax + t * ex), yn - (e.ay + t * ey));
            best_d = std::fmin(best_d, d);
            hits.push_back({e.kind, e.head, d});
        }
        const double tol = best_d + 1e-9 + 1e-7 * std::fmax(std::fabs(xn), std::fabs(yn));
        double hval = 0.0, hd = 1e300; bool has_head = false;
        for (const auto& h : hits)
            if (h.d <= tol && h.kind == FlowEdgeKind::Head && h.d < hd) {
                has_head = true; hval = h.head; hd = h.d;
            }
        if (has_head) { is_presc[node] = 1; head_val[node] = hval; }
    }
}

} // namespace katai::core
