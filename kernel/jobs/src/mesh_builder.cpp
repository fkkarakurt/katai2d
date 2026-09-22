// The model-to-mesh builder body, compiled ONCE (section 5.2).
#include <katai/jobs/mesh_builder.hpp>

#include <algorithm>
#include <cstdint>
#include <unordered_set>

#include <katai/geometry/planar_graph.hpp>
#include <katai/mesh/delaunay.hpp>

namespace katai::app {

MeshResult mesh_from_project(const model::Project& pr, double max_area,
                             int order, const MeshOptions& opt) {
    const double min_angle_deg = opt.min_angle_deg;
    MeshResult R;
    if (pr.polygons.empty()) { R.message = "No soil polygons to mesh."; return R; }
    if (max_area <= 0.0) { R.message = "Element size must be positive."; return R; }

    // 1) Build a CLEAN planar straight-line graph (PSLG). Raw input segments are the soil polygon
    // edges plus the internal lines the mesh has to follow (plates, geogrids, interfaces, line
    // loads, prescribed displacements, wells and drains). They are noded ONCE, to one tolerance
    // scaled to the model (katai/geometry/planar_graph.hpp): near-coincident points merge, a
    // point within tolerance of a segment splits it, crossings become vertices and a piece that
    // two sources share becomes one segment. The result has no crossings, no vertex inside a
    // segment and no feature below the tolerance -- the constrained-Delaunay precondition, and
    // the condition for the refinement to terminate: a corner left 1e-5 off its neighbour is a
    // sliver 1e-5 thick that Ruppert refinement would try to resolve and never finish.
    //
    // Which pieces are kept is decided by who OWNS each side of it, with the rule the rest of
    // the tree already applies to a point (material, phase activity): the last polygon that
    // contains it. A polygon piece separating two different owners is a constraint, and when one
    // side has no owner it is also the domain OUTLINE (the inside/outside test uses only those,
    // so they form the boundary of the union of the polygons, each piece once). A piece with the
    // same owner on both sides -- a boundary drawn over by a later polygon, a shared edge traced
    // twice -- bounds nothing and is dropped. Overlapping polygons therefore mesh as the later
    // one over the earlier, a lens drawn inside a layer meshes as a lens, and neither becomes a
    // hole; the old even-odd count over every polygon edge made both of them holes.
    //   - structural lines are kept where they run through soil (either side owned), so a plate
    //     poking out of the soil contributes only its in-soil part;
    //   - load/displacement/hydraulic lines are kept whole (a surface load lies on the boundary).
    enum RawKind { kPolyEdge, kStruct, kLine };
    struct Raw { double ax, ay, bx, by; RawKind kind; };
    std::vector<Raw> raw;
    std::vector<int> poly_first(pr.polygons.size(), -1), poly_count(pr.polygons.size(), 0);
    for (size_t p = 0; p < pr.polygons.size(); ++p) {
        const auto& P = pr.polygons[p];
        const int n = (int)P.x.size(); if (n < 3) continue;
        poly_first[p] = (int)raw.size(); poly_count[p] = n;
        for (int k = 0; k < n; ++k)
            raw.push_back({P.x[k], P.y[k], P.x[(k + 1) % n], P.y[(k + 1) % n], kPolyEdge});
    }
    for (const auto& s : pr.structs) {
        using SK = model::StructKind;
        // Plates/geogrids AND standalone interfaces become conforming constraints: the interface line
        // must lie on mesh edges so split_mesh_at_segment can duplicate its nodes (a slip surface).
        if (s.kind == SK::Plate || s.kind == SK::Geogrid || s.kind == SK::Interface)
            raw.push_back({s.x1, s.y1, s.x2, s.y2, kStruct});
    }
    // Distributed (line) loads become conforming constraints too, so the solver can assemble the
    // surcharge as consistent nodal forces along the resulting edge chain (build_problem).
    for (const auto& L : pr.loads)
        if (L.kind == model::LoadKind::Distributed)
            raw.push_back({L.x1, L.y1, L.x2, L.y2, kLine});
    // Prescribed-displacement lines likewise: the imposed components apply to the mesh
    // nodes ON the line, so the line must lie on mesh edges.
    for (const auto& D : pr.disps)
        raw.push_back({D.x1, D.y1, D.x2, D.y2, kLine});
    // Wells and drains for the same reason: a well prescribes a discharge along its line and a
    // drain a head at its nodes, so the line has to BE a chain of mesh edges. A hydraulic
    // condition the mesh does not follow is one whose water is applied somewhere else.
    for (const auto& H : pr.hydros)
        raw.push_back({H.x1, H.y1, H.x2, H.y2, kLine});
    if (raw.size() < 3) { R.message = "Degenerate geometry."; return R; }

    double minx = raw[0].ax, maxx = minx, miny = raw[0].ay, maxy = miny;
    for (const Raw& g : raw) {
        for (double v : {g.ax, g.ay, g.bx, g.by})
            if (!std::isfinite(v)) { R.message = "The geometry has a coordinate that is not a finite number."; return R; }
        minx = std::min({minx, g.ax, g.bx}); maxx = std::max({maxx, g.ax, g.bx});
        miny = std::min({miny, g.ay, g.by}); maxy = std::max({maxy, g.ay, g.by});
    }
    std::vector<std::array<geometry::V2, 2>> in_segs;
    in_segs.reserve(raw.size());
    for (const Raw& g : raw) in_segs.push_back({geometry::V2{g.ax, g.ay}, geometry::V2{g.bx, g.by}});
    const geometry::PlanarGraph G =
        geometry::planar_graph(in_segs, geometry::geometry_tolerance(minx, miny, maxx, maxy));

    // Each polygon re-read on the graph and split into simple loops (a pinched layer into its
    // lobes); a loop knows its orientation, so the side of a piece it owns is exact -- no probe
    // point is offset from the piece, which a narrow wedge would defeat.
    struct Loop {
        std::vector<int> ring; bool ccw;
        double x0, y0, x1, y1;
        std::unordered_set<std::uint64_t> walked;   // directed (u -> v) pieces
    };
    const auto dkey = [](int u, int v) { return ((std::uint64_t)(std::uint32_t)u << 32) | (std::uint32_t)v; };
    std::vector<std::vector<Loop>> loops(pr.polygons.size());
    for (size_t p = 0; p < pr.polygons.size(); ++p) {
        if (poly_first[p] < 0) continue;
        std::vector<int> walk;
        for (int k = 0; k < poly_count[p]; ++k)
            for (int v : G.chain_nodes(poly_first[p] + k))
                if (walk.empty() || walk.back() != v) walk.push_back(v);
        for (auto& ring : geometry::simple_loops(walk)) {
            Loop L; L.ring = std::move(ring);
            L.ccw = geometry::ring_area(G.nodes, L.ring) > 0.0;
            L.x0 = L.x1 = G.nodes[L.ring[0]].x; L.y0 = L.y1 = G.nodes[L.ring[0]].y;
            for (size_t k = 0; k < L.ring.size(); ++k) {
                const auto& q = G.nodes[L.ring[k]];
                L.x0 = std::min(L.x0, q.x); L.x1 = std::max(L.x1, q.x);
                L.y0 = std::min(L.y0, q.y); L.y1 = std::max(L.y1, q.y);
                L.walked.insert(dkey(L.ring[k], L.ring[(k + 1) % L.ring.size()]));
            }
            loops[p].push_back(std::move(L));
        }
    }
    // Owner (last containing polygon, -1 = none) on the left and right of piece a -> b.
    const auto owners = [&](int a, int b, int& left, int& right) {
        left = right = -1;
        const double mx = 0.5 * (G.nodes[a].x + G.nodes[b].x), my = 0.5 * (G.nodes[a].y + G.nodes[b].y);
        for (size_t p = 0; p < loops.size(); ++p) {
            bool cl = false, cr = false;   // even-odd over the polygon's loops
            for (const Loop& L : loops[p]) {
                if (L.walked.count(dkey(a, b)))      { cl ^= L.ccw;  cr ^= !L.ccw; }
                else if (L.walked.count(dkey(b, a))) { cl ^= !L.ccw; cr ^= L.ccw; }
                else if (mx >= L.x0 && mx <= L.x1 && my >= L.y0 && my <= L.y1 &&
                         geometry::ring_contains(G.nodes, L.ring, mx, my)) { cl = !cl; cr = !cr; }
            }
            if (cl) left = (int)p;
            if (cr) right = (int)p;
        }
    };

    // Emit the kept pieces in input order (segment by segment, along each segment), so a model
    // with nothing to node hands the mesher exactly the points and segments it always did.
    std::vector<double> px, py;
    std::vector<std::array<int, 2>> segs, outline;
    std::vector<int> pt_of(G.nodes.size(), -1);
    std::vector<char> decided(G.edges.size(), 0), keep(G.edges.size(), 0), bound(G.edges.size(), 0);
    const auto add_pt = [&](int v) {
        if (pt_of[v] < 0) { pt_of[v] = (int)px.size(); px.push_back(G.nodes[v].x); py.push_back(G.nodes[v].y); }
        return pt_of[v];
    };
    for (int i = 0; i < (int)raw.size(); ++i)
        for (const auto& st : G.chain[i]) {
            const int e = st.edge;
            if (!decided[e]) {
                decided[e] = 1;
                bool has_poly = false, has_struct = false, has_line = false;
                for (int s : G.edges[e].src) {
                    has_poly |= raw[s].kind == kPolyEdge;
                    has_struct |= raw[s].kind == kStruct;
                    has_line |= raw[s].kind == kLine;
                }
                int ol, orr;
                owners(G.edges[e].a, G.edges[e].b, ol, orr);
                keep[e] = (has_poly && ol != orr) || (has_struct && (ol >= 0 || orr >= 0)) || has_line;
                bound[e] = has_poly && ((ol < 0) != (orr < 0));
                if (keep[e]) {
                    const int u = st.forward ? G.edges[e].a : G.edges[e].b;
                    const int w = st.forward ? G.edges[e].b : G.edges[e].a;
                    const int pa = add_pt(u), pb = add_pt(w);
                    segs.push_back({pa, pb});
                    if (bound[e]) outline.push_back({pa, pb});
                }
            }
        }
    if (px.size() < 3 || outline.size() < 3) { R.message = "Degenerate geometry."; return R; }

    // 2) Local mesh density: build the sizing field from the per-object coarseness factors
    // (a factor scales the target element size). When everything is at the default (no factors,
    // no auto-refine), the field is skipped entirely and the constant-max_area mesher runs --
    // bit-identical meshes.
    const auto clampf = [](double f) { return std::clamp(f, 1.0 / 16.0, 4.0); };
    struct SizeSrc { double ax, ay, bx, by, h; };   // segment source (point: a == b)
    std::vector<SizeSrc> srcs;
    const double h0 = std::sqrt(2.0 * max_area);
    for (const auto& s : pr.structs) {
        using SK = model::StructKind;
        if (s.kind != SK::Plate && s.kind != SK::Geogrid && s.kind != SK::EmbeddedBeam) continue;
        const double f = clampf(s.coarseness * (opt.auto_refine ? opt.auto_factor : 1.0));
        if (f < 1.0 - 1e-12)
            srcs.push_back({s.x1, s.y1, s.x2, s.y2, h0 * f});
    }
    for (const auto& L : pr.loads) {
        const double f = clampf(L.coarseness * (opt.auto_refine ? opt.auto_factor : 1.0));
        if (f < 1.0 - 1e-12) {
            if (L.kind == model::LoadKind::Distributed) srcs.push_back({L.x1, L.y1, L.x2, L.y2, h0 * f});
            else srcs.push_back({L.x1, L.y1, L.x1, L.y1, h0 * f});
        }
    }
    for (const auto& D : pr.disps) {
        const double f = clampf(D.coarseness * (opt.auto_refine ? opt.auto_factor : 1.0));
        if (f < 1.0 - 1e-12) srcs.push_back({D.x1, D.y1, D.x2, D.y2, h0 * f});
    }
    for (const auto& H : pr.hydros) {
        const double f = clampf(H.coarseness * (opt.auto_refine ? opt.auto_factor : 1.0));
        if (f < 1.0 - 1e-12) srcs.push_back({H.x1, H.y1, H.x2, H.y2, h0 * f});
    }
    bool region_factors = false;
    for (const auto& P : pr.polygons)
        if (std::fabs(P.coarseness - 1.0) > 1e-12) region_factors = true;

    katai::mesh::SizeField field;   // empty => constant max_area (legacy path, bit-identical)
    if (!srcs.empty() || region_factors) {
        field = [&pr, srcs, h0, clampf, grading = opt.grading](double qx, double qy) {
            double f_region = 1.0;
            for (const auto& P : pr.polygons)
                if (point_in_polygon(qx, qy, P)) f_region = clampf(P.coarseness);   // last wins
            double h = h0 * f_region;
            for (const auto& s : srcs) {
                const double dx = s.bx - s.ax, dy = s.by - s.ay, l2 = dx * dx + dy * dy;
                double t = l2 > 1e-18 ? ((qx - s.ax) * dx + (qy - s.ay) * dy) / l2 : 0.0;
                t = std::clamp(t, 0.0, 1.0);
                const double d = std::hypot(qx - (s.ax + t * dx), qy - (s.ay + t * dy));
                h = std::fmin(h, s.h + grading * d);
            }
            return 0.5 * h * h;
        };
    }

    // 3) Ruppert quality mesh of the clean PSLG (outline = domain boundary for the in/out test).
    const katai::mesh::Triangulation T =
        field ? katai::mesh::quality_mesh(px, py, segs, min_angle_deg, field, outline)
              : katai::mesh::quality_mesh(px, py, segs, min_angle_deg, max_area, outline);

    // 4) Keep triangles whose centroid lies inside a soil polygon; tag its material.
    std::vector<std::array<int, 3>> kept; std::vector<int> kept_mat;
    for (const auto& t : T.triangles) {
        const double cx = (T.x[t[0]] + T.x[t[1]] + T.x[t[2]]) / 3.0;
        const double cy = (T.y[t[0]] + T.y[t[1]] + T.y[t[2]]) / 3.0;
        int mat = -2;
        for (const auto& P : pr.polygons) if (point_in_polygon(cx, cy, P)) mat = P.material;   // last wins
        if (mat == -2) continue;                                  // outside the soil
        kept.push_back(t); kept_mat.push_back(mat < 0 ? 0 : mat);
    }
    if (kept.empty()) { R.message = "Meshed region is empty."; return R; }

    // 5) Compact the vertices actually used (avoid orphan -> singular nodes).
    std::vector<int> remap(T.x.size(), -1);
    katai::mesh::Triangulation C;
    for (auto& t : kept) {
        std::array<int, 3> nt;
        for (int c = 0; c < 3; ++c) {
            const int v = t[c];
            if (remap[v] < 0) { remap[v] = (int)C.x.size(); C.x.push_back(T.x[v]); C.y.push_back(T.y[v]); }
            nt[c] = remap[v];
        }
        C.triangles.push_back(nt);
    }
    C.point_count = (int)C.x.size();

    // 6) Promote to a quadratic / quartic FE mesh.
    katai::mesh::Mesh m = (order == 15) ? katai::mesh::tri15_from_triangulation(C, 0)
                                        : katai::mesh::tri6_from_triangulation(C, 0);

    // 7) Per-element material from the corner-triangle centroid.
    for (int e = 0; e < m.element_count; ++e) {
        const int n0 = m.node_of(e, 0), n1 = m.node_of(e, 1), n2 = m.node_of(e, 2);
        const double cx = (m.x[n0] + m.x[n1] + m.x[n2]) / 3.0, cy = (m.y[n0] + m.y[n1] + m.y[n2]) / 3.0;
        int mat = 0;
        for (const auto& P : pr.polygons) if (point_in_polygon(cx, cy, P)) mat = (P.material < 0 ? 0 : P.material);
        m.element_material[e] = mat;
    }

    R.mesh = std::move(m); R.ok = true;
    R.quality_met = T.quality_met;
    R.min_angle_asked = min_angle_deg;
    R.corner_elements = T.corner_elements;
    R.message = "Mesh: " + std::to_string(R.mesh.node_count) + " nodes, " +
                std::to_string(R.mesh.element_count) + " elements.";
    // A mesh that did not reach its own quality bound says so IN THE MESSAGE, because the message
    // is what every front end already shows. Elements below the bound are not wrong, they are
    // worse conditioned than the run was told to accept, and which of those it is has to be the
    // reader's to make.
    if (!R.quality_met)
        R.message += " The refinement stopped at its step cap before every element met the " +
                     std::to_string((int)min_angle_deg) +
                     "-degree minimum angle, so some elements are worse conditioned than asked " +
                     "for. Coarsen the element size, or simplify the geometry where it has very " +
                     "sharp corners, and mesh again before reading the result as converged.";
    else if (R.corner_elements > 0)
        R.message += " " + std::to_string(R.corner_elements) +
                     (R.corner_elements == 1 ? " element sits" : " elements sit") +
                     " in a corner of the geometry narrower than the " +
                     std::to_string((int)min_angle_deg) + "-degree minimum angle, and keeps" +
                     " that corner's own angle: no mesh can do better there. If no corner is" +
                     " meant to be that sharp, look for a vertex slightly off a straight line.";
    return R;
}

}  // namespace katai::app
