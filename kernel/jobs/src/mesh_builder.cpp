// The model-to-mesh builder body, compiled ONCE (section 5.2).
#include <katai/jobs/mesh_builder.hpp>

#include <katai/mesh/delaunay.hpp>

namespace katai::app {

MeshResult mesh_from_project(const model::Project& pr, double max_area,
                             int order, const MeshOptions& opt) {
    const double min_angle_deg = opt.min_angle_deg;
    MeshResult R;
    if (pr.polygons.empty()) { R.message = "No soil polygons to mesh."; return R; }
    if (max_area <= 0.0) { R.message = "Element size must be positive."; return R; }

    // 1) Build a CLEAN planar straight-line graph (PSLG). Raw input segments are the soil polygon
    // edges (the domain OUTLINE) plus the internal structural lines (plates/geogrids). They are then
    // split at every pairwise crossing and at every point where one segment's endpoint touches
    // another's interior, so the result has no crossings and no vertex inside a segment -- exactly the
    // constrained-Delaunay precondition. Without this, a plate ending on/crossing the boundary (or two
    // lines crossing) violates the precondition and wrecks the mesh depending on placement.
    //   - `outline` (= split soil-polygon edges) is the domain boundary used for the inside/outside
    //     test; an internal line must never act as a boundary (or it deletes the mesh on one side).
    //   - structural sub-segments are clipped to the soil (only kept where their midpoint is inside),
    //     so a plate poking out of the soil contributes only its in-soil part.
    //   - `clip` (structural lines only) keeps just the in-soil part of a line poking out of the
    //     soil; outline edges and LOAD lines are kept whole (a surface distributed load lies on the
    //     boundary, where the soil-interior midpoint test is ambiguous).
    struct Raw { double ax, ay, bx, by; bool outline; bool clip; };
    std::vector<Raw> raw;
    for (const auto& P : pr.polygons) {
        const int n = (int)P.x.size(); if (n < 3) continue;
        for (int k = 0; k < n; ++k)
            raw.push_back({P.x[k], P.y[k], P.x[(k + 1) % n], P.y[(k + 1) % n], true, false});
    }
    for (const auto& s : pr.structs) {
        using SK = model::StructKind;
        // Plates/geogrids AND standalone interfaces become conforming constraints: the interface line
        // must lie on mesh edges so split_mesh_at_segment can duplicate its nodes (PLAXIS slip surface).
        if (s.kind == SK::Plate || s.kind == SK::Geogrid || s.kind == SK::Interface)
            raw.push_back({s.x1, s.y1, s.x2, s.y2, false, true});
    }
    // Distributed (line) loads become conforming constraints too, so the solver can assemble the
    // surcharge as consistent nodal forces along the resulting edge chain (build_problem).
    for (const auto& L : pr.loads)
        if (L.kind == model::LoadKind::Distributed)
            raw.push_back({L.x1, L.y1, L.x2, L.y2, false, false});
    // Prescribed-displacement lines likewise: the imposed components apply to the mesh
    // nodes ON the line, so the line must lie on mesh edges (same dedup rule as loads
    // when the line coincides with a boundary edge).
    for (const auto& D : pr.disps)
        raw.push_back({D.x1, D.y1, D.x2, D.y2, false, false});
    // Wells and drains for the same reason: a well prescribes a discharge along its line and a
    // drain a head at its nodes, so the line has to BE a chain of mesh edges. A hydraulic
    // condition the mesh does not follow is one whose water is applied somewhere else.
    for (const auto& H : pr.hydros)
        raw.push_back({H.x1, H.y1, H.x2, H.y2, false, false});
    if (raw.size() < 3) { R.message = "Degenerate geometry."; return R; }

    // Pairwise split parameters (segment p->p2 vs q->q2 proper/touching intersection).
    auto seg_x = [](double px0, double py0, double px1, double py1, double qx0, double qy0,
                    double qx1, double qy1, double& t, double& u) {
        const double rx = px1 - px0, ry = py1 - py0, sx = qx1 - qx0, sy = qy1 - qy0;
        const double rxs = rx * sy - ry * sx;
        if (std::fabs(rxs) < 1e-12) return false;          // parallel / collinear (overlaps ignored)
        const double qpx = qx0 - px0, qpy = qy0 - py0;
        t = (qpx * sy - qpy * sx) / rxs;
        u = (qpx * ry - qpy * rx) / rxs;
        return t >= -1e-9 && t <= 1.0 + 1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9;
    };
    const int nr = (int)raw.size();
    std::vector<std::vector<double>> param(nr);
    for (int i = 0; i < nr; ++i) { param[i].push_back(0.0); param[i].push_back(1.0); }
    // COLLINEAR OVERLAP split: seg_x only handles proper (transversal) crossings; two segments on the
    // SAME line (e.g. a surface distributed load lying on a boundary edge, or two collinear plates) are
    // parallel so seg_x skips them -> the longer segment then contains the other's endpoints in its
    // interior, which VIOLATES the constrained-Delaunay precondition (no vertex inside a segment) and
    // corrupts the mesh there (triangles spill outside). Fix: split each collinear segment at the
    // other's endpoints projected onto it, so every shared point becomes a shared vertex.
    auto split_collinear = [&](int i, int j) {
        const double rx = raw[i].bx - raw[i].ax, ry = raw[i].by - raw[i].ay, rr = rx * rx + ry * ry;
        const double sx = raw[j].bx - raw[j].ax, sy = raw[j].by - raw[j].ay, ss = sx * sx + sy * sy;
        if (rr < 1e-18 || ss < 1e-18) return;
        const double rlen = std::sqrt(rr);
        // j must lie on i's infinite line (perpendicular distance ~ 0) to be collinear.
        auto off_line = [&](double qx, double qy) {
            return std::fabs((qx - raw[i].ax) * ry - (qy - raw[i].ay) * rx) / rlen;
        };
        if (off_line(raw[j].ax, raw[j].ay) > 1e-7 || off_line(raw[j].bx, raw[j].by) > 1e-7) return;
        const double ti0 = ((raw[j].ax - raw[i].ax) * rx + (raw[j].ay - raw[i].ay) * ry) / rr;
        const double ti1 = ((raw[j].bx - raw[i].ax) * rx + (raw[j].by - raw[i].ay) * ry) / rr;
        const double uj0 = ((raw[i].ax - raw[j].ax) * sx + (raw[i].ay - raw[j].ay) * sy) / ss;
        const double uj1 = ((raw[i].bx - raw[j].ax) * sx + (raw[i].by - raw[j].ay) * sy) / ss;
        for (double t : {ti0, ti1}) if (t > 1e-9 && t < 1.0 - 1e-9) param[i].push_back(t);
        for (double u : {uj0, uj1}) if (u > 1e-9 && u < 1.0 - 1e-9) param[j].push_back(u);
    };
    for (int i = 0; i < nr; ++i)
        for (int j = i + 1; j < nr; ++j) {
            double t, u;
            if (seg_x(raw[i].ax, raw[i].ay, raw[i].bx, raw[i].by,
                      raw[j].ax, raw[j].ay, raw[j].bx, raw[j].by, t, u)) {
                param[i].push_back(std::clamp(t, 0.0, 1.0));
                param[j].push_back(std::clamp(u, 0.0, 1.0));
            } else {
                split_collinear(i, j);   // parallel: handle the collinear-overlap case
            }
        }

    // NOTE (measured 2026-08-13, kept so it is not retried): the embedded beam's CONNECTION POINT
    // is deliberately NOT injected here as an input vertex. Carrying it looked like the exact way
    // to make a hinged connection a degree-of-freedom identity, but an isolated interior point
    // cascades through Ruppert refinement: on the test_gui_solve fixture one such point cost
    // **+146 nodes, more than the +126 of a plate CONFORMED along the whole 8 m shaft**. Paying
    // more to stay non-conforming than conforming would have cost is not a trade worth making.
    // The connection is tied to the nearest existing node in the driver instead, which is the
    // convention this tree already applies to every structural point attachment (a point load,
    // an anchor end), diagnostic included. See docs/references/embedded-beam-formulation.md 7.1.

    std::vector<double> px, py;
    std::vector<std::array<int, 2>> segs, outline;
    const double tol = 1e-6;
    auto add_pt = [&](double x, double y) -> int {
        for (size_t i = 0; i < px.size(); ++i)
            if (std::fabs(px[i] - x) < tol && std::fabs(py[i] - y) < tol) return (int)i;
        px.push_back(x); py.push_back(y); return (int)px.size() - 1;
    };
    for (int i = 0; i < nr; ++i) {
        auto& pr2 = param[i];
        std::sort(pr2.begin(), pr2.end());
        pr2.erase(std::unique(pr2.begin(), pr2.end(),
                              [](double a, double b) { return std::fabs(a - b) < 1e-7; }), pr2.end());
        const Raw& g = raw[i];
        for (size_t k = 0; k + 1 < pr2.size(); ++k) {
            const double t0 = pr2[k], t1 = pr2[k + 1];
            const double x0 = g.ax + t0 * (g.bx - g.ax), y0 = g.ay + t0 * (g.by - g.ay);
            const double x1 = g.ax + t1 * (g.bx - g.ax), y1 = g.ay + t1 * (g.by - g.ay);
            if (g.clip) {   // clip internal structural sub-segments to the soil
                const double mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
                bool inside = false;
                for (const auto& P : pr.polygons) if (point_in_polygon(mx, my, P)) inside = true;
                if (!inside) continue;
            }
            const int a = add_pt(x0, y0), b = add_pt(x1, y1);
            if (a == b) continue;
            // A LOAD sub-segment that coincides with an existing constraint (e.g. a surface
            // surcharge collinear with the top outline edge) is redundant -- the outline already
            // conforms the boundary, so drop the duplicate (it would otherwise perturb the mesh).
            // Only loads are deduped: adjacent polygons legitimately add their SHARED edge twice
            // (the inside/outside ray-cast parity relies on it), and structural lines keep their
            // original handling. Internal (non-coincident) load lines are still added -> conforming.
            if (!g.outline && !g.clip) {
                bool dup = false;
                for (const auto& e : segs)
                    if ((e[0] == a && e[1] == b) || (e[0] == b && e[1] == a)) { dup = true; break; }
                if (dup) continue;
            }
            segs.push_back({a, b});
            if (g.outline) outline.push_back({a, b});
        }
    }
    if (px.size() < 3 || outline.size() < 3) { R.message = "Degenerate geometry."; return R; }

    // 2) Local mesh density: build the sizing field from the per-object coarseness factors
    // (PLAXIS semantics). When everything is at the default (no factors, no auto-refine), the
    // field is skipped entirely and the constant-max_area mesher runs -- bit-identical meshes.
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
    R.message = "Mesh: " + std::to_string(R.mesh.node_count) + " nodes, " +
                std::to_string(R.mesh.element_count) + " elements.";
    return R;
}

}  // namespace katai::app
