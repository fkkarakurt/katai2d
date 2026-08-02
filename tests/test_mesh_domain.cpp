// Mesh-domain integrity: NO element may spill outside the soil polygon. The user reports mesh
// triangles poking outside the shape when a load / structural line is added, or in post-processing.
// A correct constrained Delaunay never lets a triangle cross an outline edge, so every element's
// corners AND edge-midpoints AND centroid must lie inside (or on) some soil polygon -- regardless of
// internal constraint lines (loads / plates / geogrids) and for tri6 and tri15.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

namespace {
int g_failures = 0;

// Inside OR within tol of a polygon edge (boundary nodes are legitimately "on" the domain).
bool inside_or_on(double x, double y, const std::vector<m::SoilPolygon>& polys, double tol) {
    for (const auto& P : polys) if (katai::app::point_in_polygon(x, y, P)) return true;
    for (const auto& P : polys) {
        const int n = (int)P.x.size();
        for (int i = 0, j = n - 1; i < n; j = i++) {
            const double ax = P.x[j], ay = P.y[j], bx = P.x[i], by = P.y[i];
            const double dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy;
            double t = L2 > 1e-18 ? ((x - ax) * dx + (y - ay) * dy) / L2 : 0.0;
            t = t < 0 ? 0 : (t > 1 ? 1 : t);
            const double px = ax + t * dx, py = ay + t * dy;
            if (std::hypot(x - px, y - py) < tol) return true;
        }
    }
    return false;
}

// Mesh the project and report how many elements spill outside the domain (any corner, edge-mid or
// centroid strictly outside). Returns the count; 0 == clean.
int spill_count(const m::Project& pr, int order, const char* tag) {
    const auto M = katai::app::mesh_from_project(pr, 1.0, order);
    if (!M.ok) { std::printf("FAIL: %-40s mesh failed (%s)\n", tag, M.message.c_str()); ++g_failures; return -1; }
    const auto& mesh = M.mesh;
    const double tol = 1e-4;
    int bad = 0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        const double xa = mesh.x[a], ya = mesh.y[a], xb = mesh.x[b], yb = mesh.y[b], xc = mesh.x[c], yc = mesh.y[c];
        const double cx = (xa + xb + xc) / 3.0, cy = (ya + yb + yc) / 3.0;
        const double pts[7][2] = {{xa, ya}, {xb, yb}, {xc, yc},
            {0.5 * (xa + xb), 0.5 * (ya + yb)}, {0.5 * (xb + xc), 0.5 * (yb + yc)},
            {0.5 * (xc + xa), 0.5 * (yc + ya)}, {cx, cy}};
        for (auto& p : pts) if (!inside_or_on(p[0], p[1], pr.polygons, tol)) { ++bad; break; }
    }
    std::printf("%s %-40s tri%-2d  %d elems, %d spill\n", bad ? "FAIL:" : "ok:  ", tag, order,
                mesh.element_count, bad);
    if (bad) ++g_failures;
    return bad;
}

m::Project rect_block() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3; s.gamma_unsat = 18; pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0; P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P); pr.has_water = false; return pr;
}

// L-shaped (non-convex) domain: a square with a bite taken out of the top-right.
m::Project l_shape() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3; s.gamma_unsat = 18; pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 10, 10, 0};
    P.y = {0, 0, 6, 6, 12, 12};
    P.edge_bc.assign(6, (int)m::BCType::Free);
    P.edge_bc[0] = (int)m::BCType::FullyFixed;
    pr.polygons.push_back(P); pr.has_water = false; return pr;
}

void add_plate(m::Project& pr, double x1, double y1, double x2, double y2) {
    m::PlateMaterial pm; pm.EA = 5e6; pm.EI = 8.5e3; pr.plates.push_back(pm);
    m::StructElement e; e.kind = m::StructKind::Plate; e.material = 0; e.x1 = x1; e.y1 = y1; e.x2 = x2; e.y2 = y2;
    pr.structs.push_back(e);
}
void add_geogrid(m::Project& pr, double x1, double y1, double x2, double y2) {
    m::GeogridMaterial gm; gm.EA = 1e4; pr.geogrids.push_back(gm);
    m::StructElement e; e.kind = m::StructKind::Geogrid; e.material = 0; e.x1 = x1; e.y1 = y1; e.x2 = x2; e.y2 = y2;
    pr.structs.push_back(e);
}
void add_dist_load(m::Project& pr, double x1, double y1, double x2, double y2) {
    m::Load L; L.kind = m::LoadKind::Distributed; L.x1 = x1; L.y1 = y1; L.x2 = x2; L.y2 = y2;
    L.qx1 = 0; L.qy1 = -50; L.qx2 = 0; L.qy2 = -50; pr.loads.push_back(L);
}

}  // namespace

int main() {
    std::printf("Mesh-domain integrity: no element may spill outside the soil polygon\n\n");
    for (int order : {6, 15}) {
        { auto pr = rect_block();                                         spill_count(pr, order, "rectangle (plain)"); }
        { auto pr = rect_block(); add_dist_load(pr, 6, 10, 14, 10);       spill_count(pr, order, "rect + surface distributed load"); }
        { auto pr = rect_block(); add_dist_load(pr, 5, 7, 15, 7);         spill_count(pr, order, "rect + INTERNAL distributed load"); }
        { auto pr = rect_block(); add_plate(pr, 10, 2, 10, 10);           spill_count(pr, order, "rect + vertical plate"); }
        { auto pr = rect_block(); add_plate(pr, 10, 2, 10, 12);           spill_count(pr, order, "rect + plate poking out the top"); }
        { auto pr = rect_block(); add_plate(pr, 4, 3, 16, 8);             spill_count(pr, order, "rect + inclined plate"); }
        { auto pr = rect_block(); add_geogrid(pr, 2, 5, 18, 5);           spill_count(pr, order, "rect + horizontal geogrid"); }
        { auto pr = rect_block(); add_plate(pr, 10, 2, 10, 10); add_geogrid(pr, 2, 5, 18, 5); add_dist_load(pr, 6, 10, 14, 10);
                                                                          spill_count(pr, order, "rect + plate + geogrid + load"); }
        { auto pr = l_shape();                                            spill_count(pr, order, "L-shape (non-convex)"); }
        { auto pr = l_shape(); add_plate(pr, 5, 0, 5, 12);                spill_count(pr, order, "L-shape + vertical plate"); }
        { auto pr = l_shape(); add_dist_load(pr, 0, 12, 10, 12);          spill_count(pr, order, "L-shape + surface load"); }
    }
    if (g_failures == 0) { std::printf("\nOK: every element lies inside the domain (no spill)\n"); return 0; }
    std::fprintf(stderr, "\n%d configuration(s) spilled outside the domain\n", g_failures);
    return 1;
}
