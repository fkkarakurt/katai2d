// Regression + robustness: a structural line (plate/geogrid) must never delete or wreck the mesh,
// for ANY placement -- fully interior, ending on the boundary, crossing the boundary, or sticking
// out of the soil. The meshed area must equal the part of the block the line cannot remove.
#include <katai/jobs/mesh_builder.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

static int fails = 0;
static void chk(bool ok, const char* w) { std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", w); if (!ok) ++fails; }

static double meshed_area(const katai::mesh::Mesh& mesh) {
    double a = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int p = mesh.node_of(e, 0), q = mesh.node_of(e, 1), r = mesh.node_of(e, 2);
        a += 0.5 * std::fabs((mesh.x[q]-mesh.x[p])*(mesh.y[r]-mesh.y[p]) -
                             (mesh.x[r]-mesh.x[p])*(mesh.y[q]-mesh.y[p]));
    }
    return a;
}
static m::Project block() {
    m::Project pr;
    m::Material s; s.E = 1e4; s.nu = 0.3; pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    pr.polygons.push_back(P); pr.has_water = false;
    return pr;
}
static m::Project with_plate(double x1, double y1, double x2, double y2) {
    m::Project pr = block();
    m::PlateMaterial pm; pr.plates.push_back(pm);
    m::StructElement s; s.kind = m::StructKind::Plate; s.material = 0;
    s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2; pr.structs.push_back(s);
    return pr;
}

static void case_(const char* tag, const m::Project& pr, int order) {
    const double A = 200.0;
    const auto R = katai::app::mesh_from_project(pr, 0.5 * 2.0 * 2.0, order);
    const double a = R.ok ? meshed_area(R.mesh) : -1;
    std::printf("  [%s order %2d] ok=%d area=%.2f/%.2f nodes=%d elems=%d\n",
                tag, order, (int)R.ok, a, A, R.mesh.node_count, R.mesh.element_count);
    chk(R.ok && std::fabs(a - A) < 0.05 * A, tag);
}

int main() {
    for (int order : {6, 15}) {
        case_("interior vertical",   with_plate(10, 2, 10, 8), order);
        case_("ends on top edge",    with_plate(10, 2, 10, 10), order);
        case_("crosses right edge",  with_plate(12, 5, 25, 5), order);   // extends outside x=20
        case_("spans + both outside",with_plate(-5, 5, 25, 5), order);   // pokes out both ends
        case_("inclined corner",     with_plate(2, 1, 18, 9), order);
        case_("near boundary",       with_plate(0.3, 5, 19.7, 5), order);
    }
    std::printf(fails ? "\n%d FAIL\n" : "\nOK: structural lines never wreck the mesh\n", fails);
    return fails ? 1 : 0;
}
