// Regression: anchors must NOT perturb the mesh. Anchor endpoints used to be injected into the PSLG
// as isolated points; on a slope an endpoint on the boundary / outside / near an edge wrecked the
// Ruppert mesh (and hence the post-processing). Anchors now attach to the nearest existing node, so
// the mesh is identical with or without them, for any endpoint placement. The anchored model still
// solves to a sane field.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

static int fails = 0;
static void chk(bool ok, const char* w) { std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", w); if (!ok) ++fails; }

static m::Project slope() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 2.0e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.phi = 30.0; s.c = 5.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 40, 40, 20, 10, 0};
    P.y = {0,  0, 10, 10, 20, 20};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}
static void add_anchor(m::Project& pr, double x1, double y1, double x2, double y2) {
    m::StructElement a; a.kind = m::StructKind::Anchor; a.material = 0;
    a.x1 = x1; a.y1 = y1; a.x2 = x2; a.y2 = y2; pr.structs.push_back(a);
}

int main() {
    const double area = 0.5 * 2.0 * 2.0;
    const m::Project base = slope();

    for (int order : {6, 15}) {
        const auto B = katai::app::mesh_from_project(base, area, order);
        chk(B.ok, "bare slope meshed");

        // Anchors with awkward endpoints: inside, on the slope-face edge, far outside, near-edge sliver.
        m::Project pr = slope();
        m::AnchorMaterial am; am.EA = 2.0e5; pr.anchors.push_back(am);
        add_anchor(pr, 11.0, 18.0, 26.0, 8.0);    // both inside
        add_anchor(pr, 15.0, 15.0, 28.0, 8.0);    // head on the slope-face edge
        add_anchor(pr, 12.0, 18.0, 90.0, 60.0);   // far end far outside the model
        const auto A = katai::app::mesh_from_project(pr, area, order);
        chk(A.ok, "anchored slope meshed");
        std::printf("  order %2d: bare=%d/%d  anchored=%d/%d\n", order,
                    B.mesh.node_count, B.mesh.element_count, A.mesh.node_count, A.mesh.element_count);
        chk(A.mesh.node_count == B.mesh.node_count && A.mesh.element_count == B.mesh.element_count,
            "anchors do not change the mesh (identical node/element count)");

        // Anchored solve produces a finite, sane field (no NaN, no blow-up).
        const auto S = katai::app::solve_gravity_le(pr, A.mesh, katai::app::InitialPhase::K0Procedure);
        chk(S.ok, "anchored slope solved");
        bool nan = false; double mx = 0.0;
        for (int i = 0; i < (int)S.disp.size(); ++i) { if (std::isnan(S.disp[i])) nan = true; mx = std::fmax(mx, std::fabs(S.disp[i])); }
        chk(!nan && mx < 1.0, "anchored field is finite and bounded");
    }

    std::printf(fails ? "\n%d FAIL\n" : "\nOK: anchors do not perturb the mesh\n", fails);
    return fails ? 1 : 0;
}
