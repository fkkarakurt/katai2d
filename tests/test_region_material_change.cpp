// A soil region that changes its material in a phase keeps its stresses; the new material takes
// only the load that comes after.
//
// A one-dimensional column (rollers on the sides, a fixed base, nu = 0, weightless) is loaded at
// its top in two steps. The first, q1, is carried by the material it was drawn with, E1: the top
// settles q1 H / E1. A later phase gives the region a stiffer material, E2, and adds q2: the
// stresses the column carries are kept, so only q2 meets the new stiffness and the top settles a
// further q2 H / E2. Recomputed with the new material over the whole history it would settle
// (q1 + q2) H / E2 - q1 H / E1 in that phase instead; and the vertical stress is q1 + q2 in both.
//
// verify: KV-CST-022
//   oracle:   closed_form
//   source:   one-dimensional (oedometric, nu = 0) linear elasticity with a change of material between two load increments: the stress is a state carried across the change, the stiffness a property of each increment; KATAI 2D input contract, docs/k2d-format.md, phases[i].pmat
//   locator:  a 1 x 10 m weightless linear elastic column (E1 = 1e4 kPa), q1 = 50 kPa in phase 1; phase 2 gives the region E2 = 4e4 kPa and adds q2 = 50 kPa
//   quantity: settlement of the top in each phase [m]; vertical effective stress at mid-height [kPa]
//   expected: 0.05 m in phase 1, a further 0.0125 m in phase 2 (a recomputation over the history would heave it 0.025 m); sigma'_yy = -100 kPa at the end
//   band:     1e-6 relative (the displacement field is linear, which 6-noded elements hold exactly)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
}  // namespace

int main() {
    constexpr double W = 1.0, H = 10.0, E1 = 1.0e4, E2 = 4.0e4, q1 = 50.0, q2 = 50.0;
    m::Project pr;
    m::Material a; a.name = "Soft"; a.model = m::SoilModel::LinearElastic;
    a.E = E1; a.nu = 0.0; a.gamma_unsat = 0.0; a.gamma_sat = 0.0;
    m::Material b = a; b.name = "Improved"; b.E = E2;
    pr.materials = {a, b};
    m::SoilPolygon P; P.name = "Column"; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.x_min = 0; pr.x_max = W; pr.y_min = 0; pr.y_max = H;
    pr.has_water = false;
    m::Load L1; L1.kind = m::LoadKind::Distributed; L1.name = "q1";
    L1.x1 = 0; L1.y1 = H; L1.x2 = W; L1.y2 = H; L1.qx1 = L1.qx2 = 0; L1.qy1 = L1.qy2 = -q1;
    m::Load L2 = L1; L2.name = "q2"; L2.qy1 = L2.qy2 = -q2;
    pr.loads = {L1, L2};
    pr.initial.load_active = {0, 0};
    m::Phase first; first.name = "q1 on E1"; first.load_active = {1, 0};
    m::Phase second = first; second.name = "E2, then q2"; second.load_active = {1, 1}; second.poly_material = {1};
    pr.phases = {first, second};

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "column meshed");
    if (!M.ok) return 1;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[1].ok && res[2].ok, "the two loading phases ran");
    if (res.size() != 3 || !res[2].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str());
        return 1;
    }
    const auto top_settlement = [&](const katai::app::SolveResult& r) {
        double s = 0.0; int n = 0;
        for (int i = 0; i < r.mesh.node_count; ++i)
            if (std::fabs(r.mesh.y[i] - H) < 1e-9) { s += -r.disp[2 * i + 1]; ++n; }
        return n ? s / n : 0.0;
    };
    const double s1 = top_settlement(res[1]), s2 = top_settlement(res[2]);
    double syy = 0.0; int ns = 0;
    for (int i = 0; i < res[2].mesh.node_count; ++i)
        if (std::fabs(res[2].mesh.y[i] - 0.5 * H) < 1e-9) { syy += res[2].stress.stress[i](1); ++ns; }
    syy = ns ? syy / ns : 0.0;
    const double w1 = q1 * H / E1, w2 = q2 * H / E2, recomputed = (q1 + q2) * H / E2 - q1 * H / E1;
    std::printf("  phase 1: the top settles %.8f m (closed form %.8f)\n", s1, w1);
    std::printf("  phase 2: a further %.8f m (kept stresses: %.8f; recomputed: %.8f)\n", s2, w2, recomputed);
    std::printf("  vertical stress at mid-height %.6f kPa (closed form %.6f)\n", syy, -(q1 + q2));
    check(std::fabs(s1 - w1) < 1e-6 * w1, "the first load meets the material the region was drawn with");
    check(std::fabs(s2 - w2) < 1e-6 * w2, "after the change only the new load meets the new material");
    check(std::fabs(syy + (q1 + q2)) < 1e-6 * (q1 + q2), "and the stress carried across the change is kept");
    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all region material change checks passed\n");
    return 0;
}
