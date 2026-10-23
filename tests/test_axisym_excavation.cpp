// test_axisym_excavation -- excavation and fill in an axisymmetric model.
//
// A shaft, a cylindrical pit, the foundation of a tank dug and backfilled: until 2026-10-22 an
// axisymmetric model refused every change of the active regions, on the ground that the
// r-weighted assemblies took no element-activity mask. They did; two calls did not pass it.
//
// Two oracles, each sharing no code with the edit:
//   (1) closed form. A laterally confined column (the axis and a rolled outer radius) under its own
//       weight, K0 = nu/(1 - nu) so that the initial state is the elastic one, loses its top layer
//       of thickness h over the whole radius. The remaining ground of height H - h unloads one-
//       dimensionally by gamma h at every depth: the new surface heaves gamma h (H - h) / E_oed.
//   (2) path independence of linear elasticity. A cylindrical pit (radius a, depth d) is dug in a
//       phase. In a linear elastic soil the stress after the excavation depends on the final
//       geometry only, so it must equal the stress that gravity loading gives to the same mesh
//       with the pit inactive from the start -- a second construction of the answer that never
//       passes through an excavation.
//
// verify: KV-EXC-003
//   oracle:   closed_form
//   source:   one-dimensional elastic unloading of a laterally confined column, heave = gamma h (H - h) / E_oed with E_oed = E (1 - nu) / ((1 + nu)(1 - 2 nu)); and the path independence of a linear elastic body (the stress of the final geometry under its own weight), evaluated in this file by gravity loading the excavated geometry directly
//   locator:  the heave at the new surface on the axis and at the outer radius; the nodal stresses of the excavated phase against those of the directly loaded geometry at every node of the remaining ground
//   quantity: heave [m]; nodal stresses sigma_rr, sigma_zz, sigma_rz [kPa]
//   expected: the heave within 0.5% of the closed form, on 6- and 15-noded elements; the pit's stresses equal to the directly loaded geometry's within 1e-6 of the largest stress
//   band:     0.5% on the heave (nodal recovery on a 0.5 m mesh), 1e-6 relative on the identity (solver round-off)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

static m::Material elastic(double E, double nu, double gamma) {
    m::Material s; s.model = m::SoilModel::LinearElastic; s.E = E; s.nu = nu;
    s.gamma_unsat = s.gamma_sat = gamma; s.e_init = 0.5;
    s.k0_auto = false; s.k0 = nu / (1.0 - nu);
    return s;
}
static m::SoilPolygon rect(double x0, double y0, double x1, double y1, int mat, const char* name) {
    m::SoilPolygon P; P.name = name; P.material = mat;
    P.x = {x0, x1, x1, x0}; P.y = {y0, y0, y1, y1};
    P.edge_bc = {(int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free};
    return P;
}

int main() {
    std::printf("Excavation and fill in an axisymmetric model\n\n");
    constexpr double E = 20000.0, nu = 0.3, gamma = 18.0;
    const double Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

    // (1) the top layer of a confined column, over the whole radius
    for (int order : {6, 15}) {
        constexpr double R = 5.0, H = 10.0, h = 2.0;
        m::Project pr; pr.axisymmetric = true;
        pr.materials.push_back(elastic(E, nu, gamma));
        m::SoilPolygon low = rect(0, 0, R, H - h, 0, "Ground");
        low.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                       (int)m::BCType::HorizontallyFixed};
        m::SoilPolygon top = rect(0, H - h, R, H, 0, "Top layer");
        top.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                       (int)m::BCType::HorizontallyFixed};
        pr.polygons = {low, top};
        pr.has_water = false;
        m::Phase dig; dig.name = "Dig"; dig.type = m::PhaseType::Plastic; dig.poly_active = {1, 0};
        pr.phases.push_back(dig);
        const auto M = katai::app::mesh_from_project(pr, 0.5, order);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const std::string tag = std::to_string(order) + "-noded";
        check(M.ok && res.size() == 2 && res[1].ok, tag + ": the axisymmetric excavation solves");
        if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); continue; }
        const double want = gamma * h * (H - h) / Eoed;
        double axis = 0.0, rim = 0.0;
        for (int n = 0; n < res[1].mesh.node_count; ++n) {
            if (std::fabs(res[1].mesh.y[n] - (H - h)) > 1e-9) continue;
            if (std::fabs(res[1].mesh.x[n]) < 1e-9) axis = res[1].disp[2 * n + 1];
            if (std::fabs(res[1].mesh.x[n] - R) < 1e-9) rim = res[1].disp[2 * n + 1];
        }
        std::printf("  %s: heave on the axis %.6e m, at the rim %.6e m, closed form %.6e m\n", tag.c_str(), axis, rim, want);
        check(std::fabs(axis - want) / want < 0.005 && std::fabs(rim - want) / want < 0.005,
              tag + ": the heave is the one-dimensional closed form within 0.5%");
    }

    // (2) a cylindrical pit, against the same geometry loaded with the pit absent from the start
    {
        constexpr double R = 15.0, H = 12.0, a = 3.0, d = 4.0;
        auto build = [&](bool dig_in_phase) {
            m::Project pr; pr.axisymmetric = true;
            pr.materials.push_back(elastic(E, nu, gamma));
            // the ground outside the pit, the ground under it, and the pit
            m::SoilPolygon below = rect(0, 0, R, H - d, 0, "Below");
            below.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                             (int)m::BCType::HorizontallyFixed};
            m::SoilPolygon ring = rect(a, H - d, R, H, 0, "Ring");
            ring.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                            (int)m::BCType::Free};
            m::SoilPolygon pit = rect(0, H - d, a, H, 0, "Pit");
            pit.edge_bc = {(int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free,
                           (int)m::BCType::HorizontallyFixed};
            pr.polygons = {below, ring, pit};
            pr.has_water = false;
            if (dig_in_phase) {
                m::Phase dig; dig.name = "Dig"; dig.type = m::PhaseType::Plastic; dig.poly_active = {1, 1, 0};
                pr.phases.push_back(dig);
            } else {
                pr.initial.poly_active = {1, 1, 0};
            }
            return pr;
        };
        const m::Project staged = build(true), direct = build(false);
        const auto M = katai::app::mesh_from_project(staged, 0.75, 6);
        const auto rs = katai::app::solve_phases(staged, M.mesh, InitialPhase::GravityLoading);
        const auto rd = katai::app::solve_phases(direct, M.mesh, InitialPhase::GravityLoading);
        const bool ok = rs.size() == 2 && rs[1].ok && rd.size() == 1 && rd[0].ok;
        check(ok, "the pit dug in a phase and the geometry loaded without it both solve");
        if (!ok) {
            if (!rs.empty()) std::printf("  (%s)\n", rs.back().message.c_str());
            if (!rd.empty()) std::printf("  (%s)\n", rd.back().message.c_str());
        } else {
            const auto& A = rs[1];
            const auto& B = rd[0];
            double worst = 0.0, scale = 0.0;
            int nodes = 0;
            for (int n = 0; n < A.mesh.node_count; ++n) {
                const double r = A.mesh.x[n], z = A.mesh.y[n];
                if (r < a - 1e-9 && z > H - d + 1e-9) continue;   // inside the pit
                const auto& sa = A.stress.stress[(size_t)n];
                const auto& sb = B.stress.stress[(size_t)n];
                scale = std::max(scale, sb.cwiseAbs().maxCoeff());
                worst = std::max(worst, (sa - sb).cwiseAbs().maxCoeff());
                ++nodes;
            }
            std::printf("  pit: largest stress difference %.3e kPa over %d nodes (largest stress %.2f kPa)\n",
                        worst, nodes, scale);
            check(nodes > 100 && worst < 1e-6 * scale,
                  "the excavated pit's stress is the directly loaded geometry's (path independence)");
        }
    }

    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
