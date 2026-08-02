// test_solve — GUI model -> mesh -> linear-elastic self-weight solve (build_problem.hpp).
//
// Analytic check: a laterally-confined soil column (base fully fixed, sides on
// normal rollers, top free) under self-weight is 1-D oedometric:
//   sigma_yy(y) = -gamma (H - y),   E_oed = E(1-v) / ((1+v)(1-2v)),
//   u_y,top = -gamma H^2 / (2 E_oed).
// tri6 represents the quadratic u_y(y) exactly, so the FE result should match
// closely. Also checks sigma_yy at the base and that u_x ~ 0 (truly 1-D).

#include <cmath>
#include <cstdio>
#include <vector>

#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>

namespace m = katai::model;

static int failures = 0;
#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) { std::printf("FAIL: %s\n", msg); ++failures; }             \
        else std::printf("ok:   %s\n", msg);                                     \
    } while (0)

int main() {
    constexpr double W = 2.0, H = 10.0, E = 10000.0, nu = 0.3, gamma = 20.0;

    m::Project pr;
    m::Material mat; mat.model = m::SoilModel::LinearElastic;   // LE oedometric / K0 analytic checks
    mat.E = E; mat.nu = nu; mat.gamma_unsat = gamma; pr.materials.push_back(mat);
    m::SoilPolygon P;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H}; P.material = 0;
    // edges: 0 bottom, 1 right, 2 top, 3 left
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);

    const auto M = katai::app::mesh_from_project(pr, /*max_area*/ 0.5, /*tri6*/ 6);
    CHECK(M.ok, "column meshed");
    // Gravity loading: self-weight from a stress-free start -> 1-D oedometric settlement.
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::GravityLoading);
    CHECK(R.ok, "column solved");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return 1; }

    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double u_top_exact = -gamma * H * H / (2.0 * Eoed);

    // settlement at the top (most negative u_y among top nodes)
    double u_top = 0.0, max_ux = 0.0, min_syy = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        max_ux = std::fmax(max_ux, std::fabs(R.disp[n * 2]));
        if (M.mesh.y[n] > H - 1e-6) u_top = std::fmin(u_top, R.disp[n * 2 + 1]);
        if (M.mesh.y[n] < 1e-6) min_syy = std::fmin(min_syy, R.stress.stress[n].y());
    }
    std::printf("      u_top  FE = %.6f   exact = %.6f   (E_oed = %.1f)\n", u_top, u_top_exact, Eoed);
    std::printf("      sigma_yy base FE = %.3f   exact = %.3f\n", min_syy, -gamma * H);
    std::printf("      max |u_x| = %.3e\n", max_ux);

    CHECK(std::fabs(u_top - u_top_exact) < 0.02 * std::fabs(u_top_exact), "top settlement matches 1-D oedometric (<2%)");
    CHECK(std::fabs(min_syy - (-gamma * H)) < 0.05 * gamma * H, "base sigma_yy matches -gamma*H (<5%)");
    CHECK(max_ux < 1e-3, "u_x ~ 0 (confined, 1-D)");

    // K0 procedure: the same column starts in geostatic equilibrium -> self-weight gives ~zero
    // displacement, and the recovered effective stress is the K0 field (base sigma_yy = -gamma*H).
    const auto RK = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    CHECK(RK.ok, "K0 column solved");
    double k0_max_u = 0.0, k0_base_syy = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        k0_max_u = std::fmax(k0_max_u, std::hypot(RK.disp[n * 2], RK.disp[n * 2 + 1]));
        if (M.mesh.y[n] < 1e-6) k0_base_syy = std::fmin(k0_base_syy, RK.stress.stress[n].y());
    }
    std::printf("      K0: max|u| = %.3e   base sigma_yy = %.3f (exact %.3f)\n",
                k0_max_u, k0_base_syy, -gamma * H);
    CHECK(k0_max_u < 1e-6, "K0 procedure: undisturbed self-weight gives ~zero displacement");
    CHECK(std::fabs(k0_base_syy - (-gamma * H)) < 0.05 * gamma * H, "K0 base sigma_yy matches -gamma*H");

    // ---- Structural plate coupling: a stiff buried plate spreads a point load ---------------------
    // 10x5 soil block, central surface point load. A stiff horizontal plate buried under the load
    // distributes the pressure, so the settlement at the load point is smaller than without it.
    auto settlement_at_load = [](bool with_plate) -> double {
        m::Project pr;
        m::Material mt; mt.model = m::SoilModel::LinearElastic;
        mt.E = 5000.0; mt.nu = 0.3; mt.gamma_unsat = 0.0; pr.materials.push_back(mt);
        m::SoilPolygon Q; Q.x = {0, 10, 10, 0}; Q.y = {0, 0, 5, 5}; Q.material = 0;
        Q.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                     (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
        pr.polygons.push_back(Q);
        if (with_plate) {
            m::PlateMaterial pm; pm.EA = 1e7; pm.EI = 1e6; pm.nu = 0.15; pr.plates.push_back(pm);
            m::StructElement s; s.kind = m::StructKind::Plate; s.x1 = 2; s.y1 = 4; s.x2 = 8; s.y2 = 4; s.material = 0;
            pr.structs.push_back(s);
        }
        m::Load L; L.kind = m::LoadKind::Point; L.x1 = 5; L.y1 = 5; L.qx1 = 0; L.qy1 = -100; pr.loads.push_back(L);
        const auto MM = katai::app::mesh_from_project(pr, 0.25, 6);
        const auto RR = katai::app::solve_gravity_le(pr, MM.mesh);
        if (!RR.ok) { std::printf("  (plate-case solve failed: %s)\n", RR.message.c_str()); return 0.0; }
        int nn = -1; double bd = 1e300;
        for (int n = 0; n < MM.mesh.node_count; ++n) {
            const double d = std::hypot(MM.mesh.x[n] - 5.0, MM.mesh.y[n] - 5.0);
            if (d < bd) { bd = d; nn = n; }
        }
        return -RR.disp[nn * 2 + 1];
    };
    const double d_noplate = settlement_at_load(false);
    const double d_plate = settlement_at_load(true);
    std::printf("      settlement at load: no plate = %.5e, with plate = %.5e\n", d_noplate, d_plate);
    CHECK(d_noplate > 0.0 && d_plate > 0.0, "plate case: both solves produced settlement");
    CHECK(d_plate < d_noplate, "stiff buried plate reduces the load-point settlement (coupling works)");

    // ---- Distributed (surcharge) load: full-top uniform pressure on a confined column ----------
    // No gravity; a uniform surface pressure q over the whole top gives uniform sigma_yy = -q and
    // 1-D oedometric settlement u_top = -q H / E_oed. This exercises the distributed-load path end
    // to end (mesh conforms to the load line -> consistent nodal forces in build_problem).
    {
        constexpr double Wc = 2.0, Hc = 10.0, Ec = 12000.0, nuc = 0.3, q = 60.0;
        m::Project pr2;
        m::Material mt; mt.model = m::SoilModel::LinearElastic;
        mt.E = Ec; mt.nu = nuc; mt.gamma_unsat = 0.0; pr2.materials.push_back(mt);
        m::SoilPolygon Q; Q.x = {0, Wc, Wc, 0}; Q.y = {0, 0, Hc, Hc}; Q.material = 0;
        Q.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                     (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
        pr2.polygons.push_back(Q);
        m::Load L; L.kind = m::LoadKind::Distributed;
        L.x1 = 0; L.y1 = Hc; L.x2 = Wc; L.y2 = Hc;     // along the top edge
        L.qx1 = 0; L.qy1 = -q; L.qx2 = 0; L.qy2 = -q;  // uniform downward surcharge
        pr2.loads.push_back(L);

        const auto M2 = katai::app::mesh_from_project(pr2, 0.5, 6);
        CHECK(M2.ok, "surcharge column meshed (load line conforms)");
        const auto R2 = katai::app::solve_gravity_le(pr2, M2.mesh, katai::app::InitialPhase::GravityLoading);
        CHECK(R2.ok, "surcharge column solved");
        const double Eoed2 = Ec * (1.0 - nuc) / ((1.0 + nuc) * (1.0 - 2.0 * nuc));
        const double u_exact = -q * Hc / Eoed2;
        double u_top2 = 0.0, max_ux2 = 0.0;
        for (int n = 0; n < M2.mesh.node_count; ++n) {
            max_ux2 = std::fmax(max_ux2, std::fabs(R2.disp[n * 2]));
            if (M2.mesh.y[n] > Hc - 1e-6) u_top2 = std::fmin(u_top2, R2.disp[n * 2 + 1]);
        }
        std::printf("      surcharge u_top FE = %.6f   exact -qH/E_oed = %.6f   max|u_x|=%.2e\n",
                    u_top2, u_exact, max_ux2);
        CHECK(std::fabs(u_top2 - u_exact) < 0.02 * std::fabs(u_exact),
              "distributed surcharge gives 1-D oedometric settlement -qH/E_oed (<2%)");
        CHECK(max_ux2 < 1e-3, "surcharge: u_x ~ 0 (confined, 1-D)");
    }

    std::printf(failures ? "\n%d CHECK(S) FAILED\n" : "\nAll checks passed.\n", failures);
    return failures ? 1 : 0;
}
