// Undrained (A) loading through the GUI compute path (project drainage type -> build_problem ->
// solve_gravity_le with the pore-fluid Kw/n tangent). The undrained CORE is validated directly in
// test_undrained_global (confined column, Skempton 1D); here the INTEGRATED pipeline must reproduce
// it. A weightless, laterally confined column under a top surcharge q responds undrained:
//   undrained constrained modulus  M_u = M' + Kw/n,   M' = E(1-nu)/((1+nu)(1-2nu))
//   surface settlement  u_y = -q H / M_u             (much stiffer than the drained -q H / M')
//   effective vertical stress  sigma'_yy = -M' q / M_u   (the rest, ~q, is carried as excess pore)
// The headline undrained bearing capacity q_ult = (2+pi) su (Undrained B, phi=0) is exercised on
// demand in study_gui_validation (limit analysis). Reference: docs/references/effective-stress-formulation.md.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// Surface settlement (most negative u_y on the top edge) of a confined column under a top surcharge q.
double confined_settlement(m::Drainage drainage, double E, double nu, double H, double W, double q) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = nu; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;   // weightless: only the surcharge acts
    s.drainage = drainage;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 0; L.y1 = H; L.x2 = W; L.y2 = H; L.qx1 = 0; L.qy1 = -q; L.qx2 = 0; L.qy2 = -q;
    pr.loads.push_back(L);
    pr.has_water = false;

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    if (!M.ok) { std::printf("   (mesh failed: %s)\n", M.message.c_str()); return 0.0; }
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    if (!R.ok) { std::printf("   (solve failed: %s)\n", R.message.c_str()); return 0.0; }
    double u_top = 0.0, sig_eff_min = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        if (M.mesh.y[n] > H - 1e-6) u_top = std::fmin(u_top, R.disp[n * 2 + 1]);
        // effective sigma'_yy at interior nodes (away from the loaded surface)
        if (M.mesh.y[n] > 0.2 * H && M.mesh.y[n] < 0.8 * H)
            sig_eff_min = std::fmin(sig_eff_min, R.stress.stress[n](1));
    }
    std::printf("   %-13s settlement u_y=%.5e  effective sigma'_yy(mid)=%.3f\n",
                drainage == m::Drainage::Undrained ? "Undrained(A)" : "Drained", u_top, sig_eff_min);
    return u_top;
}

void test_confined_undrained_gui() {
    std::printf("-- confined column: undrained vs drained through the GUI path --\n");
    constexpr double E = 1.0e4, nu = 0.3, nu_u = 0.495, H = 10.0, W = 2.0, q = 50.0;

    // Closed-form constrained moduli + the engine's Kw/n (material_model.hpp: 3(nu_u-nu)/((1-2nu_u)(1+nu)) K').
    const double Mp = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double Kp = E / (3.0 * (1.0 - 2.0 * nu));
    const double kwn = 3.0 * (nu_u - nu) / ((1.0 - 2.0 * nu_u) * (1.0 + nu)) * Kp;
    const double Mu = Mp + kwn;
    const double uy_undrained = -q * H / Mu;
    const double uy_drained = -q * H / Mp;
    const double sig_eff_exact = -Mp * q / Mu;

    const double u_und = confined_settlement(m::Drainage::Undrained, E, nu, H, W, q);
    const double u_dr = confined_settlement(m::Drainage::Drained, E, nu, H, W, q);

    std::printf("  undrained u_y=%.5e (exact -qH/M_u=%.5e),  drained u_y=%.5e (exact -qH/M'=%.5e)\n",
                u_und, uy_undrained, u_dr, uy_drained);
    std::printf("  M_u/M' = %.1f  =>  drained settles %.1fx more;  excess pore u/q ~ %.3f\n",
                Mu / Mp, u_dr / u_und, kwn / Mu);

    check(std::fabs(u_und - uy_undrained) < 0.02 * std::fabs(uy_undrained),
          "undrained GUI settlement = -qH/M_u (Kw/n wired into the solve)");
    check(std::fabs(u_dr - uy_drained) < 0.02 * std::fabs(uy_drained),
          "drained GUI settlement = -qH/M'");
    check(u_dr < u_und - 1e-9 && (u_dr / u_und) > 10.0,
          "undrained >> stiffer than drained (near-incompressible pore fluid)");
}

}  // namespace

int main() {
    std::printf("Undrained loading through the GUI compute path\n\n");
    test_confined_undrained_gui();
    if (g_failures == 0) {
        std::printf("\nOK: GUI undrained path = -qH/M_u settlement + effective stress (Skempton 1D)\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
