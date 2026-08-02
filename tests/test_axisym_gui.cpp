// Axisymmetric analysis through the GUI compute path (build_problem with project.axisymmetric=true).
// PLAXIS's 2nd analysis mode -- circular footings, piles, shafts, cylinders. The core axisym element
// is validated directly (test_axisym_*, Lame %1e-5, cylinder collapse %0.3); here we validate the
// INTEGRATED GUI path: r-weighted gravity / internal force / traction + Kinematics::Axisymmetric.
//
//   (1) K0 procedure, confined cylinder: undisturbed geostatic state -> max|u|=0 and the recovered
//       effective stress is sigma_zz=-gamma(H-z), sigma_rr=K0 sigma_zz (axisym K0 baseline = r-weighted
//       assemble_axisym_internal_force balances the seeded state on any mesh).
//   (2) Gravity loading, confined cylinder: 1-D oedometric settlement -gamma H^2/(2 E_oed) (eps_rr=
//       eps_theta=0). Exercises the r-weighted body force assemble_axisym_gravity.
//   (3) Free-sided cylinder under a top surcharge q: UNIAXIAL stress (sigma_rr=sigma_theta=0), so the
//       settlement is qH/E and the RADIAL expansion is u_r(r)=nu q r/E -- the hoop strain eps_theta=
//       u_r/r is the defining axisymmetric feature, distinguishing it from plane strain. Exercises
//       assemble_axisym_traction_varying.
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

// Rectangular r-z domain [0,R]x[0,H], axisymmetric. edge_bc order: {bottom, right, top, left(axis)}.
m::Project cylinder(double R, double H, double E, double nu, double gamma, double phi,
                    m::BCType bottom, m::BCType right) {
    m::Project pr;
    pr.axisymmetric = true;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = nu; s.gamma_unsat = gamma; s.phi = phi; s.c = 1.0; s.k0_auto = true;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, R, R, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)bottom, (int)right, (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

void test_axisym_k0_cylinder() {
    const double R = 5.0, H = 10.0, E = 1.0e4, nu = 0.3, gamma = 18.0, phi = 30.0;
    const double K0 = 1.0 - std::sin(phi * 3.14159265358979 / 180.0);
    m::Project pr = cylinder(R, H, E, nu, gamma, phi, m::BCType::FullyFixed, m::BCType::HorizontallyFixed);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "axisym K0 cylinder meshed");
    const auto Rk = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(Rk.ok, "axisym K0 solve ok");
    double max_sv = 0.0, max_sh = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        const double y = M.mesh.y[n];
        const double sv_exact = -gamma * (H - y);
        max_sv = std::fmax(max_sv, std::fabs(Rk.stress.stress[n](1) - sv_exact));
        max_sh = std::fmax(max_sh, std::fabs(Rk.stress.stress[n](0) - K0 * sv_exact));
    }
    std::printf("  (1) axisym K0: max|u|=%.2e  max|sig_zz err|=%.3e  max|sig_rr err|=%.3e\n",
                Rk.max_disp, max_sv, max_sh);
    check(Rk.max_disp < 1e-6, "axisym K0: undisturbed cylinder does not move (geostatic equilibrium)");
    check(max_sv < 0.5 && max_sh < 0.5, "axisym K0: recovered sigma_zz=-gamma(H-z), sigma_rr=K0 sigma_zz");
}

void test_axisym_gravity_oedometer() {
    const double R = 5.0, H = 10.0, E = 1.0e4, nu = 0.3, gamma = 18.0;
    m::Project pr = cylinder(R, H, E, nu, gamma, 30.0, m::BCType::FullyFixed, m::BCType::HorizontallyFixed);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto R2 = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    check(R2.ok, "axisym gravity-loading solve ok");
    const double Eoed = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double u_exact = -gamma * H * H / (2.0 * Eoed);
    double u_top = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n)
        if (M.mesh.y[n] > H - 1e-6) u_top = std::fmin(u_top, R2.disp[n * 2 + 1]);
    std::printf("  (2) axisym gravity oedometer: u_top=%.6f  exact -gamma H^2/2E_oed=%.6f\n", u_top, u_exact);
    check(std::fabs(u_top - u_exact) < 0.02 * std::fabs(u_exact),
          "axisym confined cylinder settles -gamma H^2/(2 E_oed) (r-weighted gravity)");
}

void test_axisym_uniaxial_radial() {
    // Free-sided solid cylinder, top surcharge q, weightless -> uniaxial: settle qH/E, u_r(R)=nu q R/E.
    const double R = 4.0, H = 8.0, E = 2.0e4, nu = 0.3, q = 50.0;
    m::Project pr = cylinder(R, H, E, nu, 0.0, 30.0, m::BCType::VerticallyFixed, m::BCType::Free);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 0; L.y1 = H; L.x2 = R; L.y2 = H; L.qx1 = 0; L.qy1 = -q; L.qx2 = 0; L.qy2 = -q;  // top pressure
    pr.loads.push_back(L);
    const auto M = katai::app::mesh_from_project(pr, 0.25, 6);
    const auto Rr = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    check(Rr.ok, "axisym uniaxial solve ok");
    const double u_z_exact = -q * H / E;          // free sides -> uniaxial (NOT E_oed)
    const double u_r_exact = nu * q * R / E;       // radial expansion u_r(R) = nu q R / E (hoop strain)
    double u_top = 0.0, u_r_outer = 0.0; double bestd = 1e300;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        if (M.mesh.y[n] > H - 1e-6) u_top = std::fmin(u_top, Rr.disp[n * 2 + 1]);
        const double d = std::hypot(M.mesh.x[n] - R, M.mesh.y[n] - H / 2);  // outer wall, mid-height
        if (d < bestd) { bestd = d; u_r_outer = Rr.disp[n * 2 + 0]; }
    }
    std::printf("  (3) axisym uniaxial: u_z=%.6f (exact %.6f)  u_r(R)=%.6f (exact %.6f)\n",
                u_top, u_z_exact, u_r_outer, u_r_exact);
    check(std::fabs(u_top - u_z_exact) < 0.03 * std::fabs(u_z_exact), "axisym uniaxial settlement qH/E");
    check(std::fabs(u_r_outer - u_r_exact) < 0.05 * std::fabs(u_r_exact),
          "axisym RADIAL expansion u_r(R)=nu q R/E (hoop strain -- distinguishes axisymmetry)");
}

}  // namespace

int main() {
    std::printf("Axisymmetric analysis through the GUI compute path (build_problem)\n");
    test_axisym_k0_cylinder();
    test_axisym_gravity_oedometer();
    test_axisym_uniaxial_radial();
    if (g_failures == 0) {
        std::printf("\nOK: axisymmetric GUI path verified (K0, gravity oedometer, uniaxial radial)\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
