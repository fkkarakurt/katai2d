// Embedded barrier wall integration (P2.4): split_mesh_at_wall + plate on INDEPENDENT
// (extra) DOFs + an interface on each side (build_embedded_wall). Wiring check: with the
// soil DEACTIVATED, the wall (plate on its own extra DOFs, fixed at the toe) is a free
// Timoshenko cantilever under a horizontal tip load -> delta = PL^3/(3EI) + PL/(kGA').
// This isolates and validates the plate-on-extra-DOFs path + the wall-builder DOF
// allocation (independent of soil; the earth-pressure coupling is exercised separately).
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/embedded_wall.hpp>
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::SeamPair;
using katai::core::WallBuild;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(b)); }
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(k);
    return s->solve(r);
}

void test_wall_cantilever_limit() {
    constexpr double W = 4.0, L = 8.0, xw = 2.0, P = 50.0;
    const double Ep = 3.0e7, d = 0.4, nu = 0.15;     // stiff wall
    const double EI = Ep * d * d * d / 12.0;
    const double kGA = (5.0 / 6.0) * (Ep * d) / (2.0 * (1.0 + nu));
    const double delta_cant = P * L * L * L / (3.0 * EI) + P * L / kGA;

    RectangularDomain domain{0.0, 0.0, W, L, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 16);   // xw=2 a column line
    // Split the wall line above the base (toe at y=0, base row shared/fixed).
    std::vector<SeamPair> seam = katai::core::split_mesh_at_wall(mesh, xw, 0.0, L);
    // Toe node = the shared base node on the wall line (x=xw, y=0).
    int toe = -1;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - xw) < 1e-9 && std::fabs(mesh.y[n]) < 1e-9) toe = n;

    DofMap dofs(mesh.node_count, 2);
    katai::core::plate::PlateProps pp;
    pp.EA = Ep * d; pp.EI = EI; pp.nu = nu;
    katai::core::iface::InterfaceProps ip;  // (interfaces unused here -- soil deactivated)
    ip.kn = 1.0e2; ip.ks = 1.0e2; ip.c_i = 1.0e9; ip.phi_i = 0.0;
    const WallBuild wall = katai::core::build_embedded_wall(mesh, seam, toe, dofs, pp, ip);

    // Soil DEACTIVATED -> isolate the wall (plate on extra DOFs). Orphan soil nodes fixed.
    const std::vector<char> active(mesh.element_count, 0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    katai::core::fix_inactive_nodes(mesh, active, dofs);  // all soil nodes -> fixed (incl. toe)
    dofs.fix(wall.dof_phi.front());   // toe rotation fixed (cantilever fixed base)
    dofs.finalize();

    // Horizontal tip load at the wall top (last position) -- the wall's own trans-x DOF.
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const int top_x = dofs.equation(wall.dof_x.back());
    f(top_x) = P;

    const MaterialModel soft{MaterialType::LinearElastic, 1.0, 0.3};  // (unused; soil inactive)
    const auto r = katai::core::solve_nonlinear(
        mesh, dofs, {soft}, f, solve_spd, {1, 50, 1e-8}, {}, active,
        katai::core::Structures{wall.plates, {}, {}, {}});  // plates only (free cantilever)
    check(r.converged, "embedded-wall solve converged");

    for (size_t i = 0; i < wall.y.size(); ++i)
        std::printf("    pos y=%.2f ux=%.4e (cantilever ux=%.4e)\n", wall.y[i],
                    r.displacement[wall.dof_x[i]],
                    P * wall.y[i] * wall.y[i] * (3 * L - wall.y[i]) / (6.0 * EI));
    const double tip = std::fabs(r.displacement[wall.dof_x.back()]);
    std::printf("  wall cantilever (soft soil): tip=%.5e  free-cantilever=%.5e  (%.1f%%)  npos=%zu\n",
                tip, delta_cant, 100.0 * std::fabs(tip - delta_cant) / delta_cant, wall.y.size());
    check(close(tip, delta_cant, 0.08),
          "soft-soil limit = free Timoshenko cantilever PL^3/3EI + PL/kGA'");
}

} // namespace

int main() {
    test_wall_cantilever_limit();
    if (g_failures == 0) {
        std::printf("OK: embedded barrier wall (split + plate-extra-DOF + 2 interfaces) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
