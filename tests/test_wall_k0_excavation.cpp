// Embedded barrier wall + K0 interface seeding + staged excavation (P2.4). The wall is a true
// barrier: the mesh is split along the wall line (split_mesh_at_wall) and the wall is a plate on
// independent extra DOFs connected to each soil side by an interface (build_embedded_wall). On a
// split mesh the K0 horizontal stress at the wall line is discontinuous, so without seeding the
// interface starts at zero normal stress and the wall "installs" with a large spurious movement
// that swamps the excavation effect. seed_interface_k0 sets sigma_n0 = K0*sigma_v on each
// interface so the wished-in-place wall is in equilibrium at u=0 (zero deformation).
//
//   (A) No excavation: the seeded wall does NOT move (K0 install equilibrium, exact).
//   (B) Excavation (one side): the staged release (baseline B = soil sigma0 + interface sigma_n0
//       at u=0, ramped to active gravity) makes the wall deflect toward the cut -- a genuine,
//       excavation-driven response, NOT the installation artifact.
//
// Reference: PLAXIS K0 procedure -> interface initial stress; interface-formulation.md §6.
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/embedded_wall.hpp>
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace katai::core;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// Staged-construction baseline B = f_int at u=0 = integral_active B^T sigma0 (soil) + interface
// sigma_n0 forces (tau=0 at u=0). target(lambda) = B + lambda*(grav_active - B): residual 0 at
// lambda=0 (equilibrium with the seeded state), full release at lambda=1.
Eigen::VectorXd baseline_force(const Mesh& mesh, const DofMap& dofs,
                               const std::vector<GaussState>& init,
                               const std::vector<char>& active,
                               const std::vector<InterfaceElement>& ifaces) {
    namespace tri6 = katai::core::tri6;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto gp = tri6::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active[e]) continue;
        tri6::NodeCoords X;
        for (int k = 0; k < 6; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        Eigen::Matrix<double, 12, 1> fe = Eigen::Matrix<double, 12, 1>::Zero();
        for (int g = 0; g < 3; ++g) {
            const auto sg = tri6::strain_displacement(X, gp[g].xi, gp[g].eta);
            fe.noalias() += (gp[g].weight * sg.det_jacobian) * sg.B.transpose() * init[e * 3 + g].stress;
        }
        for (int k = 0; k < 6; ++k)
            for (int c = 0; c < 2; ++c) {
                const int eq = dofs.equation(dofs.global_dof(mesh.node_of(e, k), c));
                if (eq >= 0) F(eq) += fe(2 * k + c);
            }
    }
    const auto ncp = iface::nc_points();
    for (const auto& ie : ifaces) {
        iface::NodeCoords Xe;
        for (int k = 0; k < 3; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
        for (int q = 0; q < 3; ++q) {
            const int nd = ncp[q].node;
            const auto fr = iface::edge_frame(Xe, ncp[q].xi);
            const double b[4] = {fr.s, -fr.c, -fr.s, fr.c};   // d(Du_n)/d(dof): [soil_x,soil_y,struct_x,struct_y]
            const double wJ = ncp[q].w * fr.J;
            const int idx[4] = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                dofs.equation(ie.struct_dof[2 * nd + 0]),
                                dofs.equation(ie.struct_dof[2 * nd + 1])};
            for (int i = 0; i < 4; ++i)
                if (idx[i] >= 0) F(idx[i]) += wJ * b[i] * ie.sigma_n0[q];
        }
    }
    return F;
}

struct Out { std::vector<double> y, ux; bool converged = false; };

Out run(bool excavate) {
    constexpr double W = 20.0, Htot = 12.0, xw = 10.0, y_exc = 6.0, gamma = 18.0;
    const double pi = std::acos(-1.0);
    const double phi = 30.0 * pi / 180.0, c = 2.0, psi = phi;
    const double K0 = 1.0 - std::sin(phi);
    constexpr double E = 2.0e4, nu = 0.3;
    const double Ep = 2.0e7, d = 0.4;
    const double EA = Ep * d, EI = Ep * d * d * d / 12.0;

    RectangularDomain domain{0.0, 0.0, W, Htot, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 10, 6);
    std::vector<SeamPair> seam = split_mesh_at_wall(mesh, xw, 0.0, Htot);
    int toe = -1;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - xw) < 1e-9 && std::fabs(mesh.y[n]) < 1e-9) toe = n;

    DofMap dofs(mesh.node_count, 2);
    plate::PlateProps pp; pp.EA = EA; pp.EI = EI; pp.nu = 0.15;
    iface::InterfaceProps ip; ip.kn = 1e5; ip.ks = 1e5;
    ip.c_i = 2.0;                                 // R_inter ~ 0.67 (PLAXIS strength reduction)
    ip.phi_i = std::atan(0.67 * std::tan(phi));
    WallBuild wall = build_embedded_wall(mesh, seam, toe, dofs, pp, ip);

    const K0Options k0o{Htot, gamma, K0};
    seed_interface_k0(wall, mesh, k0o);  // sigma_n0 = K0*sigma_v (wished-in-place equilibrium)

    std::vector<char> active(mesh.element_count, 1);
    for (int e = 0; e < mesh.element_count && excavate; ++e) {
        double xc = 0.0, yc = 0.0;
        for (int k = 0; k < 6; ++k) { xc += mesh.x[mesh.node_of(e, k)]; yc += mesh.y[mesh.node_of(e, k)]; }
        if (xc / 6.0 < xw && yc / 6.0 > y_exc) active[e] = 0;
    }

    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    fix_inactive_nodes(mesh, active, dofs);
    dofs.fix(wall.dof_phi.front());
    dofs.finalize();

    const std::vector<GaussState> init = compute_k0_initial_stress(mesh, k0o);
    Eigen::VectorXd grav = Eigen::VectorXd::Zero(dofs.equation_count());
    assemble_gravity(mesh, dofs, {gamma}, grav, active);

    const Eigen::VectorXd B = baseline_force(mesh, dofs, init, active, wall.interfaces);
    const Eigen::VectorXd f_ramp = grav - B;
    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, E, nu, c, phi, psi}};
    const auto r = solve_nonlinear(mesh, dofs, mm, f_ramp, solve_unsym, {100, 60, 1e-6}, init,
                                   active, Structures{wall.plates, {}, {}, wall.interfaces}, {}, B);

    Out o; o.converged = r.converged;
    for (size_t i = 0; i < wall.y.size(); ++i) {
        o.y.push_back(wall.y[i]);
        o.ux.push_back(r.displacement[wall.dof_x[i]]);
    }
    return o;
}

void test_wall_k0_excavation() {
    // (A) Wished-in-place wall, no excavation: K0 install equilibrium -> wall must not move.
    const Out a = run(false);
    check(a.converged, "K0 install (no excavation) converged");
    double max_a = 0.0;
    for (double u : a.ux) max_a = std::fmax(max_a, std::fabs(u));
    std::printf("  K0 install (no excavation): max|wall ux| = %.3e\n", max_a);
    check(max_a < 1e-8, "seeded wall does not move at K0 (no spurious installation deflection)");

    // (B) Excavation: the flexible barrier (interface allows slip/separation) deflects toward the
    // cut over the loaded zone -- a genuine, excavation-driven response (NOT the install artifact).
    const Out b = run(true);
    check(b.converged, "excavation phase converged");
    double max_b = 0.0, net = 0.0;
    int imax = 0;
    for (size_t i = 0; i < b.y.size(); ++i) {
        net += b.ux[i];
        if (std::fabs(b.ux[i]) > max_b) { max_b = std::fabs(b.ux[i]); imax = static_cast<int>(i); }
        if (std::fmod(b.y[i], 2.0) < 1e-9)
            std::printf("     y=%5.2f  ux=%+.4e\n", b.y[i], b.ux[i]);
    }
    std::printf("  excavation: max|wall ux|=%.4e at y=%.2f, net sum=%.4e\n", max_b, b.y[imax], net);
    check(b.ux[imax] < 0.0, "max wall deflection is toward the cut (ux < 0)");
    check(net < 0.0, "net wall deflection is toward the cut");
    check(max_b > 1e-3, "excavation produces a meaningful wall deflection");
    check(max_b > 100.0 * max_a, "deflection is excavation-driven (>> install artifact)");
}

} // namespace

int main() {
    test_wall_k0_excavation();
    if (g_failures == 0) {
        std::printf("OK: barrier wall + K0 interface seeding + staged excavation verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
