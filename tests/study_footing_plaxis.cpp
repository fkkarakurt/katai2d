// CAPSTONE: PLAXIS 2D Tutorial Lesson 1 "Settlement of a circular footing on sand", Case A
// (rigid footing), reproduced from scratch and compared to the documented PLAXIS result.
//
// Axisymmetric, Mohr-Coulomb (drained) sand: gamma=17/20, E'=13 MPa, nu=0.3, c'=1 kPa, phi'=30,
// psi=0. 4 m sand over rigid base, model radius 5 m, footing radius 1 m. Water table at y=2.
// K0 procedure (K0=1-sin30=0.5) for the initial stress; Phase 1 = a prescribed uniform vertical
// settlement of 0.05 m over the footing (rigid footing). PLAXIS result: total footing reaction
// Force-Y x 2*pi ~ 588 kN.
//
// Exercises: axisymmetry + Mohr-Coulomb plasticity + K0 initial stress (with water table) +
// prescribed displacement. Reports ACCURACY (vs 588 kN) and PERFORMANCE (DOFs, iterations, time).
// Build/run: cmake --build build/msvc-rwdi --target study_footing_plaxis && bin/study_footing_plaxis
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace katai::core;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

static Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// Axisymmetric consistent nodal internal force F = integral B^T sigma * (w * detJ * r), full DOF.
static Eigen::VectorXd axisym_nodal_force(const Mesh& mesh, const std::vector<GaussState>& gs) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const auto gp = Tri6Element::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        Tri6Element::NodeCoords X;
        for (int k = 0; k < 6; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        Eigen::Matrix<double, 12, 1> fe = Eigen::Matrix<double, 12, 1>::Zero();
        for (int g = 0; g < 3; ++g) {
            const auto sg = katai::core::axisym::strain_displacement<Tri6Element>(X, gp[g].xi, gp[g].eta);
            Eigen::Vector4d sig;
            const auto& s = gs[e * 3 + g];
            sig << s.stress(0), s.stress(1), s.stress(2), s.stress_zz;
            fe.noalias() += (gp[g].weight * sg.det_jacobian * sg.radius) * sg.B.transpose() * sig;
        }
        for (int k = 0; k < 6; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
    return F;
}

int main() {
    constexpr double R_model = 5.0, H = 4.0, R_foot = 1.0, y_wat = 2.0;
    constexpr double gamma_d = 17.0, gamma_sat = 20.0, gamma_w = 10.0;
    constexpr double E = 13.0e3, nu = 0.3, c = 1.0, gw = 0.0; (void)gw;
    const double pi = std::acos(-1.0);
    const double phi = 30.0 * pi / 180.0, psi = 0.0;     // PLAXIS: psi=0 (non-associated)
    const double K0 = 1.0 - std::sin(phi);                // 0.5
    constexpr double settle = 0.05;                       // prescribed footing settlement [m]
    const double plaxis_total = 588.0;                    // documented total reaction [kN]

    // Axisymmetric mesh: x = radius. r=1 (footing), y=2 (water), y=4 (top) land on node lines.
    const RectangularDomain domain{0.0, 0.0, R_model, H, 0};
    const int nx = 20, ny = 16;                           // fine 0.125 m
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, ny);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);   // axis r=0: roller
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);   // far field roller
    // Footing: rigid -> prescribe uniform vertical settlement over r <= R_foot at the top (y=H).
    int n_foot = 0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) {
            dofs.fix_node_component(n, 1);                          // vertical prescribed
            presc[dofs.global_dof(n, 1)] = -settle;                // downward
            ++n_foot;
        }
    dofs.finalize();

    // K0 initial effective stress (two-layer: dry above the water table, buoyant below).
    const auto gpts = Tri6Element::gauss_points();
    std::vector<GaussState> init(static_cast<size_t>(mesh.element_count) * 3);
    for (int e = 0; e < mesh.element_count; ++e)
        for (int g = 0; g < 3; ++g) {
            const auto N = Tri6Element::shape_functions(gpts[g].xi, gpts[g].eta);
            double y = 0.0;
            for (int k = 0; k < 6; ++k) y += N(k) * mesh.y[mesh.node_of(e, k)];
            const double sv = (y >= y_wat) ? -gamma_d * (H - y)
                                           : -(gamma_d * (H - y_wat) + (gamma_sat - gamma_w) * (y_wat - y));
            GaussState& s = init[e * 3 + g];
            s.stress << K0 * sv, sv, 0.0; s.stress_zz = K0 * sv;
        }

    // Staged-release baseline: constant_force = f_int0 of the K0 seed (holds the geostatic state,
    // residual 0 at u=0); the ramped prescribed footing settlement develops the reaction on top.
    const Eigen::VectorXd F0_full = axisym_nodal_force(mesh, init);
    Eigen::VectorXd cf = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int n = 0; n < mesh.node_count; ++n)
        for (int c2 = 0; c2 < 2; ++c2) {
            const int eq = dofs.equation(dofs.global_dof(n, c2));
            if (eq >= 0) cf(eq) = F0_full(2 * n + c2);
        }

    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, E, nu, c, phi, psi}};
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), presc.size());
    NewtonOptions opt{60, 60, 1e-6}; opt.kinematics = Kinematics::Axisymmetric;
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());

    const auto t0 = std::chrono::steady_clock::now();
    const auto r = solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, opt, init, {}, {}, pv, cf);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // Total footing reaction = (final - initial) vertical nodal force at footing nodes, x 2*pi.
    const Eigen::VectorXd Ffin = axisym_nodal_force(mesh, r.gauss_states);
    double Ry = 0.0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) Ry += (Ffin(2 * n + 1) - F0_full(2 * n + 1));
    const double total = std::fabs(Ry) * 2.0 * pi;

    std::printf("\n=== PLAXIS Tutorial Lesson 1: circular footing on sand (Case A) ===\n");
    std::printf("  mesh: %d elements, %d nodes, %d DOFs (axisym tri6); footing nodes=%d\n",
                mesh.element_count, mesh.node_count, dofs.equation_count(), n_foot);
    std::printf("  ACCURACY: KATAI total footing force = %.1f kN   PLAXIS = %.1f kN   (%.1f%%)\n",
                total, plaxis_total, 100.0 * (total - plaxis_total) / plaxis_total);
    std::printf("  PERFORMANCE: converged=%d  load_factor=%.3f  total iters=%d  solve time=%.2f s\n",
                (int)r.converged, r.load_factor, r.total_iterations, secs);
    return 0;
}
