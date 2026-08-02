// Lateral earth pressure on a rigid wall -- Rankine active/passive (P2.4 / supported
// excavation foundation). A smooth vertical wall (the left boundary x=0) retains a
// cohesionless Mohr-Coulomb soil under its K0 geostatic state. The wall is translated
// horizontally (prescribed displacement, ramped while gravity is held constant); the
// recovered horizontal thrust is compared to the classical Rankine resultants:
//   at rest  P0 = 1/2 K0 gamma H^2 ,  K0 = 1 - sin(phi)         (Jaky)
//   active   Pa = 1/2 Ka gamma H^2 ,  Ka = tan^2(45 - phi/2)
//   passive  Pp = 1/2 Kp gamma H^2 ,  Kp = tan^2(45 + phi/2) = 1/Ka
// Smooth wall (uy free) -> no wall friction -> Rankine is exact. Exercises the new
// prescribed-displacement + constant-force solver path, K0 initial stress, and MC
// plasticity at the active/passive limit. Reference: Rankine (1857); Das, Craig.
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::GaussState;
using katai::core::K0Options;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace tri6 = katai::core::tri6;

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

// Consistent nodal internal force F = sum_e integral B^T sigma over elements (full DOF).
Eigen::VectorXd nodal_internal_force(const Mesh& mesh,
                                     const std::vector<GaussState>& gs) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const auto gp = tri6::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        tri6::NodeCoords X;
        for (int k = 0; k < 6; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        Eigen::Matrix<double, 12, 1> fe = Eigen::Matrix<double, 12, 1>::Zero();
        for (int g = 0; g < 3; ++g) {
            const auto sg = tri6::strain_displacement(X, gp[g].xi, gp[g].eta);
            const Eigen::Vector3d sigma = gs[e * 3 + g].stress;
            fe.noalias() += (gp[g].weight * sg.det_jacobian) * sg.B.transpose() * sigma;
        }
        for (int k = 0; k < 6; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
    return F;
}

// Translate the wall (left boundary) by delta (m); return the horizontal thrust on it.
double wall_thrust(double delta, const char* label) {
    constexpr double W = 12.0, H = 10.0, gamma = 18.0;
    const double phi = 30.0 * std::acos(-1.0) / 180.0, c = 0.0, psi = phi;  // associated (robust)
    constexpr double E = 2.0e4, nu = 0.3;
    const double K0 = 1.0 - std::sin(phi);

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 16, 14);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) dofs.fix_node_component(n, 1);          // vertical support
    for (int n : mesh.right_nodes)  dofs.fix_node_component(n, 0);          // far-field roller
    for (int n : mesh.left_nodes) {                                        // the wall
        dofs.fix_node_component(n, 0);                                      // ux prescribed
        presc[dofs.global_dof(n, 0)] = delta;                              // (uy free -> smooth)
    }
    dofs.finalize();

    const std::vector<GaussState> init =
        katai::core::compute_k0_initial_stress(mesh, K0Options{H, gamma, K0});
    Eigen::VectorXd grav = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {gamma}, grav);   // constant (geostatic) force

    Eigen::Map<const Eigen::VectorXd> presc_v(presc.data(), presc.size());
    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, E, nu, c, phi, psi}};
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f0, solve_unsym,
                                                {100, 60, 1e-6}, init, {}, {},
                                                presc_v, grav);

    const Eigen::VectorXd F = nodal_internal_force(mesh, r.gauss_states);
    double thrust = 0.0;
    for (int n : mesh.left_nodes) thrust += F(2 * n + 0);  // horizontal reaction on wall
    std::printf("  %-8s delta=%+.3f  conv=%d lf=%.2f  thrust=%.1f\n",
                label, delta, (int)r.converged, r.load_factor, std::fabs(thrust));
    return std::fabs(thrust);
}

void test_rankine_earth_pressure() {
    constexpr double H = 10.0, gamma = 18.0;
    const double pi = std::acos(-1.0);
    const double phi = 30.0 * pi / 180.0;
    const double Ka = std::tan(pi / 4 - phi / 2) * std::tan(pi / 4 - phi / 2);
    const double Kp = 1.0 / Ka;
    const double K0 = 1.0 - std::sin(phi);
    const double P0 = 0.5 * K0 * gamma * H * H;
    const double Pa = 0.5 * Ka * gamma * H * H;
    const double Pp = 0.5 * Kp * gamma * H * H;
    std::printf("  Rankine: Ka=%.3f Kp=%.3f K0=%.3f -> Pa=%.1f P0=%.1f Pp=%.1f\n",
                Ka, Kp, K0, Pa, P0, Pp);

    const double t0 = wall_thrust(0.0, "at-rest");
    const double ta = wall_thrust(-0.05, "active");
    const double tp = wall_thrust(+0.40, "passive");

    check(std::fabs(t0 - P0) < 0.02 * P0, "at-rest thrust = 1/2 K0 gamma H^2 (Jaky)");
    check(std::fabs(ta - Pa) < 0.02 * Pa, "active thrust = 1/2 Ka gamma H^2 (Rankine)");
    check(std::fabs(tp - Pp) < 0.03 * Pp, "passive thrust = 1/2 Kp gamma H^2 (Rankine)");
    check(ta < t0 && t0 < tp, "active < at-rest < passive (monotone)");
}

} // namespace

int main() {
    test_rankine_earth_pressure();
    if (g_failures == 0) {
        std::printf("OK: Rankine active/passive earth pressure verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
