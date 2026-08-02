// Hardening Soil in AXISYMMETRY -- PLAXIS supports HS in both plane strain and axisymmetric
// analysis, so parity ("istisnasiz") requires HS at the axisymmetric material point and BVP.
// HS was previously wired only into plane strain (integrate_point); integrate_point_axisym now
// reuses the SAME kinematics-agnostic principal return (hs_return_core): the (r,z) block is the
// in-plane Mohr circle and the hoop sigma_theta is the third principal, with the consistent 4x4
// tangent D_T = algo_jacobian * D_e_axisym (mirroring Mohr-Coulomb).
//
// (A) MATERIAL POINT -- axisymmetric oedometer (eps_r = eps_theta = 0, axial eps_z): the axial
//     stress grows and the radial/hoop follow with K0 = sigma_r/sigma_z. This is EXACTLY the
//     stress path of hs_oedometer_probe (principal axis 1 loaded, 2=3 confined), so the developed
//     K0 must match the calibrated K0^NC -- proving the principal return is identical across
//     kinematics. Hardening is monotone and the deviator stays bounded by qf(sigma3).
// (B) BVP -- a circular footing on HS sand (axisymmetric), the PLAXIS Tutorial L1 geometry, under
//     a prescribed rigid settlement: the analytic 4x4 consistent tangent must drive the global
//     Newton to convergence with a physical, bounded footing reaction.
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
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
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

HardeningSoilParams sand() {
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.nu_ur = 0.2; p.friction = 32 * kDeg;
    p.dilatancy = 2 * kDeg; p.cohesion = 5.0; p.Rf = 0.9;
    hs_calibrate_cap(p, 1.0 - std::sin(p.friction));
    return p;
}

MaterialModel hs_model(const HardeningSoilParams& hs) {
    MaterialModel m;
    m.type = MaterialType::HardeningSoil;
    m.youngs_modulus = hs.Eur_ref; m.poisson_ratio = hs.nu_ur;
    m.friction_angle = hs.friction; m.dilatancy_angle = hs.dilatancy; m.cohesion = hs.cohesion;
    m.hs = hs;
    return m;
}

// (A) Axisymmetric oedometer at a single material point -> developed K0 must match the
// reference probe (and hence the calibrated K0^NC). Strain order [r, z, rz, theta].
void test_axisym_oedometer_K0() {
    const HardeningSoilParams hs = sand();
    const MaterialModel m = hs_model(hs);

    // Reference K0 from the principal-space oedometer probe (axis 1 loaded, 2=3 confined).
    double Eoed_ref = 0.0, K0_ref = 0.0;
    hs_oedometer_probe(hs, Eoed_ref, K0_ref);

    // Material point: start from a small isotropic compressive stress on the cap, seed history.
    const double p0 = 0.02 * hs.p_ref;
    GaussState s;
    s.stress << -p0, -p0, 0.0;  // [sigma_r, sigma_z, sigma_rz] (tension-positive)
    s.stress_zz = -p0;          // sigma_theta
    const Eigen::Vector3d sig_cp(p0, p0, p0);
    s.gamma_p = hs_initial_gamma_p(hs, sig_cp);
    s.pp = hs_initial_pp(hs, sig_cp);

    const double de = 2.0e-4;  // axial (z) compression increment; r and theta confined
    double K0 = 0.0; bool monotone = true, bounded = true, symmetric = true;
    double prev_gp = s.gamma_p;
    for (int i = 0; i < 2000; ++i) {
        GaussState tr; Eigen::Matrix4d tan;
        integrate_point_axisym(m, s, Eigen::Vector4d(0.0, -de, 0.0, 0.0), tr, tan);
        if (tr.gamma_p < prev_gp - 1e-10) monotone = false;
        prev_gp = tr.gamma_p;
        s = tr;
        // q = sigma_z - sigma_r (compression-positive deviator); qf(sigma3=sigma_r).
        const double sz = -s.stress(1), sr = -s.stress(0);
        const double q = sz - sr, qf = hs.q_failure(sr);
        if (q > qf * 1.001) bounded = false;
        // Radial (in-plane) and hoop (out-of-plane) confinement must stay equal: the oedometer
        // is the sigma2=sigma3 triaxial edge. (The bare-gradient drift used to break this.)
        if (std::fabs(s.stress(0) - s.stress_zz) > 1e-6 * (1.0 + std::fabs(s.stress(0))))
            symmetric = false;
        if (-s.stress(1) >= hs.p_ref) {
            K0 = s.stress(0) / s.stress(1);
            std::printf("    [at capture: sigma_r=%.3f sigma_z=%.3f sigma_theta=%.3f]\n",
                        s.stress(0), s.stress(1), s.stress_zz);
            break;
        }
    }
    std::printf("  axisym oedometer: K0(FE) = %.4f  K0(probe) = %.4f  (calib K0nc = %.4f)\n",
                K0, K0_ref, 1.0 - std::sin(hs.friction));
    check(monotone, "axisym HS oedometer: shear hardening is monotone");
    check(bounded, "axisym HS oedometer: deviator stays bounded by qf(sigma3)");
    check(symmetric, "axisym HS oedometer keeps sigma_r == sigma_theta (sigma2=sigma3 edge)");
    check(std::fabs(K0 - K0_ref) < 0.03,
          "axisym HS oedometer develops the same K0 as the principal-space probe");
}

// Axisymmetric consistent nodal internal force F = integral B^T sigma * (w detJ r), full DOF.
Eigen::VectorXd axisym_force(const Mesh& mesh, const std::vector<GaussState>& gs) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const auto gp = Tri6Element::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        Tri6Element::NodeCoords X;
        for (int k = 0; k < 6; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        Eigen::Matrix<double, 12, 1> fe = Eigen::Matrix<double, 12, 1>::Zero();
        for (int g = 0; g < 3; ++g) {
            const auto sg = axisym::strain_displacement<Tri6Element>(X, gp[g].xi, gp[g].eta);
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

// (B) Axisymmetric circular-footing BVP on HS sand: must converge with a physical reaction.
void test_axisym_footing_bvp() {
    constexpr double R_model = 5.0, H = 4.0, R_foot = 1.0, gamma = 18.0;
    constexpr double settle = 0.02;  // prescribed rigid footing settlement [m]
    const HardeningSoilParams hs = sand();
    const double K0 = 1.0 - std::sin(hs.friction);

    const RectangularDomain domain{0.0, 0.0, R_model, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 12, 10);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);   // axis r=0
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    int n_foot = 0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) {
            dofs.fix_node_component(n, 1);
            presc[dofs.global_dof(n, 1)] = -settle;
            ++n_foot;
        }
    dofs.finalize();

    // K0 initial effective stress + admissible HS history seed.
    const auto gpts = Tri6Element::gauss_points();
    std::vector<GaussState> init(static_cast<size_t>(mesh.element_count) * 3);
    for (int e = 0; e < mesh.element_count; ++e)
        for (int g = 0; g < 3; ++g) {
            const auto N = Tri6Element::shape_functions(gpts[g].xi, gpts[g].eta);
            double y = 0.0;
            for (int k = 0; k < 6; ++k) y += N(k) * mesh.y[mesh.node_of(e, k)];
            const double sv = -gamma * (H - y);
            GaussState& s = init[e * 3 + g];
            s.stress << K0 * sv, sv, 0.0; s.stress_zz = K0 * sv;
            const Eigen::Vector3d sig_cp(-K0 * sv, -sv, -K0 * sv);
            s.gamma_p = hs_initial_gamma_p(hs, sig_cp);
            s.pp = hs_initial_pp(hs, sig_cp);
        }

    // Unified K0: hold the geostatic internal force constant, ramp only the prescribed settlement.
    const Eigen::VectorXd F0_full = axisym_force(mesh, init);
    Eigen::VectorXd cf = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int n = 0; n < mesh.node_count; ++n)
        for (int c = 0; c < 2; ++c) {
            const int eq = dofs.equation(dofs.global_dof(n, c));
            if (eq >= 0) cf(eq) = F0_full(2 * n + c);
        }

    const std::vector<MaterialModel> mm = {hs_model(hs)};
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), presc.size());
    NewtonOptions opt{40, 80, 1e-2}; opt.kinematics = Kinematics::Axisymmetric;
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());

    const auto r = solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, opt, init, {}, {}, pv, cf);

    const Eigen::VectorXd Ffin = axisym_force(mesh, r.gauss_states);
    double Ry = 0.0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) Ry += (Ffin(2 * n + 1) - F0_full(2 * n + 1));
    const double total = std::fabs(Ry) * 2.0 * kPi;
    std::printf("  axisym HS footing: converged=%d load_factor=%.3f iters=%d  reaction=%.1f kN\n",
                r.converged, r.load_factor, r.total_iterations, total);
    check(r.converged && r.load_factor > 0.999, "axisym HS footing BVP converges");
    check(total > 0.0 && total < 1.0e5, "axisym HS footing reaction is positive and bounded");
}

}  // namespace

int main() {
    test_axisym_oedometer_K0();
    test_axisym_footing_bvp();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil in axisymmetry (material point + BVP) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
