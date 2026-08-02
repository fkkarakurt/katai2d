// OFFICIAL PLAXIS VALIDATION MANUAL cases, reproduced in KATAI (Benchmark Wave-1 / ROADMAP
// Track 7). Source: PLAXIS Validation Manual (Version 8), the company's own published benchmark
// set -- each case states BOTH the analytical solution and PLAXIS's published number, so the table
// docs/validation/plaxis-comparison-table.md can show analytic | PLAXIS | KATAI side by side.
//
//   (2.1) Smooth rigid strip footing on elastic soil (Giroud 1972).
//         G=500 kPa, nu=0.333, H=4 m, B=2 m, prescribed uy=10 mm, weightless.
//         Analytic F = 2(1+nu) G B s / rho, rho=0.88 -> 15.15 kN/m. PLAXIS: 15.24 (+0.6%).
//   (2.2) Strip load on elastic GIBSON soil (Gibson 1967): G=100 z (E=299 z), nu=0.495,
//         q=10 kPa on B=2 m, H=4 m. Exact (half-space): s = q/(2 dG/dz) = 0.050 m -- uniform
//         under the load. PLAXIS (finite 4 m layer): 0.047 m (-6%). The SAME finite-layer bias
//         must appear in KATAI: this pins our E(y) machinery against their Advanced E-increment.
//   (3.1) Bearing capacity of a smooth rigid CIRCULAR footing (Cox 1962): axisymmetric MC,
//         E=2400, nu=0.2, c=1.6, phi=30, psi=0, gamma=16 (gamma R/c = 10), K0=0.5, R=1, H=4.
//         Analytic p_max = 141 c = 225.6 kPa. PLAXIS: 220 (-2.5%).
//   (3.2) Strip footing on clay with strength increasing with depth (Davis & Booker 1973):
//         c = 1 + 2 z kPa, E = 299 + 498 z, weightless, B=2 m, smooth.
//         Analytic p_max = rho[(2+pi) c0 + B c_inc/4] = 1.27 * 6.1416 = 7.80 kPa. PLAXIS: 7.86.
//         Pins the c(y) + E(y) profile machinery against a THEORETICAL collapse load.
//
// Method notes: displacement-controlled (prescribed footing settlement, reactions from the
// committed-stress nodal internal force -- the study_footing_plaxis pattern), tri15 elements
// (PLAXIS's own 15-node default). Bands stated per case and honest.
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace katai::core;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// Full-DOF nodal internal force from committed Gauss stresses (reactions live on the fixed DOFs).
// Plane strain: F = integral B^T sigma dA; axisym: r-weighted per radian, with the hoop component.
template <class E>
Eigen::VectorXd nodal_force_ps(const Mesh& mesh, const std::vector<GaussState>& gs) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        typename E::NodeCoords X;
        for (int k = 0; k < E::kNodeCount; ++k) {
            X(k, 0) = mesh.x[mesh.node_of(e, k)];
            X(k, 1) = mesh.y[mesh.node_of(e, k)];
        }
        Eigen::Matrix<double, E::kDofCount, 1> fe = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sg = E::strain_displacement(X, gp[g].xi, gp[g].eta);
            fe.noalias() += (gp[g].weight * sg.det_jacobian) *
                            sg.B.transpose() * gs[(size_t)e * E::kGaussCount + g].stress;
        }
        for (int k = 0; k < E::kNodeCount; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
    return F;
}
template <class E>
Eigen::VectorXd nodal_force_axi(const Mesh& mesh, const std::vector<GaussState>& gs) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        typename E::NodeCoords X;
        for (int k = 0; k < E::kNodeCount; ++k) {
            X(k, 0) = mesh.x[mesh.node_of(e, k)];
            X(k, 1) = mesh.y[mesh.node_of(e, k)];
        }
        Eigen::Matrix<double, E::kDofCount, 1> fe = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sg = axisym::strain_displacement<E>(X, gp[g].xi, gp[g].eta);
            const auto& s = gs[(size_t)e * E::kGaussCount + g];
            Eigen::Vector4d sig;
            sig << s.stress(0), s.stress(1), s.stress(2), s.stress_zz;
            fe.noalias() += (gp[g].weight * sg.det_jacobian * sg.radius) * sg.B.transpose() * sig;
        }
        for (int k = 0; k < E::kNodeCount; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
    return F;
}

// ============================================================================================
// verify: KV-FND-001
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Giroud (1972)
//   locator:  Section 2.1, smooth rigid strip footing on elastic soil
//   quantity: footing force F at prescribed settlement u_y = 10 mm [kN/m]
//   expected: 15.15 (analytic, F = 2(1+nu) G B s / rho with rho = 0.88); PLAXIS publishes 15.24
//   band:     2% vs analytic and 3% vs PLAXIS, as asserted below -- measured +1.4% / +0.8% (KATAI 15.35)
// (2.1) Smooth rigid strip footing on elastic soil -- Giroud 15.15, PLAXIS 15.24.
// ============================================================================================
void case_2_1() {
    std::printf("-- (2.1) smooth rigid strip footing on elastic soil (Giroud) --\n");
    constexpr double G = 500.0, nu = 1.0 / 3.0, s0 = 0.010, Bhalf = 1.0;
    const double E = 2.0 * G * (1.0 + nu);
    const RectangularDomain dom{0.0, 0.0, 7.0, 4.0, 0};
    Mesh mesh = katai::mesh::generate_structured_tri15(dom, 14, 8);   // 0.5 m grid, x=1 on a line

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);    // symmetry axis
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);   // far-field roller
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= Bhalf + 1e-9) {                            // smooth: uy prescribed, ux free
            dofs.fix_node_component(n, 1);
            presc[dofs.global_dof(n, 1)] = -s0;
        }
    dofs.finalize();

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu, 0.0, 0.0, 0.0}};
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), presc.size());
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto r = solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, {1, 30, 1e-10}, {}, {}, {}, pv);
    check(r.converged, "elastic indentation solved");

    const Eigen::VectorXd F = nodal_force_ps<Tri15Element>(mesh, r.gauss_states);
    double Ry = 0.0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= Bhalf + 1e-9) Ry += F(2 * n + 1);
    const double F_katai = 2.0 * std::fabs(Ry);                     // half-model -> full footing
    const double F_exact = 15.15, F_plaxis = 15.24;
    std::printf("   F: analytic %.2f | PLAXIS %.2f (+0.6%%) | KATAI %.2f (%+.1f%% vs analytic, "
                "%+.1f%% vs PLAXIS)\n", F_exact, F_plaxis, F_katai,
                100 * (F_katai - F_exact) / F_exact, 100 * (F_katai - F_plaxis) / F_plaxis);
    check(std::fabs(F_katai - F_exact) < 0.02 * F_exact,
          "KATAI within 2% of the Giroud analytic footing force");
    check(std::fabs(F_katai - F_plaxis) < 0.03 * F_plaxis,
          "KATAI within 3% of the published PLAXIS number");
}

// ============================================================================================
// verify: KV-FND-002
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Gibson (1967)
//   locator:  Section 2.2, strip load on incompressible Gibson soil (E = 299 z, nu = 0.495)
//   quantity: settlement under a q = 10 kPa strip load on a 4 m layer [m]
//   expected: 0.047 (PLAXIS, same finite layer); the half-space closed form gives 0.050
//   band:     5% vs the PLAXIS finite-layer value, as asserted below -- measured -4.0% (KATAI 0.0451); the finite-layer bias vs the half-space is shared by construction and pinned separately below
// (2.2) Strip load on Gibson soil -- exact 0.050 m (half-space), PLAXIS 0.047 m (finite layer).
// ============================================================================================
void case_2_2() {
    std::printf("-- (2.2) strip load on elastic Gibson soil (E = 299 z, nu = 0.495) --\n");
    constexpr double q = 10.0, Bhalf = 1.0, H = 4.0;
    const RectangularDomain dom{0.0, 0.0, 7.0, H, 0};
    // 0.125 m grid: with E -> 0 at the surface the compliance is concentrated in a thin surface
    // band, and nu = 0.495 adds mild volumetric locking -- a 0.5 m grid under-settles by ~7%,
    // 0.25 m by ~5.5%; the elastic solve is a single factorization, so the fine grid is cheap.
    Mesh mesh = katai::mesh::generate_structured_tri15(dom, 56, 32);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // E(y) = E_ref + E_inc (y_ref - y): tiny at the surface, 299 per metre of depth (their
    // "Advanced E-increment" input, our MaterialProfile -- the same modelling decision).
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, 0.01, 0.495, 0.0, 0.0, 0.0}};
    const std::vector<MaterialProfile> prof = {{299.0, 0.0, H}};   // {E_inc, c_inc, y_ref}

    // q over x <= 1 on the top edge (ordered chain of the loaded part).
    std::vector<int> chain;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= Bhalf + 1e-9) chain.push_back(n);
    std::sort(chain.begin(), chain.end(),
              [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    assemble_surface_traction(mesh, dofs, chain, 0.0, -q, f);

    const auto r = solve_nonlinear(mesh, dofs, mm, f, solve_unsym, {1, 30, 1e-10}, {}, {}, {}, {},
                                   {}, prof);
    check(r.converged, "Gibson strip load solved");
    int nc = -1;
    for (int n : mesh.top_nodes)
        if (std::fabs(mesh.x[n]) < 1e-9) nc = n;
    const double s_katai = -r.displacement(2 * nc + 1);
    const double s_exact = 0.050, s_plaxis = 0.047;
    std::printf("   settlement: exact (half-space) %.4f | PLAXIS (4 m layer) %.4f (-6%%) | KATAI "
                "%.4f (%+.1f%% vs exact, %+.1f%% vs PLAXIS)\n", s_exact, s_plaxis, s_katai,
                100 * (s_katai - s_exact) / s_exact, 100 * (s_katai - s_plaxis) / s_plaxis);
    // The finite 4 m layer sits BELOW the half-space value -- KATAI must show the SAME bias as
    // PLAXIS (same layer, same E(z) modelling), i.e. land near 0.047, not near 0.050.
    check(std::fabs(s_katai - s_plaxis) < 0.05 * s_plaxis,
          "KATAI within 5% of the published PLAXIS settlement (same finite-layer model)");
    check(s_katai < s_exact, "the finite layer settles less than the half-space (shared bias)");
}

// ============================================================================================
// verify: KV-FND-003
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); slip-line solution Cox (1962)
//   locator:  Section 3.1, bearing capacity of a smooth rigid circular footing (axisymmetric Mohr-Coulomb)
//   quantity: limit pressure p_max [kPa], associated flow (psi = phi -- the slip-line solution is the associated limit load)
//   expected: 225.6 (analytic, 141 c); PLAXIS publishes 220.0
//   band:     5% vs analytic, as asserted below -- measured +3.9% at the 0.25 m mesh (KATAI 234.4); the 0.5 m mesh over-predicts by +9% (two elements across the radius; recorded)
// (3.1) Bearing capacity of a smooth circular footing -- Cox 225.6 kPa, PLAXIS 220 kPa.
// ============================================================================================
void case_3_1() {
    std::printf("-- (3.1) bearing capacity of a circular footing (Cox, axisym MC) --\n");
    constexpr double E = 2400.0, nu = 0.20, c = 1.6, gamma = 16.0, K0 = 0.5;
    constexpr double R_foot = 1.0, H = 4.0, settle = 0.35;
    const double phi = 30.0 * kPi / 180.0;
    // ASSOCIATED flow (psi = phi): Cox (1962) is a slip-line solution, i.e. the ASSOCIATED limit
    // load -- comparing it against a non-associated run would mix a modelling difference into a
    // verification number (and psi = 0 converges far more slowly near the limit; the wall-benchmark
    // lesson). The soft soil (E = 2400) also needs a large indentation to reach the plateau.
    const double psi = phi;
    const RectangularDomain dom{0.0, 0.0, 5.0, H, 0};
    // 0.25 m grid: a 0.5 m grid puts only two elements across the footing radius and OVER-predicts
    // the collapse by ~9% (the classic coarse-mesh bearing-capacity bias PLAXIS's own manual warns
    // about for low-order elements; even tri15 needs a few elements across the punch).
    Mesh mesh = katai::mesh::generate_structured_tri15(dom, 20, 16);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    int n_foot = 0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) {                           // smooth: uy only
            dofs.fix_node_component(n, 1);
            presc[dofs.global_dof(n, 1)] = -settle;
            ++n_foot;
        }
    dofs.finalize();

    // K0 geostatic seed (dry) + its internal force as the constant baseline.
    const auto gp = Tri15Element::gauss_points();
    std::vector<GaussState> init((size_t)mesh.element_count * Tri15Element::kGaussCount);
    for (int e = 0; e < mesh.element_count; ++e)
        for (int g = 0; g < Tri15Element::kGaussCount; ++g) {
            typename Tri15Element::NodeCoords X;
            for (int k = 0; k < 15; ++k) {
                X(k, 0) = mesh.x[mesh.node_of(e, k)];
                X(k, 1) = mesh.y[mesh.node_of(e, k)];
            }
            const auto N = Tri15Element::shape_functions(gp[g].xi, gp[g].eta);
            double y = 0.0;
            for (int k = 0; k < 15; ++k) y += N(k) * X(k, 1);
            const double sv = -gamma * (H - y);
            GaussState& s = init[(size_t)e * Tri15Element::kGaussCount + g];
            s.stress << K0 * sv, sv, 0.0;
            s.stress_zz = K0 * sv;
        }
    const Eigen::VectorXd F0 = nodal_force_axi<Tri15Element>(mesh, init);
    Eigen::VectorXd cf = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int n = 0; n < mesh.node_count; ++n)
        for (int c2 = 0; c2 < 2; ++c2) {
            const int eq = dofs.equation(dofs.global_dof(n, c2));
            if (eq >= 0) cf(eq) = F0(2 * n + c2);
        }

    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, E, nu, c, phi, psi}};
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), presc.size());
    NewtonOptions opt{50, 60, 1e-6};
    opt.kinematics = Kinematics::Axisymmetric;
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto r = solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, opt, init, {}, {}, pv, cf);
    check(r.converged && r.load_factor > 0.999, "pushed to 0.30 m, displacement-controlled");

    const Eigen::VectorXd Ff = nodal_force_axi<Tri15Element>(mesh, r.gauss_states);
    double Ry = 0.0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= R_foot + 1e-9) Ry += (Ff(2 * n + 1) - F0(2 * n + 1));
    const double p_katai = 2.0 * std::fabs(Ry) / (R_foot * R_foot);   // p = 2 F_rad / R^2
    const double p_exact = 225.6, p_plaxis = 220.0;
    std::printf("   p_max: Cox %.1f | PLAXIS %.1f (-2.5%%) | KATAI %.1f (%+.1f%% vs Cox, "
                "%+.1f%% vs PLAXIS)   [%d iters]\n", p_exact, p_plaxis, p_katai,
                100 * (p_katai - p_exact) / p_exact, 100 * (p_katai - p_plaxis) / p_plaxis,
                r.total_iterations);
    check(std::fabs(p_katai - p_exact) < 0.05 * p_exact,
          "KATAI within 5% of the Cox exact collapse pressure");
}

// ============================================================================================
// verify: KV-FND-004
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Davis & Booker (1973)
//   locator:  Section 3.2, smooth strip footing on clay with strength increasing with depth
//   quantity: limit pressure p_max [kPa] with c(z) = c0 + c_inc z
//   expected: 7.80 (analytic, rho [(2 + pi) c0 + B c_inc / 4]); PLAXIS publishes 7.86
//   band:     5% vs analytic, as asserted below -- measured +2.8% (KATAI 8.02)
// (3.2) Strip footing on clay with strength increasing with depth -- Davis & Booker 7.80 (smooth),
//       PLAXIS 7.86. Pins the c(y) + E(y) profiles against a theoretical collapse load.
// ============================================================================================
void case_3_2() {
    std::printf("-- (3.2) strip footing on c(z) clay (Davis & Booker, smooth) --\n");
    constexpr double Bhalf = 1.0, H = 4.0, settle = 0.03;
    const RectangularDomain dom{0.0, 0.0, 5.0, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri15(dom, 10, 8);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= Bhalf + 1e-9) {
            dofs.fix_node_component(n, 1);
            presc[dofs.global_dof(n, 1)] = -settle;
        }
    dofs.finalize();

    // Tresca (phi=0), c = 1 + 2 z, E = 299 + 498 z (their Advanced increments = our profiles).
    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, 299.0, 0.3, 1.0, 0.0, 0.0}};
    const std::vector<MaterialProfile> prof = {{498.0, 2.0, H}};   // {E_inc, c_inc, y_ref = top}

    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), presc.size());
    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto r = solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, {30, 60, 1e-8}, {}, {}, {}, pv,
                                   {}, prof);
    check(r.converged && r.load_factor > 0.999, "pushed to 0.03 m, displacement-controlled");

    const Eigen::VectorXd F = nodal_force_ps<Tri15Element>(mesh, r.gauss_states);
    double Ry = 0.0;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] <= Bhalf + 1e-9) Ry += F(2 * n + 1);
    const double p_katai = std::fabs(Ry) / Bhalf;                   // average stress under footing
    const double p_exact = 7.80, p_plaxis = 7.86;
    std::printf("   p_max: Davis-Booker %.2f | PLAXIS %.2f (+0.8%%) | KATAI %.2f (%+.1f%% vs "
                "analytic, %+.1f%% vs PLAXIS)   [%d iters]\n", p_exact, p_plaxis, p_katai,
                100 * (p_katai - p_exact) / p_exact, 100 * (p_katai - p_plaxis) / p_plaxis,
                r.total_iterations);
    check(std::fabs(p_katai - p_exact) < 0.05 * p_exact,
          "KATAI within 5% of the Davis-Booker exact collapse pressure");
}

}  // namespace

int main() {
    std::printf("OFFICIAL PLAXIS VALIDATION MANUAL cases reproduced in KATAI\n"
                "(analytic | PLAXIS published | KATAI -- the comparison-table backbone)\n\n");
    case_2_1();
    std::printf("\n");
    case_2_2();
    std::printf("\n");
    case_3_1();
    std::printf("\n");
    case_3_2();
    if (g_failures == 0) {
        std::printf("\nOK: KATAI reproduces the official PLAXIS validation set -- elastic footing, "
                    "Gibson E(z), Cox circular bearing, Davis-Booker c(z) bearing\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
