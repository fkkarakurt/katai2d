// EC7 (EN 1997-1) material-factoring verification -- v0.3 workstream B1.2.
//
// Verifies that the design-code MATERIAL factoring (design_code.hpp) produces the correct,
// EXPECTED design result on a real boundary-value problem, cross-checked against:
//   (a) the natural analytical result -- an EXACT composition identity;
//   (b) a published slope-stability benchmark (Griffiths & Lane / Rocscience);
//   (c) the algorithm PLAXIS 2D uses for its "Design Approaches" facility.
//
// The identity. Strength reduction (SRM) defines the factor of safety FoS as the factor S by which
// BOTH c and tan(phi) are divided to reach incipient collapse:  strength / FoS  at collapse.
// EC7 DA3 pre-divides c and tan(phi) by the M2 partial factor gamma_M = 1.25 (gamma_c' = gamma_phi'),
// then the residual "over-design factor" ODF is again an SRM factor on the ALREADY-factored strength.
// Because the solver sees the identical effective (c, tan phi) at a given TOTAL factor regardless of
// how it is reached, collapse occurs at the same total factor:
//        gamma_M * ODF  =  FoS_characteristic       =>   ODF = FoS_char / 1.25   (EXACT, mesh-free).
// PLAXIS documents the same mechanism: "cohesion, friction angle and dilatancy are reduced using the
// partial factor" (Bentley/PLAXIS Design Approaches). So KATAI's result agrees with PLAXIS by
// construction of the identical algorithm. EC7 criterion: the design is safe iff ODF >= 1.0.
//
// Benchmark (from test_slope): homogeneous 1:2 foundation slope, gamma = 20.2, c = 3 kPa,
// phi = 19.6 deg, psi = 0. LEM FoS ~ 0.987-0.988 (Bishop/Spencer), Phase2 T6 0.997.
#include <katai/analysis/design_code.hpp>
#include <katai/analysis/strength_reduction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DesignApproach;
using katai::core::design_factors;
using katai::core::DofMap;
using katai::core::factor_material_strength;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::StrengthReductionOptions;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
using katai::mesh::Triangulation;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
constexpr double kDeg = 3.14159265358979323846 / 180.0;
}  // namespace

int main() {
    // --- Slope model (identical geometry to test_slope) --------------------------------------
    const std::vector<double> px = {20.0, 70.0, 70.0, 50.0, 30.0, 20.0};
    const std::vector<double> py = {20.0, 20.0, 35.0, 35.0, 25.0, 25.0};
    const std::vector<std::array<int, 2>> segs = {
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 4}}, {{4, 5}}, {{5, 0}}};
    const Triangulation tg = katai::mesh::quality_mesh(px, py, segs, 20.0, 5.0);
    const Mesh mesh = katai::mesh::tri6_from_triangulation(tg, 0);
    std::printf("EC7 material-factoring verification (DA3 slope stability)\n");
    std::printf("  mesh: %d elements, %d nodes\n", mesh.element_count, mesh.node_count);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.boundary_nodes) {
        const double x = mesh.x[n], y = mesh.y[n];
        if (std::fabs(y - 20.0) < 1e-7) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        if (std::fabs(x - 20.0) < 1e-7 || std::fabs(x - 70.0) < 1e-7) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();

    Eigen::VectorXd gravity = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {20.2}, gravity);

    const MaterialModel base{MaterialType::MohrCoulomb, 1.0e5, 0.3, 3.0, 19.6 * kDeg, 0.0};

    auto solve = [](const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(k);
        return s->solve(r);
    };

    StrengthReductionOptions opt;
    opt.srf_min = 0.3;   // ODF ~ 0.79 must stay inside the bracket
    opt.srf_max = 2.5;
    opt.bisection_iterations = 13;  // resolution ~2.2/2^13 ~= 2.7e-4 (tight enough for the identity)
    opt.newton.load_steps = 4;
    opt.newton.max_iterations = 150;
    opt.newton.tolerance = 1e-3;

    // --- (1) Characteristic factor of safety (unfactored strength) -----------------------------
    const double fos_char = katai::core::factor_of_safety(mesh, dofs, gravity, base, solve, opt);

    // --- (2) EC7 DA3 over-design factor: SRM on the M2-factored strength -----------------------
    MaterialModel da3 = base;
    const auto f = design_factors(DesignApproach::EC7_DA3);
    factor_material_strength(da3, f);  // c/1.25, atan(tan phi/1.25); gamma_G = 1.0 -> gravity unchanged
    const double odf = katai::core::factor_of_safety(mesh, dofs, gravity, da3, solve, opt);

    const double gamma_M = 1.25;  // M2 (drained: gamma_c' = gamma_phi')
    const double predicted_odf = fos_char / gamma_M;
    const double id_err = std::fabs(odf - predicted_odf) / predicted_odf;

    std::printf("\n  FoS (characteristic)      = %.4f   (LEM benchmark ~0.99)\n", fos_char);
    std::printf("  ODF (EC7 DA3, M2 factored) = %.4f\n", odf);
    std::printf("  identity  FoS/gamma_M      = %.4f   (gamma_M = 1.25)\n", predicted_odf);
    std::printf("  identity error             = %.3f%%   (expected ~0, mesh-free)\n", 100.0 * id_err);
    std::printf("  EC7 verdict: ODF %s 1.0  ->  the slope %s the DA3 ULS check\n",
                odf >= 1.0 ? ">=" : "<", odf >= 1.0 ? "SATISFIES" : "FAILS");

    // --- Checks --------------------------------------------------------------------------------
    // The composition identity is the core correctness proof (mesh-independent, analytical).
    check(id_err < 0.015, "ODF * gamma_M == FoS_characteristic (EC7 DA3 material-factoring identity)");
    // Bonus sanity: the characteristic FoS reproduces the published benchmark.
    check(std::fabs(fos_char - 0.99) / 0.99 < 0.05, "characteristic FoS within 5% of LEM benchmark");
    // The design conclusion must be the correct one for these (low) characteristic strengths.
    check(odf < 1.0, "EC7 DA3 correctly reports this slope as unsafe (ODF < 1)");

    if (g_failures) { std::fprintf(stderr, "\n%d check(s) FAILED\n", g_failures); return 1; }
    std::printf("\nOK: EC7 DA3 material-factoring verified (analytical identity + benchmark + PLAXIS-consistent).\n");
    return 0;
}
