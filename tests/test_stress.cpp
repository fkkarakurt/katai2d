// Stress-recovery verification (P0.7 + P0.8a).
//
// Laterally confined column (exx=0). In this case:
//     syy = E_oed * eyy            (vertical)
//     sxx = (v/(1-v)) * syy        (yatay; elastik at-rest K0 = v/(1-v))
//     sxy = 0
// Surface pressure p → syy = -p (uniform). Since the stress field is constant/linear,
// the Gauss→node extrapolation is EXACT; the FEM result must match to round-off.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

constexpr double kE = 10000.0, kNu = 0.3, kH = 10.0, kW = 2.0, kPressure = 100.0;

double oedometer_modulus() {
    return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));
}

Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k,
                          const Eigen::VectorXd& f) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(f);
}

void test_confined_column_stress() {
    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 8);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) {
        dofs.fix_node_component(n, 0);
        dofs.fix_node_component(n, 1);
    }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{kE, kNu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    const katai::math::CsrMatrix k = builder.build();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0,
                                           -kPressure, f);

    const Eigen::VectorXd full =
        katai::core::expand_to_full(dofs, solve_spd(k, f));
    const auto field = katai::core::recover_nodal_stresses(mesh, materials, full);

    // Analitik: syy = -p, sxx = (v/(1-v))(-p), sxy = 0.
    const double syy_exact = -kPressure;
    const double sxx_exact = (kNu / (1.0 - kNu)) * syy_exact;

    double max_syy_err = 0.0, max_sxx_err = 0.0, max_sxy = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const Eigen::Vector3d& s = field.stress[n];
        max_sxx_err = std::fmax(max_sxx_err, std::fabs(s(0) - sxx_exact));
        max_syy_err = std::fmax(max_syy_err, std::fabs(s(1) - syy_exact));
        max_sxy = std::fmax(max_sxy, std::fabs(s(2)));
    }
    const double scale = std::fabs(syy_exact);
    check(max_syy_err < 1e-9 * scale, "syy = -p (uniform)");
    check(max_sxx_err < 1e-9 * scale, "sxx = K0 * syy");
    check(max_sxy < 1e-9 * scale, "sxy = 0");
    std::printf("  kolon: syy_err=%.2e  sxx_err=%.2e  sxy=%.2e  (E_oed=%.3f)\n",
                max_syy_err, max_sxx_err, max_sxy, oedometer_modulus());
}

} // namespace

int main() {
    test_confined_column_stress();
    if (g_failures == 0) {
        std::printf("OK: stress recovery matched the analytic solution\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
