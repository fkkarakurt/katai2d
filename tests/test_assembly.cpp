// Pipeline validation (P0.6): mesh -> DOF -> assembly -> BC -> load -> solve.
//
// Laterally confined (oedometric) column: base clamped, side rollers (ux=0), top free.
// Lateral confinement → exx=0 → uniaxial compression. The closed-form solution is known
// for this case and the displacement field is represented EXACTLY in the tri6 polynomial
// space, so the FEM result must match the analytic one to round-off.
//
//   confined (oedometric) modulus:  E_oed = E (1-v) / ((1+v)(1-2v))
//   surface pressure p:             uy(y) = -p y / E_oed                (linear)
//   self-weight γ:                  uy(y) = -(γ/E_oed)(H y - y²/2)      (quadratic)
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/math/sparse_matrix.hpp>
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
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_failures;
    }
}

constexpr double kYoungs = 10000.0;
constexpr double kPoisson = 0.3;
constexpr double kHeight = 10.0;
constexpr double kWidth = 2.0;

double oedometer_modulus() {
    const double v = kPoisson;
    return kYoungs * (1.0 - v) / ((1.0 + v) * (1.0 - 2.0 * v));
}

// The confined-column DOF map: base ux+uy clamped, sides ux=0.
DofMap make_confined_dofs(const Mesh& mesh) {
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) {
        dofs.fix_node_component(n, 0);
        dofs.fix_node_component(n, 1);
    }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    return dofs;
}

Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k,
                          const Eigen::VectorXd& f) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(f);
}

// Compare the full displacement field against the analytic uy(y); ux≈0 expected.
void verify_field(const Mesh& mesh, const DofMap& dofs,
                  const Eigen::VectorXd& full,
                  double (*uy_exact)(double), const char* label) {
    const double reference = std::fabs(uy_exact(kHeight));  // top-displacement scale
    double max_uy_err = 0.0;
    double max_ux = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double ux = full(dofs.global_dof(n, 0));
        const double uy = full(dofs.global_dof(n, 1));
        max_ux = std::fmax(max_ux, std::fabs(ux));
        max_uy_err = std::fmax(max_uy_err, std::fabs(uy - uy_exact(mesh.y[n])));
    }
    check(max_ux < 1e-10 * reference, "yanal deplasman ux ≈ 0");
    check(max_uy_err < 1e-8 * reference, "the uy field matches the analytic one");
    std::printf("  %s: |uy_top|=%.6e  max_uy_err=%.2e  max_ux=%.2e\n", label,
                reference, max_uy_err, max_ux);
}

double g_pressure = 0.0;
double uy_pressure(double y) { return -g_pressure * y / oedometer_modulus(); }

double g_unit_weight = 0.0;
double uy_gravity(double y) {
    return -(g_unit_weight / oedometer_modulus()) * (kHeight * y - y * y / 2.0);
}

// Case A: uniform pressure p on the top surface.
void test_surface_pressure() {
    const RectangularDomain domain{0.0, 0.0, kWidth, kHeight, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 8);
    const DofMap dofs = make_confined_dofs(mesh);
    const std::vector<LinearElastic> materials = {{kYoungs, kPoisson}};

    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    const katai::math::CsrMatrix k = builder.build();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    g_pressure = 100.0;
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0,
                                           -g_pressure, f);

    const Eigen::VectorXd full =
        katai::core::expand_to_full(dofs, solve_spd(k, f));
    verify_field(mesh, dofs, full, uy_pressure, "surface pressure");
}

// Case B: self-weight (body force).
void test_self_weight() {
    const RectangularDomain domain{0.0, 0.0, kWidth, kHeight, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 8);
    const DofMap dofs = make_confined_dofs(mesh);
    const std::vector<LinearElastic> materials = {{kYoungs, kPoisson}};

    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    const katai::math::CsrMatrix k = builder.build();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    g_unit_weight = 20.0;
    katai::core::assemble_gravity(mesh, dofs, {g_unit_weight}, f);

    const Eigen::VectorXd full =
        katai::core::expand_to_full(dofs, solve_spd(k, f));
    verify_field(mesh, dofs, full, uy_gravity, "self-weight");
}

} // namespace

int main() {
    test_surface_pressure();
    test_self_weight();

    if (g_failures == 0) {
        std::printf("OK: the assembly pipeline matched the analytic solution\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
