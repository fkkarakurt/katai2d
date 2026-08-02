// tri6 element + linear elastic plane-strain unit tests.
// Verified: partition of unity, area integration, patch test (constant strain),
// Ke symmetry, rigid-body modes (zero energy), exactly 3 zero eigenvalues of Ke
// (2 translations + 1 rotation).
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/linear_elastic.hpp>

#include <cmath>
#include <cstdio>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using katai::core::LinearElastic;
namespace tri6 = katai::core::tri6;

namespace {

int g_failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_failures;
    }
}

// Straight-edged test triangle; mid-edge nodes at the true edge midpoints.
//   1=(0,0) 2=(2,0) 3=(0,3)  →  alan = 3
tri6::NodeCoords make_test_triangle() {
    tri6::NodeCoords nodes;
    nodes.row(0) << 0.0, 0.0;  // 1
    nodes.row(1) << 2.0, 0.0;  // 2
    nodes.row(2) << 0.0, 3.0;  // 3
    nodes.row(3) << 1.0, 0.0;  // 4 (1-2 orta)
    nodes.row(4) << 1.0, 1.5;  // 5 (2-3 orta)
    nodes.row(5) << 0.0, 1.5;  // 6 (3-1 orta)
    return nodes;
}

// Shape functions: Σ N_i = 1, Σ dN = 0 at every point.
void test_partition_of_unity() {
    const double pts[][2] = {{0.2, 0.3}, {1.0 / 6, 1.0 / 6}, {0.0, 0.0}};
    for (const auto& p : pts) {
        const auto n = tri6::shape_functions(p[0], p[1]);
        check(std::fabs(n.sum() - 1.0) < 1e-12, "Σ N_i = 1");
        const auto dn = tri6::shape_derivatives_natural(p[0], p[1]);
        check(dn.col(0).sum() < 1e-12 && dn.col(1).sum() < 1e-12, "Σ dN = 0");
    }
}

// Gauss integration must reproduce the element area: Σ w_g detJ_g = area.
void test_area_integration() {
    const auto nodes = make_test_triangle();
    double area = 0.0;
    for (const auto& gp : tri6::gauss_points()) {
        const auto g = tri6::strain_displacement(nodes, gp.xi, gp.eta);
        area += gp.weight * g.det_jacobian;
    }
    check(std::fabs(area - 3.0) < 1e-12, "integrasyonla alan = 3");
}

// Patch test: a linear displacement field → constant strain (B correctness).
//   u = a1 x + a2 y,  v = b1 x + b2 y
//   exx=a1, eyy=b2, gxy=a2+b1
void test_constant_strain_patch() {
    const auto nodes = make_test_triangle();
    const double a1 = 1e-3, a2 = 2e-3, b1 = 3e-3, b2 = 4e-3;

    Eigen::Matrix<double, 12, 1> u;
    for (int i = 0; i < tri6::kNodeCount; ++i) {
        const double x = nodes(i, 0), y = nodes(i, 1);
        u(2 * i) = a1 * x + a2 * y;       // ux
        u(2 * i + 1) = b1 * x + b2 * y;   // uy
    }
    const Eigen::Vector3d expected(a1, b2, a2 + b1);

    for (const auto& gp : tri6::gauss_points()) {
        const auto g = tri6::strain_displacement(nodes, gp.xi, gp.eta);
        const Eigen::Vector3d strain = g.B * u;
        check((strain - expected).cwiseAbs().maxCoeff() < 1e-12,
              "constant strain (patch)");
    }
}

// Ke must be symmetric.
void test_stiffness_symmetry(const tri6::ElementMatrix& ke) {
    check((ke - ke.transpose()).cwiseAbs().maxCoeff() < 1e-8 * ke.cwiseAbs().maxCoeff(),
          "Ke simetrik");
}

// Rigid-body motion → zero nodal force (zero strain energy).
void test_rigid_body_modes(const tri6::NodeCoords& nodes,
                           const tri6::ElementMatrix& ke) {
    const double rel = ke.norm();

    // Translation: all nodes (1, 0) then (0, 1).
    for (int axis = 0; axis < 2; ++axis) {
        Eigen::Matrix<double, 12, 1> u = Eigen::Matrix<double, 12, 1>::Zero();
        for (int i = 0; i < tri6::kNodeCount; ++i) u(2 * i + axis) = 1.0;
        check((ke * u).norm() < 1e-10 * rel, "translation → zero force");
    }

    // Infinitesimal rotation: u = -θ y, v = θ x.
    const double theta = 1e-3;
    Eigen::Matrix<double, 12, 1> u = Eigen::Matrix<double, 12, 1>::Zero();
    for (int i = 0; i < tri6::kNodeCount; ++i) {
        u(2 * i) = -theta * nodes(i, 1);
        u(2 * i + 1) = theta * nodes(i, 0);
    }
    check((ke * u).norm() < 1e-10 * rel * u.norm(), "rotation → zero force");
}

// Ke must have exactly 3 zero eigenvalues (2 translations + 1 rotation); the other 9 positive.
void test_zero_energy_modes(const tri6::ElementMatrix& ke) {
    Eigen::SelfAdjointEigenSolver<tri6::ElementMatrix> solver(ke);
    const auto eigenvalues = solver.eigenvalues();
    const double max_eig = eigenvalues.cwiseAbs().maxCoeff();
    const double tol = 1e-9 * max_eig;

    int zeros = 0;
    for (int i = 0; i < eigenvalues.size(); ++i)
        if (std::fabs(eigenvalues(i)) < tol) ++zeros;

    check(zeros == 3, "exactly 3 zero-energy modes");
    check(eigenvalues(eigenvalues.size() - 1) > 0.0, "largest eigenvalue positive");
}

} // namespace

int main() {
    test_partition_of_unity();
    test_area_integration();
    test_constant_strain_patch();

    const auto nodes = make_test_triangle();
    const LinearElastic material{/*E=*/10000.0, /*v=*/0.25};
    const auto ke = tri6::element_stiffness(nodes, material.plane_strain_matrix());

    test_stiffness_symmetry(ke);
    test_rigid_body_modes(nodes, ke);
    test_zero_energy_modes(ke);

    if (g_failures == 0) {
        std::printf("OK: tri6 passed all checks\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
