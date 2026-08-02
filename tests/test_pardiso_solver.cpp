// PardisoSolver wrapper unit test — small SPD and nonsymmetric systems.
// CTest + assert style: 0 = passed, 1 = failed.
#include <katai/math/pardiso_solver.hpp>
#include <katai/math/sparse_matrix.hpp>

#include <cmath>
#include <cstdio>

#include <Eigen/Core>

using katai::math::CsrMatrix;
using katai::math::Index;
using katai::math::PardisoSolver;
using katai::math::Scalar;
using katai::math::SparseMatrixBuilder;

namespace {

int g_failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_failures;
    }
}

Scalar residual(const CsrMatrix& a, const Eigen::VectorXd& x,
                const Eigen::VectorXd& b) {
    return (a * x - b).norm();
}

// 1D Laplace: tridiag(-1, 2, -1), symmetric positive definite.
CsrMatrix make_laplacian(Index n) {
    SparseMatrixBuilder builder(n);
    for (Index i = 0; i < n; ++i) {
        builder.add_entry(i, i, 2.0);
        if (i > 0) builder.add_entry(i, i - 1, -1.0);
        if (i + 1 < n) builder.add_entry(i, i + 1, -1.0);
    }
    return builder.build();
}

// SPD system: factor once, solve again for two different right-hand sides.
void test_spd_solve_and_refactor_reuse() {
    const Index n = 5;
    const CsrMatrix a = make_laplacian(n);

    PardisoSolver solver(PardisoSolver::MatrixType::RealSymmetricPositiveDefinite);
    solver.factorize(a);
    check(solver.is_factorized(), "is_factorized after factorize");
    check(solver.rows() == n, "rows() = n");

    // 1) Known solution x = 1: b = A·1 = [1,0,0,0,1].
    Eigen::VectorXd x_expected = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b1 = a * x_expected;
    Eigen::VectorXd x1 = solver.solve(b1);
    check(residual(a, x1, b1) < 1e-10, "SPD rezidual (rhs #1)");
    check((x1 - x_expected).cwiseAbs().maxCoeff() < 1e-10, "SPD solution = the ones vector");

    // 2) Reuse the same factorization for a different right-hand side.
    Eigen::VectorXd x2_expected(n);
    x2_expected << 1, 2, 3, 4, 5;
    Eigen::VectorXd b2 = a * x2_expected;
    Eigen::VectorXd x2 = solver.solve(b2);
    check((x2 - x2_expected).cwiseAbs().maxCoeff() < 1e-9,
          "SPD solve reuse (rhs #2)");
}

// Nonsimetrik sistem (mtype = 11): a[2][1] != a[1][2].
void test_nonsymmetric_solve() {
    SparseMatrixBuilder builder(3);
    builder.add_entry(0, 0, 4.0);
    builder.add_entry(0, 1, 1.0);
    builder.add_entry(1, 0, 1.0);
    builder.add_entry(1, 1, 3.0);
    builder.add_entry(1, 2, 1.0);
    builder.add_entry(2, 1, 2.0);  // breaks symmetry (2, not 1)
    builder.add_entry(2, 2, 5.0);
    const CsrMatrix a = builder.build();

    Eigen::VectorXd x_expected(3);
    x_expected << 1, 1, 1;
    Eigen::VectorXd b = a * x_expected;  // = [5, 5, 7]

    PardisoSolver solver(PardisoSolver::MatrixType::RealNonsymmetric);
    solver.factorize(a);
    Eigen::VectorXd x = solver.solve(b);
    check(residual(a, x, b) < 1e-10, "nonsimetrik rezidual");
    check((x - x_expected).cwiseAbs().maxCoeff() < 1e-10, "nonsymmetric solve");
}

} // namespace

int main() {
    test_spd_solve_and_refactor_reuse();
    test_nonsymmetric_solve();

    if (g_failures == 0) {
        std::printf("OK: pardiso_solver passed all checks\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
