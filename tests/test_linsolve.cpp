// Sparse direct solver: the linked backend and the Eigen backend must agree with
// each other and with an independent dense factorization.
//
// Oracle discipline. The reference is Eigen's dense FullPivLU, which shares no
// code path with either sparse solver -- a different algorithm on a different
// storage layout. Cross-backend agreement alone would only prove the two are
// consistent, not that either is right, so both comparisons are made:
//
//   1. every backend against the dense reference   (is it correct?)
//   2. backend against backend                    (does the configuration change the answer?)
//
// The second is the mechanical form of the claim that a result never depends on
// which solver a program was built with.
//
// CTest style: 0 = pass, 1 = fail.
//
// verify: KV-NUM-002
//   oracle:   independent_path
//   source:   Eigen dense FullPivLU factorization -- a different algorithm on a different storage layout, sharing no code with either sparse backend
//   locator:  this file: SPD and nonsymmetric fixtures solved by every backend and by the dense reference
//   quantity: relative solution error against the dense reference, and cross-backend agreement [-]
//   expected: agreement to round-off; measured 2.66e-15 (SPD) and 4.44e-16 (nonsymmetric) vs dense, 4.44e-16 cross-backend
//   band:     1e-9 on max|dx| vs dense and cross-backend, 1e-12 for the in-place identity, as asserted below -- round-off headroom; a singular system must be refused by every backend (residual guard)

#include <katai/linsolve/direct_solver.hpp>
#include <katai/math/sparse_matrix.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

using katai::linsolve::DirectSolver;
using katai::linsolve::MatrixType;
using katai::linsolve::SolveError;
using katai::math::CsrMatrix;
using katai::math::Index;
using katai::math::Scalar;
using katai::math::SparseMatrixBuilder;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what.c_str());
        ++g_failures;
    }
}

Eigen::MatrixXd to_dense(const CsrMatrix& m) {
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(m.rows, m.cols);
    for (Index r = 0; r < m.rows; ++r)
        for (Index p = m.row_ptr[r]; p < m.row_ptr[r + 1]; ++p)
            d(r, m.col_indices[p]) = m.values[p];
    return d;
}

// A symmetric positive definite band matrix: the 1D stiffness pattern, which is
// what an elastic tangent looks like, plus a diagonal boost so it stays SPD.
CsrMatrix spd_system(Index n, Scalar diag_scale = 1.0) {
    SparseMatrixBuilder b(n, n);
    for (Index i = 0; i < n; ++i) {
        b.add_entry(i, i, diag_scale * (2.0 + 0.1 * i));
        if (i + 1 < n) {
            b.add_entry(i, i + 1, -1.0);
            b.add_entry(i + 1, i, -1.0);
        }
    }
    return b.build();
}

// Nonsymmetric, as a non-associated plastic tangent or a Coulomb interface gives.
CsrMatrix nonsymmetric_system(Index n) {
    SparseMatrixBuilder b(n, n);
    for (Index i = 0; i < n; ++i) {
        b.add_entry(i, i, 3.0 + 0.2 * i);
        if (i + 1 < n) b.add_entry(i, i + 1, -1.0);
        if (i > 0) b.add_entry(i, i - 1, -0.4);   // asymmetric off-diagonals
    }
    return b.build();
}

Eigen::VectorXd ramp(Index n) {
    Eigen::VectorXd f(n);
    for (Index i = 0; i < n; ++i) f[i] = 1.0 + 0.5 * i;
    return f;
}

// Solve with one backend and report the worst deviation from the dense reference.
double deviation_from_dense(DirectSolver& solver, const CsrMatrix& a, const Eigen::VectorXd& f) {
    solver.factorize(a);
    const Eigen::VectorXd x = solver.solve(f);
    const Eigen::VectorXd reference = to_dense(a).fullPivLu().solve(f);
    return (x - reference).cwiseAbs().maxCoeff();
}

void test_agrees_with_dense_reference() {
    struct Case {
        const char* name;
        MatrixType type;
        CsrMatrix a;
    };
    const std::vector<Case> cases = {
        {"SPD band, n=40", MatrixType::RealSymmetricPositiveDefinite, spd_system(40)},
        {"nonsymmetric band, n=40", MatrixType::RealNonsymmetric, nonsymmetric_system(40)},
    };

    for (const Case& c : cases) {
        const Eigen::VectorXd f = ramp(c.a.rows);
        auto linked = katai::linsolve::make_direct_solver(c.type);
        const double d_linked = deviation_from_dense(*linked, c.a, f);
        check(d_linked < 1e-9,
              std::string(c.name) + ": linked backend (" + katai::linsolve::backend_name() +
                  ") vs dense, max|dx| = " + std::to_string(d_linked));

        auto eigen = katai::linsolve::make_eigen_solver(c.type);
        const double d_eigen = deviation_from_dense(*eigen, c.a, f);
        check(d_eigen < 1e-9,
              std::string(c.name) + ": eigen backend vs dense, max|dx| = " + std::to_string(d_eigen));

        // Backend against backend: the configuration must not change the answer.
        auto a_solver = katai::linsolve::make_direct_solver(c.type);
        auto b_solver = katai::linsolve::make_eigen_solver(c.type);
        a_solver->factorize(c.a);
        b_solver->factorize(c.a);
        const double cross = (a_solver->solve(f) - b_solver->solve(f)).cwiseAbs().maxCoeff();
        check(cross < 1e-9, std::string(c.name) + ": backend-to-backend agreement, max|dx| = " +
                                std::to_string(cross));

        // The interface verifies every solve against the system; on a healthy one
        // the reported residual must be far below the tolerance, otherwise the
        // guard is protecting nothing and would fire on legitimate models.
        check(a_solver->last_relative_residual() < 1e-10,
              std::string(c.name) + ": reported residual is tight on a healthy system, " +
                  std::to_string(a_solver->last_relative_residual()));

        std::printf("  %-26s dense: linked %.2e  eigen %.2e   cross-backend %.2e   residual %.2e\n",
                    c.name, d_linked, d_eigen, cross, a_solver->last_relative_residual());
    }
}

// Newton iterations refactorize the same pattern with new values every step, so
// the second factorization must be correct while skipping symbolic analysis.
void test_pattern_reuse_refactorizes_correctly() {
    const CsrMatrix first = spd_system(30, 1.0);
    const CsrMatrix second = spd_system(30, 4.0);   // identical pattern, different values
    check(first.row_ptr == second.row_ptr && first.col_indices == second.col_indices,
          "pattern-reuse fixture: the two matrices share a pattern");

    const Eigen::VectorXd f = ramp(first.rows);
    auto solver = katai::linsolve::make_direct_solver(MatrixType::RealSymmetricPositiveDefinite);

    solver->factorize(first);
    const Eigen::VectorXd x1 = solver->solve(f);
    solver->factorize(second);                       // symbolic analysis skipped
    const Eigen::VectorXd x2 = solver->solve(f);

    const double d1 = (x1 - to_dense(first).fullPivLu().solve(f)).cwiseAbs().maxCoeff();
    const double d2 = (x2 - to_dense(second).fullPivLu().solve(f)).cwiseAbs().maxCoeff();
    check(d1 < 1e-9, "pattern reuse: first factorization correct, max|dx| = " + std::to_string(d1));
    check(d2 < 1e-9, "pattern reuse: refactorization correct, max|dx| = " + std::to_string(d2));

    // The values really did change, so a stale factorization would be visible.
    check((x1 - x2).cwiseAbs().maxCoeff() > 1e-6,
          "pattern reuse: the refactorized solution differs, so the check is not vacuous");
    std::printf("  pattern reuse              first %.2e  refactorized %.2e\n", d1, d2);
}

// rhs and solution may alias; the interface promises it and callers rely on it.
// Checked on both backends, because this is precisely where they behaved
// differently: Eigen writes the result after solving, while PARDISO reads the
// right-hand side as it writes, so it needs a scratch copy when the two overlap.
// Without that guard an aliased solve returned a wrong answer in silence.
void test_solve_in_place() {
    const CsrMatrix a = spd_system(20);
    const Eigen::VectorXd f = ramp(a.rows);

    struct Backend {
        const char* name;
        std::unique_ptr<DirectSolver> solver;
    };
    std::vector<Backend> backends;
    backends.push_back({katai::linsolve::backend_name(),
                        katai::linsolve::make_direct_solver(MatrixType::RealSymmetricPositiveDefinite)});
    backends.push_back({"eigen",
                        katai::linsolve::make_eigen_solver(MatrixType::RealSymmetricPositiveDefinite)});

    for (Backend& b : backends) {
        b.solver->factorize(a);
        const Eigen::VectorXd out_of_place = b.solver->solve(f);
        Eigen::VectorXd x = f;
        b.solver->solve(x.data(), x.data());        // in place
        const double d = (x - out_of_place).cwiseAbs().maxCoeff();
        check(d < 1e-12, std::string("in-place solve matches out-of-place on backend '") + b.name +
                             "', max|dx| = " + std::to_string(d));
    }
}

// A singular system must be refused, not answered: a plausible wrong number is
// worse than a failure, because it does not announce itself.
//
// The fixture verifies its own singularity with a dense rank computation. The
// first version of this test was not singular at all -- determinant -1 -- and so
// would have passed only by accident, or as here failed and pointed at the solver
// instead of at itself. A fixture that asserts the property it is built to
// exercise cannot rot into a vacuous test.
void test_singular_is_refused() {
    SparseMatrixBuilder b(4, 4);
    b.add_entry(0, 0, 1.0); b.add_entry(0, 1, 1.0);
    b.add_entry(1, 0, 1.0); b.add_entry(1, 1, 1.0);   // row 1 == row 0 -> rank deficient
    b.add_entry(2, 2, 1.0);
    b.add_entry(3, 3, 1.0);
    const CsrMatrix singular = b.build();

    const Eigen::FullPivLU<Eigen::MatrixXd> reference(to_dense(singular));
    check(reference.rank() < singular.rows,
          "singular fixture: the dense rank is deficient (" + std::to_string(reference.rank()) +
              " < " + std::to_string(singular.rows) + ")");

    for (const char* which : {"linked", "eigen"}) {
        bool refused = false;
        std::string how;
        try {
            auto solver = (std::string(which) == "eigen")
                              ? katai::linsolve::make_eigen_solver(MatrixType::RealNonsymmetric)
                              : katai::linsolve::make_direct_solver(MatrixType::RealNonsymmetric);
            solver->factorize(singular);
            const Eigen::VectorXd x = solver->solve(ramp(4));
            // A factorization may also "succeed" and return non-finite values; that
            // is still a detectable refusal, unlike a finite wrong answer.
            refused = !x.allFinite();
            how = refused ? "non-finite solution" : "returned a finite answer";
        } catch (const std::exception& e) {
            refused = true;
            how = std::string("threw: ") + e.what();
        }
        check(refused, std::string("singular system refused by the ") + which + " backend (" + how + ")");
        std::printf("  singular, %-6s backend  %s\n", which, how.c_str());
    }
}

void test_solve_before_factorize_is_refused() {
    bool refused = false;
    try {
        auto solver = katai::linsolve::make_eigen_solver(MatrixType::RealSymmetricPositiveDefinite);
        Eigen::VectorXd x(3);
        solver->solve(x.data(), x.data());
    } catch (const std::exception&) {
        refused = true;
    }
    check(refused, "solving before factorizing is refused");
}

} // namespace

int main() {
    // Unbuffered: a solver backend that dies inside a native library takes the
    // process down without unwinding, and a buffered report would be lost with it.
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("linsolve backend linked into this build: %s\n", katai::linsolve::backend_name());
    test_agrees_with_dense_reference();
    test_pattern_reuse_refactorizes_correctly();
    test_solve_in_place();
    test_singular_is_refused();
    test_solve_before_factorize_is_refused();

    if (g_failures == 0) std::printf("test_linsolve: OK\n");
    return g_failures == 0 ? 0 : 1;
}
