#pragma once
// Sparse direct solver interface, with the backend chosen at link time.
//
// The engine solves Ax = b through this interface and never names a backend. Two
// implementations exist:
//
//   eigen    SimplicialLDLT for symmetric positive definite systems, SparseLU
//            otherwise. Always available, since Eigen is vendored.
//   pardiso  Intel oneMKL PARDISO, available when the MKL backend is linked.
//
// Which one a program contains is decided by CMake, by compiling exactly one
// selection translation unit. There is deliberately no preprocessor conditional
// in any consumer: an optional capability is a separate target, not an #ifdef,
// because conditional compilation scattered through shared code multiplies the
// number of programs that exist while testing only one of them.
//
// Typical use — factorize once, solve for many right-hand sides:
//
//     auto solver = katai::linsolve::make_direct_solver(
//         katai::linsolve::MatrixType::RealSymmetricPositiveDefinite);
//     solver->factorize(K);
//     Eigen::VectorXd u = solver->solve(f);
//
// Pattern reuse: handing the same sparsity pattern (row_ptr and col_indices) to
// the same solver object again skips the symbolic analysis and refactorizes
// numerically only. Newton iterations change values on a fixed mesh pattern, so
// that is the common path; a changed pattern is detected and re-analysed. This is
// transparent to the caller and both backends honour it.

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include <katai/math/solve_error.hpp>
#include <katai/math/sparse_matrix.hpp>

namespace katai::linsolve {

using math::CsrMatrix;
using math::Index;
using math::Scalar;

enum class MatrixType {
    RealNonsymmetric,
    RealSymmetricIndefinite,
    RealSymmetricPositiveDefinite,
};

// Both failure types come from katai::math, so katai::core can name them without
// depending on this module. SingularSystem means no usable solution exists and the
// caller may recover; SolveError means the solve could not be attempted at all. See
// katai/math/solve_error.hpp for why that is the axis the split runs along.
using math::SingularSystem;
using math::SolveError;

// Every solve is verified against the system it claims to solve.
//
// This is not belt-and-braces, it is the contract. A direct solver can complete a
// factorization and return a finite answer that satisfies nothing: PARDISO
// perturbs tiny pivots to stay robust and reports no error, so a rank-deficient
// matrix -- what an insufficiently restrained model produces -- yields a plausible
// displacement field. The perturbation count is not a usable signal either:
// measured on this suite, 3 of 119 tests perturb pivots while producing correct,
// validated results. What distinguishes a good answer from a bad one is whether it
// satisfies the equations, so that is what is checked, and it is checked
// identically for every backend.
//
// Cost: one sparse matrix-vector product per solve, against a factorization that
// is orders of magnitude more expensive, plus one copy of the matrix.
constexpr double kDefaultResidualTolerance = 1e-6;

class DirectSolver {
public:
    virtual ~DirectSolver() = default;

    DirectSolver(const DirectSolver&) = delete;
    DirectSolver& operator=(const DirectSolver&) = delete;

    // Symbolic analysis when the pattern is new, then numerical factorization.
    // The caller passes the full matrix; a backend that needs only one triangle
    // extracts it itself. Throws SolveError on failure.
    //
    // Not virtual on purpose: it keeps the system for verification and then calls
    // the backend, so a backend cannot be written that forgets to be verifiable.
    void factorize(const CsrMatrix& matrix) {
        do_factorize(matrix);
        system_ = matrix;
    }

    // Solve the factorized system for one right-hand side. `rhs` and `solution`
    // must be at least rows() long and may alias. Throws SolveError if the result
    // does not satisfy the system to within the residual tolerance.
    void solve(const Scalar* rhs, Scalar* solution) {
        std::vector<Scalar> b(rhs, rhs + rows());   // rhs may alias solution
        do_solve(rhs, solution);
        verify(b.data(), solution);
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) {
        Eigen::VectorXd x(rows());
        solve(rhs.data(), x.data());
        return x;
    }

    virtual Index rows() const = 0;
    bool is_factorized() const { return rows() > 0; }

    // Relative infinity-norm residual of the last solve, ||Ax-b|| / max(||b||,1).
    // Kept so a caller can report solution quality rather than only react to a
    // refusal, and so a benchmark can show what the tolerance is protecting.
    double last_relative_residual() const { return last_residual_; }

    void set_residual_tolerance(double tolerance) { tolerance_ = tolerance; }
    double residual_tolerance() const { return tolerance_; }

protected:
    DirectSolver() = default;

    virtual void do_factorize(const CsrMatrix& matrix) = 0;
    virtual void do_solve(const Scalar* rhs, Scalar* solution) = 0;

private:
    void verify(const Scalar* rhs, const Scalar* solution) {
        const Index n = rows();
        std::vector<Scalar> ax(static_cast<std::size_t>(n), 0.0);
        system_.multiply(solution, ax.data());

        double residual = 0.0, scale = 1.0;
        for (Index i = 0; i < n; ++i) {
            residual = std::max(residual, std::abs(ax[i] - rhs[i]));
            scale = std::max(scale, std::abs(rhs[i]));
        }
        last_residual_ = residual / scale;

        if (!(last_residual_ <= tolerance_)) {
            throw SingularSystem(
                "linsolve: the computed solution does not satisfy the system - relative residual " +
                std::to_string(last_residual_) + " exceeds " + std::to_string(tolerance_) +
                " (the matrix is singular or nearly so, commonly an insufficiently restrained model)");
        }
    }

    CsrMatrix system_;
    double tolerance_ = kDefaultResidualTolerance;
    double last_residual_ = 0.0;
};

// Create a solver of the linked backend. Never returns null: failure to build a
// solver is a programming error, not a run-time condition.
std::unique_ptr<DirectSolver> make_direct_solver(
    MatrixType type = MatrixType::RealSymmetricPositiveDefinite);

// The backend compiled into this program: "eigen" or "pardiso". Shown in the
// application's About panel today. It is deliberately a short stable token rather
// than prose, because it is meant to be written into result files and printed by
// the CLI so a number can be traced to the configuration that produced it -- both
// of those are still to come, and neither is claimed here as done.
const char* backend_name();

// Always available regardless of which backend was selected, so a test can
// compare the two on the same system rather than trusting either alone.
std::unique_ptr<DirectSolver> make_eigen_solver(MatrixType type);

} // namespace katai::linsolve
