#pragma once
// PARDISO direct sparse solver wrapper (Intel oneMKL — Decision D6).
//
// Typical use (factor once, solve again for many right-hand sides):
//     PardisoSolver solver(PardisoSolver::MatrixType::RealSymmetricPositiveDefinite);
//     solver.factorize(K);                 // analysis + numerical factorization
//     Eigen::VectorXd u = solver.solve(f);  // reused for every load
//
// This header does not include <mkl.h>; the MKL details are encapsulated in the .cpp.
// PARDISO's internal handle and control array are stored opaquely.

#include <array>
#include <vector>

#include <Eigen/Core>

#include <katai/math/sparse_matrix.hpp>

namespace katai::math {

class PardisoSolver {
public:
    // The values are PARDISO's `mtype` codes directly.
    enum class MatrixType : int {
        RealNonsymmetric = 11,
        RealSymmetricIndefinite = -2,
        RealSymmetricPositiveDefinite = 2,
    };

    explicit PardisoSolver(
        MatrixType type = MatrixType::RealSymmetricPositiveDefinite);
    ~PardisoSolver();

    // The handle + internal factorization are uniquely owned → no copy/move.
    // (The solver is constructed where it is used; a move can be added later if needed.)
    PardisoSolver(const PardisoSolver&) = delete;
    PardisoSolver& operator=(const PardisoSolver&) = delete;
    PardisoSolver(PardisoSolver&&) = delete;
    PardisoSolver& operator=(PardisoSolver&&) = delete;

    // Analysis (reordering/symbolic) + numerical factorization. For a symmetric type
    // the matrix's upper triangle is extracted internally → the caller passes the full
    // matrix.
    //
    // PATTERN REUSE: if the same solver object is given a matrix whose sparsity PATTERN
    // (row_ptr + col_indices) is identical to the previous one, the symbolic analysis is
    // NOT repeated — only the numerical factorization (PARDISO phase 22) runs. Newton
    // iterations change only values on the same mesh/DOF pattern, so this is the common
    // path; if the pattern changes (new mesh, phase activation) the full analysis is
    // redone automatically. Transparent to the caller.
    void factorize(const CsrMatrix& matrix);

    // Solves the factored system for one right-hand side (Ax = b).
    // rhs and solution must be at least `rows()` long; they may alias the same array.
    void solve(const Scalar* rhs, Scalar* solution);
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs);

    Index rows() const { return dimension_; }
    bool is_factorized() const { return dimension_ > 0; }

    // Pivots PARDISO had to perturb in the last factorization (iparm(14)).
    // Diagnostic only: perturbation is routine numerical robustness and not
    // evidence of singularity -- measured, 3 of 119 tests in this suite perturb
    // pivots while producing correct validated results. Acceptability is judged
    // from the residual of the solution, not from this count.
    Index perturbed_pivots() const;

private:
    void release() noexcept;  // frees PARDISO's internal memory

    MatrixType type_;
    Index dimension_ = 0;

    // PARDISO's internal handle (64 pointers) + control array (iparm, 64 integers).
    // Kept as void*/Index so the header does not depend on <mkl.h>.
    std::array<void*, 64> handle_{};
    std::array<Index, 64> control_{};

    // The CSR handed to PARDISO (upper triangle for symmetric). The solve phase needs
    // these arrays too, so the solver owns them.
    std::vector<Index> row_ptr_;
    std::vector<Index> col_indices_;
    std::vector<Scalar> values_;

    // The pattern of the last factored FULL matrix (a copy of the input
    // row_ptr/col_indices) — compared on the next factorize call to decide whether the
    // symbolic analysis can be skipped.
    std::vector<Index> pattern_row_ptr_;
    std::vector<Index> pattern_cols_;
};

} // namespace katai::math
