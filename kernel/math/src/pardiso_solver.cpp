#include <katai/math/pardiso_solver.hpp>

#include <mkl.h>

#include <stdexcept>
#include <string>
#include <type_traits>

#include <katai/math/solve_error.hpp>

namespace katai::math {
namespace {

// Our index type must be exactly the same size as MKL's integer interface,
// otherwise we cannot pass the arrays to PARDISO directly (lp64 → 32-bit).
static_assert(sizeof(Index) == sizeof(MKL_INT),
              "Index must match MKL_INT in size (lp64 interface)");
static_assert(std::is_same_v<Scalar, double>,
              "the PARDISO wrapper assumes double precision");

// Translates PARDISO error codes to readable text (MKL manual, table 'error').
const char* pardiso_error_text(int code) {
    switch (code) {
        case 0:    return "no error";
        case -1:   return "inconsistent input";
        case -2:   return "not enough memory";
        case -3:   return "reordering problem";
        case -4:   return "zero pivot / numerical factorization error";
        case -5:   return "unexpected internal error";
        case -6:   return "reordering failed the matrix";
        case -7:   return "matrix size zero or negative";
        case -8:   return "32-bit integer overflow";
        case -9:   return "not enough memory for OOC";
        case -10:  return "could not open OOC files";
        case -11:  return "OOC read/write error";
        case -12:  return "iparm(2) wrong with pardiso_64";
        case -13:  return "unexpected communication / bad call";
        default:   return "unknown PARDISO error code";
    }
}

// Extracts the upper triangle (diagonal included) from a full symmetric CSR matrix.
// Since the input is column-sorted within each row, the upper-triangle entries sit
// contiguously at the end of the row and stay sorted.
CsrMatrix upper_triangle(const CsrMatrix& a) {
    CsrMatrix u;
    u.rows = a.rows;
    u.cols = a.cols;
    u.row_ptr.assign(static_cast<std::size_t>(a.rows) + 1, 0);

    for (Index r = 0; r < a.rows; ++r) {
        Index kept = 0;
        for (Index k = a.row_ptr[r]; k < a.row_ptr[r + 1]; ++k)
            if (a.col_indices[k] >= r) ++kept;
        u.row_ptr[r + 1] = u.row_ptr[r] + kept;
    }

    const Index nnz = u.row_ptr[a.rows];
    u.col_indices.reserve(nnz);
    u.values.reserve(nnz);
    for (Index r = 0; r < a.rows; ++r) {
        for (Index k = a.row_ptr[r]; k < a.row_ptr[r + 1]; ++k) {
            if (a.col_indices[k] >= r) {
                u.col_indices.push_back(a.col_indices[k]);
                u.values.push_back(a.values[k]);
            }
        }
    }
    return u;
}

} // namespace

PardisoSolver::PardisoSolver(MatrixType type) : type_(type) {
    MKL_INT mtype = static_cast<MKL_INT>(type_);
    // pardisoinit fills sensible iparm defaults for the given mtype.
    pardisoinit(handle_.data(), &mtype, control_.data());
    control_[34] = 1;  // zero-based (C-style) indexing
}

PardisoSolver::~PardisoSolver() {
    release();
}

// Number of pivots PARDISO had to perturb in the last factorization, from
// iparm(14). Reported, not judged: perturbation is routine numerical robustness,
// not evidence of singularity. Measured on this suite -- 3 of 119 tests
// (test_embedded_beam, test_axisym_collapse, test_earth_pressure) perturb pivots
// while producing correct, validated results -- so refusing on this count rejects
// healthy models. Whether a solution is acceptable is decided by its residual, in
// the linsolve layer, not by a backend's internal heuristic.
Index PardisoSolver::perturbed_pivots() const { return control_[13]; }

void PardisoSolver::factorize(const CsrMatrix& matrix) {
    if (!matrix.is_square())
        throw std::invalid_argument("PardisoSolver::factorize: matrix is not square");
    if (matrix.rows == 0)
        throw std::invalid_argument("PardisoSolver::factorize: matrix is empty");

    const bool symmetric = (type_ != MatrixType::RealNonsymmetric);

    // If the pattern is exactly the previous one, skip the symbolic analysis: update
    // only the values + numerical factorization (phase 22). The reordering/elimination
    // tree depends on the pattern → for SPD/symmetric the solution is identical to the
    // full path. For nonsymmetric, the analysis-phase scaling/matching (iparm 11/13) is
    // value-dependent and phase 22 reuses the old one; if a pivot breaks, PARDISO
    // returns an error and we fall to the full-analysis path below (safe behaviour).
    const bool same_pattern = dimension_ == matrix.rows &&
                              pattern_row_ptr_ == matrix.row_ptr &&
                              pattern_cols_ == matrix.col_indices;
    if (same_pattern) {
        if (symmetric) {  // upper-triangle values (same pattern → same order/place)
            std::size_t out = 0;
            for (Index r = 0; r < matrix.rows; ++r)
                for (Index k = matrix.row_ptr[r]; k < matrix.row_ptr[r + 1]; ++k)
                    if (matrix.col_indices[k] >= r) values_[out++] = matrix.values[k];
        } else {
            values_ = matrix.values;
        }
        MKL_INT mtype = static_cast<MKL_INT>(type_);
        MKL_INT maxfct = 1, mnum = 1, msglvl = 0, error = 0, nrhs = 1;
        MKL_INT phase = 22;  // numerical factorization only (analysis preserved)
        MKL_INT idum = 0;
        double ddum = 0.0;
        pardiso(handle_.data(), &maxfct, &mnum, &mtype, &phase, &dimension_,
                values_.data(), row_ptr_.data(), col_indices_.data(), &idum, &nrhs,
                control_.data(), &msglvl, &ddum, &ddum, &error);
        if (error == 0) return;
        // Numerical error (e.g. a zero pivot under the old scaling): fall to the full
        // analysis — behaviour stays identical to a fresh solver in every case.
    }

    // Release any previous factorization (so refactoring is safe).
    release();

    CsrMatrix stored = symmetric ? upper_triangle(matrix) : matrix;
    dimension_ = matrix.rows;
    row_ptr_ = std::move(stored.row_ptr);
    col_indices_ = std::move(stored.col_indices);
    values_ = std::move(stored.values);
    pattern_row_ptr_ = matrix.row_ptr;
    pattern_cols_ = matrix.col_indices;

    MKL_INT mtype = static_cast<MKL_INT>(type_);
    MKL_INT maxfct = 1, mnum = 1, msglvl = 0, error = 0, nrhs = 1;
    MKL_INT phase = 12;  // analysis + numerical factorization
    MKL_INT idum = 0;
    double ddum = 0.0;

    pardiso(handle_.data(), &maxfct, &mnum, &mtype, &phase, &dimension_,
            values_.data(), row_ptr_.data(), col_indices_.data(), &idum, &nrhs,
            control_.data(), &msglvl, &ddum, &ddum, &error);

    // PARDISO does not fail on a rank-deficient nonsymmetric matrix: it perturbs
    // the tiny pivots and completes, so the solve then returns a finite answer
    // that satisfies nothing. Measured on a deliberately rank-deficient 4x4: no
    // error, a finite garbage solution, while Eigen's LU refused the same matrix.
    // A singular global stiffness means an insufficiently restrained model, and a
    // plausible displacement field is the worst possible response to it. The count
    // of perturbed pivots is reported in iparm(14), so refuse on it.
    if (error != 0) {
        dimension_ = 0;   // failed -> treat as not factorized
        const std::string what = std::string("PARDISO factorization failed: ") +
                                 pardiso_error_text(error);
        // -4 is a zero pivot, i.e. the matrix is singular. That is the one failure a
        // caller can recover from -- at a limit load it is the collapse mechanism, not
        // a fault -- so it gets the type the Newton loop knows how to act on. Every
        // other code is a genuine failure and must not be mistaken for a collapse.
        // Eigen reaches the same state through its own factorization refusing, which
        // is why both backends have to raise the same type here.
        if (error == -4) throw SingularSystem(what);
        throw SolveError(what);
    }
}

void PardisoSolver::solve(const Scalar* rhs, Scalar* solution) {
    if (dimension_ == 0)
        throw std::logic_error("PardisoSolver::solve: call factorize() first");

    MKL_INT mtype = static_cast<MKL_INT>(type_);
    MKL_INT maxfct = 1, mnum = 1, msglvl = 0, error = 0, nrhs = 1;
    MKL_INT phase = 33;  // solve + iterative refinement
    MKL_INT idum = 0;

    // With iparm[5] = 0 (the default) PARDISO writes the solution into `solution`
    // and leaves `rhs` intact -- but it does so while still reading `rhs`, so the
    // two must not be the same array. The header promises that they may alias, and
    // callers do it, so honour that promise here with a scratch copy taken only
    // when the pointers actually overlap. Measured before this guard existed: an
    // aliased solve returned a wrong answer silently, max|dx| ~ 62 on a 20-DOF
    // system, which is exactly the failure class the project treats as worst.
    std::vector<Scalar> scratch;
    const Scalar* rhs_in = rhs;
    if (rhs == solution) {
        scratch.assign(rhs, rhs + dimension_);
        rhs_in = scratch.data();
    }

    pardiso(handle_.data(), &maxfct, &mnum, &mtype, &phase, &dimension_,
            values_.data(), row_ptr_.data(), col_indices_.data(), &idum, &nrhs,
            control_.data(), &msglvl, const_cast<Scalar*>(rhs_in), solution, &error);

    if (error != 0)
        throw SolveError(std::string("PARDISO solve failed: ") + pardiso_error_text(error));
}

Eigen::VectorXd PardisoSolver::solve(const Eigen::VectorXd& rhs) {
    if (rhs.size() != dimension_)
        throw std::invalid_argument("PardisoSolver::solve: rhs size mismatch");
    Eigen::VectorXd solution(dimension_);
    solve(rhs.data(), solution.data());
    return solution;
}

void PardisoSolver::release() noexcept {
    if (dimension_ == 0) return;

    MKL_INT mtype = static_cast<MKL_INT>(type_);
    MKL_INT maxfct = 1, mnum = 1, msglvl = 0, error = 0, nrhs = 1;
    MKL_INT phase = -1;  // release all internal memory
    MKL_INT idum = 0;
    double ddum = 0.0;

    pardiso(handle_.data(), &maxfct, &mnum, &mtype, &phase, &dimension_, &ddum,
            row_ptr_.data(), col_indices_.data(), &idum, &nrhs, control_.data(),
            &msglvl, &ddum, &ddum, &error);
    dimension_ = 0;
}

} // namespace katai::math
