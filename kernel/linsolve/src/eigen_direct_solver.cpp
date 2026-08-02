// Eigen-backed sparse direct solver: the portable backend.
//
// SimplicialLDLT is used for symmetric positive definite systems -- the elastic
// and SPD tangent case, and the effective-stiffness matrix of the Newmark scheme.
// SparseLU is used otherwise, because an LDLT of an indefinite matrix is not
// reliable and the nonsymmetric tangent (non-associated Mohr-Coulomb or Hardening
// Soil, Coulomb interfaces, embedded-beam skin plasticity) is not a candidate for
// it at all.
//
// Both Eigen solvers separate analyzePattern from factorize, so the pattern-reuse
// contract of the interface maps onto them directly.

#include <string>
#include <vector>

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <katai/linsolve/direct_solver.hpp>

namespace katai::linsolve {
namespace {

using SparseColMajor = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index>;

// CSR is row-major with sorted, deduplicated column indices, which is exactly an
// Eigen row-major sparse matrix; assigning it to a column-major matrix performs
// the transpose-copy Eigen's direct solvers want. Mapping avoids a first copy.
SparseColMajor to_eigen(const CsrMatrix& m) {
    const Eigen::Map<const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>> row_major(
        m.rows, m.cols, m.nonzeros(),
        m.row_ptr.data(), m.col_indices.data(), m.values.data());
    return SparseColMajor(row_major);
}

class EigenDirectSolver final : public DirectSolver {
public:
    explicit EigenDirectSolver(MatrixType type) : type_(type) {}

    void do_factorize(const CsrMatrix& matrix) override {
        if (!matrix.is_square())
            throw SolveError("linsolve: matrix is not square (" +
                             std::to_string(matrix.rows) + "x" + std::to_string(matrix.cols) + ")");
        if (matrix.rows == 0) throw SolveError("linsolve: matrix is empty");

        const SparseColMajor a = to_eigen(matrix);
        const bool same_pattern = pattern_matches(matrix);
        if (!same_pattern) remember_pattern(matrix);

        // info() is queried once, after factorize, for both solvers. It cannot be
        // queried between the two steps: SimplicialLDLT marks itself initialized in
        // analyzePattern, but SparseLU only does so in factorize, so an info() call
        // in between trips an assertion there -- which in a release build appears as
        // an unexplained fail-fast rather than a message.
        if (type_ == MatrixType::RealSymmetricPositiveDefinite) {
            if (!same_pattern || !analysed_) ldlt_.analyzePattern(a);
            ldlt_.factorize(a);
            if (ldlt_.info() != Eigen::Success)
                throw SingularSystem("linsolve/eigen: LDLT failed - the matrix is singular, not "
                                     "positive definite, or its pattern could not be analysed");
        } else {
            if (!same_pattern || !analysed_) lu_.analyzePattern(a);
            lu_.factorize(a);
            if (lu_.info() != Eigen::Success)
                throw SingularSystem("linsolve/eigen: LU failed - the matrix is singular or its "
                                     "pattern could not be analysed");
        }
        analysed_ = true;
        dimension_ = matrix.rows;
    }

    void do_solve(const Scalar* rhs, Scalar* solution) override {
        if (dimension_ <= 0) throw SolveError("linsolve: solve before factorize");
        const Eigen::Map<const Eigen::VectorXd> b(rhs, dimension_);
        Eigen::VectorXd x;
        if (type_ == MatrixType::RealSymmetricPositiveDefinite) {
            x = ldlt_.solve(b);
            if (ldlt_.info() != Eigen::Success)
                throw SingularSystem("linsolve/eigen: solve failed (LDLT)");
        } else {
            x = lu_.solve(b);
            if (lu_.info() != Eigen::Success)
                throw SingularSystem("linsolve/eigen: solve failed (LU)");
        }
        // Written after the solve so that rhs and solution may alias.
        Eigen::Map<Eigen::VectorXd>(solution, dimension_) = x;
    }

    Index rows() const override { return dimension_; }

private:
    bool pattern_matches(const CsrMatrix& m) const {
        return pattern_row_ptr_ == m.row_ptr && pattern_cols_ == m.col_indices;
    }
    void remember_pattern(const CsrMatrix& m) {
        pattern_row_ptr_ = m.row_ptr;
        pattern_cols_ = m.col_indices;
    }

    MatrixType type_;
    Index dimension_ = 0;
    bool analysed_ = false;

    Eigen::SimplicialLDLT<SparseColMajor> ldlt_;
    Eigen::SparseLU<SparseColMajor> lu_;

    std::vector<Index> pattern_row_ptr_;
    std::vector<Index> pattern_cols_;
};

} // namespace

std::unique_ptr<DirectSolver> make_eigen_solver(MatrixType type) {
    return std::make_unique<EigenDirectSolver>(type);
}

} // namespace katai::linsolve
