// Backend selection: Intel oneMKL PARDISO.
//
// Compiled instead of select_eigen.cpp when MKL is available. It adapts the
// existing PardisoSolver to the DirectSolver interface rather than reimplementing
// it: that wrapper already owns the handle, the control array and the
// pattern-reuse logic, and it is covered by its own tests.

#include <katai/linsolve/direct_solver.hpp>
#include <katai/math/pardiso_solver.hpp>

namespace katai::linsolve {
namespace {

math::PardisoSolver::MatrixType to_pardiso(MatrixType type) {
    switch (type) {
        case MatrixType::RealNonsymmetric:
            return math::PardisoSolver::MatrixType::RealNonsymmetric;
        case MatrixType::RealSymmetricIndefinite:
            return math::PardisoSolver::MatrixType::RealSymmetricIndefinite;
        case MatrixType::RealSymmetricPositiveDefinite:
            break;
    }
    return math::PardisoSolver::MatrixType::RealSymmetricPositiveDefinite;
}

class PardisoDirectSolver final : public DirectSolver {
public:
    explicit PardisoDirectSolver(MatrixType type) : impl_(to_pardiso(type)) {}

    void do_factorize(const CsrMatrix& matrix) override { impl_.factorize(matrix); }
    void do_solve(const Scalar* rhs, Scalar* solution) override { impl_.solve(rhs, solution); }
    Index rows() const override { return impl_.rows(); }

private:
    math::PardisoSolver impl_;
};

} // namespace

std::unique_ptr<DirectSolver> make_direct_solver(MatrixType type) {
    return std::make_unique<PardisoDirectSolver>(type);
}

const char* backend_name() { return "pardiso"; }

} // namespace katai::linsolve
