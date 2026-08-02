// Backend selection: Eigen.
//
// CMake compiles this translation unit when no accelerated backend is available,
// and compiles select_pardiso.cpp instead when Intel oneMKL is. Exactly one of
// them defines make_direct_solver and backend_name, so the choice is a link-time
// fact with no preprocessor conditional anywhere in the tree.

#include <katai/linsolve/direct_solver.hpp>

namespace katai::linsolve {

std::unique_ptr<DirectSolver> make_direct_solver(MatrixType type) {
    return make_eigen_solver(type);
}

const char* backend_name() { return "eigen"; }

} // namespace katai::linsolve
