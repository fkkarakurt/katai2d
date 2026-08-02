#include <katai/math/katai_math.hpp>

#include <Eigen/Dense>

namespace katai::math {

const char* version() {
    return "katai_math 0.0.1";
}

double smoke_solve_residual() {
    Eigen::Matrix2d A;
    A << 4.0, 1.0,
         1.0, 3.0;
    Eigen::Vector2d b;
    b << 1.0, 2.0;

    const Eigen::Vector2d x = A.ldlt().solve(b);   // SPD → LDLᵀ
    return (A * x - b).norm();
}

} // namespace katai::math
