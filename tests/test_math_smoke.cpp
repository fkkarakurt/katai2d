// katai_math smoke test — verifies the Eigen solve really works.
#include <katai/math/katai_math.hpp>

#include <cmath>
#include <cstdio>

int main() {
    const double res = katai::math::smoke_solve_residual();
    if (!(std::abs(res) < 1e-12)) {
        std::fprintf(stderr, "FAIL: rezidual cok buyuk: %.3e\n", res);
        return 1;
    }
    std::printf("OK: %s, rezidual=%.3e\n", katai::math::version(), res);
    return 0;
}
