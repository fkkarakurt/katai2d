// Verification of the 15-node quartic triangle (tri15, P1.5) in isolation, before
// it is wired into the assembly/solver pipeline. Each property is checked against
// a closed-form expectation so any error in the (generically generated) quartic
// shape functions, their derivatives, or the 12-point Gauss rule is caught early.
#include <katai/fem/elements/tri15.hpp>

#include <array>
#include <cmath>
#include <cstdio>

using katai::core::tri15::ElementMatrix;
using katai::core::tri15::gauss_points;
using katai::core::tri15::kNodeCount;
using katai::core::tri15::NodeCoords;
using katai::core::tri15::shape_derivatives_natural;
using katai::core::tri15::shape_functions;
using katai::core::tri15::strain_displacement;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol = 1e-12) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

// Reference (xi, eta) coordinates of the 15 nodes (see tri15.hpp node ordering).
constexpr double kRef[15][2] = {
    {0.00, 0.00}, {1.00, 0.00}, {0.00, 1.00},
    {0.25, 0.00}, {0.50, 0.00}, {0.75, 0.00},
    {0.75, 0.25}, {0.50, 0.50}, {0.25, 0.75},
    {0.00, 0.75}, {0.00, 0.50}, {0.00, 0.25},
    {0.25, 0.25}, {0.50, 0.25}, {0.25, 0.50}};

// Shape functions: partition of unity, the Kronecker-delta (interpolation)
// property N_i(node_j) = delta_ij, and the derivative partition Sum dN = 0.
void test_shape_functions() {
    const double pts[][2] = {{0.2, 0.3}, {0.5, 0.1}, {0.05, 0.6}, {0.33, 0.33}};
    for (const auto& p : pts) {
        const auto n = shape_functions(p[0], p[1]);
        check(close(n.sum(), 1.0), "tri15: partition of unity (sum N = 1)");
        const auto d = shape_derivatives_natural(p[0], p[1]);
        check(close(d.col(0).sum(), 0.0) && close(d.col(1).sum(), 0.0),
              "tri15: derivative partition (sum dN = 0)");
    }
    for (int j = 0; j < kNodeCount; ++j) {
        const auto n = shape_functions(kRef[j][0], kRef[j][1]);
        bool delta = true;
        for (int i = 0; i < kNodeCount; ++i)
            if (!close(n(i), i == j ? 1.0 : 0.0, 1e-10)) delta = false;
        check(delta, "tri15: N_i(node_j) = delta_ij");
    }
}

// The 12-point Dunavant rule must integrate every monomial xi^i eta^j with
// i + j <= 6 exactly: Integral over the reference triangle = i! j! / (i+j+2)!.
void test_gauss_exactness() {
    auto fact = [](int k) { double f = 1.0; for (int i = 2; i <= k; ++i) f *= i; return f; };
    const auto gp = gauss_points();
    double wsum = 0.0;
    for (const auto& g : gp) wsum += g.weight;
    check(close(wsum, 0.5), "tri15: Gauss weights sum to 1/2 (triangle area)");

    bool exact = true;
    double worst = 0.0;
    for (int i = 0; i + 0 <= 6; ++i)
        for (int j = 0; i + j <= 6; ++j) {
            double q = 0.0;
            for (const auto& g : gp) q += g.weight * std::pow(g.xi, i) * std::pow(g.eta, j);
            const double exactv = fact(i) * fact(j) / fact(i + j + 2);
            worst = std::max(worst, std::fabs(q - exactv));
            if (std::fabs(q - exactv) > 1e-12) exact = false;
        }
    std::printf("  [tri15] Gauss degree-6 monomial worst error = %.2e\n", worst);
    check(exact, "tri15: 12-point rule integrates all monomials up to degree 6");
}

// Build a straight-sided element from three vertices: node i at the lattice point
// P0*L1 + P1*L2 + P2*L3, (L1,L2,L3) = (1-xi-eta, xi, eta) at the node's reference
// coordinate. The isoparametric map is then affine (constant Jacobian).
NodeCoords element_from(double P[3][2]) {
    NodeCoords nodes;
    for (int i = 0; i < kNodeCount; ++i) {
        const double xi = kRef[i][0], eta = kRef[i][1];
        const double L1 = 1.0 - xi - eta, L2 = xi, L3 = eta;
        nodes(i, 0) = P[0][0] * L1 + P[1][0] * L2 + P[2][0] * L3;
        nodes(i, 1) = P[0][1] * L1 + P[1][1] * L2 + P[2][1] * L3;
    }
    return nodes;
}

// Patch tests on a distorted (but straight-sided) element.
void test_patch() {
    double P[3][2] = {{0.0, 0.0}, {2.0, 0.3}, {0.4, 1.7}};
    const NodeCoords nodes = element_from(P);

    // (a) Constant-strain patch test: a linear displacement field must yield the
    //     exact constant strain at every Gauss point (FE completeness / convergence
    //     requirement). u = 2 + 0.5x - 0.3y, v = -1 + 0.2x + 0.7y.
    Eigen::Matrix<double, 30, 1> d;
    for (int i = 0; i < kNodeCount; ++i) {
        const double x = nodes(i, 0), y = nodes(i, 1);
        d(2 * i) = 2.0 + 0.5 * x - 0.3 * y;
        d(2 * i + 1) = -1.0 + 0.2 * x + 0.7 * y;
    }
    const Eigen::Vector3d expected(0.5, 0.7, -0.3 + 0.2);  // [exx, eyy, gxy]
    bool ok = true;
    for (const auto& g : gauss_points()) {
        const auto sg = strain_displacement(nodes, g.xi, g.eta);
        if ((sg.B * d - expected).cwiseAbs().maxCoeff() > 1e-10) ok = false;
        if (sg.det_jacobian <= 0.0) ok = false;  // positive area
    }
    check(ok, "tri15: constant-strain patch test (exact constant strain)");

    // (b) Completeness to degree 4: a quartic physical field sampled at the nodes
    //     is reproduced exactly by the interpolation at arbitrary points.
    auto quartic = [](double x, double y) {
        return 1.0 + 2.0 * x - y + 0.5 * x * y + x * x - 0.3 * y * y +
               x * x * y - 0.2 * y * y * y + 0.1 * x * x * y * y;
    };
    Eigen::Matrix<double, 15, 1> f;
    for (int i = 0; i < kNodeCount; ++i) f(i) = quartic(nodes(i, 0), nodes(i, 1));
    bool reproduced = true;
    const double q[][2] = {{0.2, 0.3}, {0.5, 0.2}, {0.1, 0.05}, {0.3, 0.3}};
    for (const auto& r : q) {
        const auto n = shape_functions(r[0], r[1]);
        const double xi = r[0], eta = r[1], L1 = 1.0 - xi - eta;
        const double x = P[0][0] * L1 + P[1][0] * xi + P[2][0] * eta;
        const double y = P[0][1] * L1 + P[1][1] * xi + P[2][1] * eta;
        if (!close(n.dot(f), quartic(x, y), 1e-10)) reproduced = false;
    }
    check(reproduced, "tri15: reproduces a quartic field exactly (completeness deg 4)");
}

} // namespace

int main() {
    test_shape_functions();
    test_gauss_exactness();
    test_patch();
    if (g_failures == 0) {
        std::printf("OK: tri15 quartic element kernels (P1.5) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
