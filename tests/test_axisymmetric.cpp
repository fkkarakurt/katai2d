// Axisymmetric element kernels (P1.6) verified in isolation, before global use.
// The defining new term is the hoop strain eps_theta = u_r / r (B's 4th row). Two
// uniform-strain fields exercise it on a distorted element, for BOTH tri6 and tri15:
//   - radial expansion u_r = c*r, u_z = 0  -> strain = [c, 0, 0, c]
//   - axial stretch    u_r = 0,   u_z = d*z -> strain = [0, d, 0, 0]
// (A linear u_r = c*r is the only radial field giving a constant hoop strain, since
// eps_theta = u_r/r; this is the axisymmetric constant-strain patch test.)
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/linear_elastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::LinearElastic;
using katai::core::Tri15Element;
using katai::core::Tri6Element;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

constexpr double kTri6Ref[6][2] = {
    {0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5}};
constexpr double kTri15Ref[15][2] = {
    {0.00, 0.00}, {1.00, 0.00}, {0.00, 1.00}, {0.25, 0.00}, {0.50, 0.00},
    {0.75, 0.00}, {0.75, 0.25}, {0.50, 0.50}, {0.25, 0.75}, {0.00, 0.75},
    {0.00, 0.50}, {0.00, 0.25}, {0.25, 0.25}, {0.50, 0.25}, {0.25, 0.50}};

template <class E>
void test_patch(const double ref[][2], const char* name) {
    // Distorted element in r > 0 (vertices in (r, z)); nodes at the lattice points.
    const double P[3][2] = {{1.0, 0.0}, {2.4, 0.4}, {1.3, 1.8}};
    typename E::NodeCoords nodes;
    for (int i = 0; i < E::kNodeCount; ++i) {
        const double l1 = 1.0 - ref[i][0] - ref[i][1], l2 = ref[i][0], l3 = ref[i][1];
        nodes(i, 0) = P[0][0] * l1 + P[1][0] * l2 + P[2][0] * l3;
        nodes(i, 1) = P[0][1] * l1 + P[1][1] * l2 + P[2][1] * l3;
    }

    // Field (a): radial expansion u_r = c*r -> strain [c, 0, 0, c].
    constexpr double c = 0.0007;
    Eigen::Matrix<double, E::kDofCount, 1> ua = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
    for (int i = 0; i < E::kNodeCount; ++i) ua(2 * i) = c * nodes(i, 0);  // u_r
    // Field (b): axial stretch u_z = d*z -> strain [0, d, 0, 0].
    constexpr double d = -0.0005;
    Eigen::Matrix<double, E::kDofCount, 1> ub = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
    for (int i = 0; i < E::kNodeCount; ++i) ub(2 * i + 1) = d * nodes(i, 1);  // u_z

    const Eigen::Vector4d ea(c, 0, 0, c), eb(0, d, 0, 0);
    bool ok_a = true, ok_b = true, pos_r = true;
    for (const auto& gp : E::gauss_points()) {
        const auto g = katai::core::axisym::strain_displacement<E>(nodes, gp.xi, gp.eta);
        if (g.radius <= 0.0 || g.det_jacobian <= 0.0) pos_r = false;
        if ((g.B * ua - ea).cwiseAbs().maxCoeff() > 1e-10) ok_a = false;
        if ((g.B * ub - eb).cwiseAbs().maxCoeff() > 1e-10) ok_b = false;
    }
    check(pos_r, "axisym: positive radius and Jacobian at Gauss points");
    check(ok_a, "axisym: radial expansion -> eps_r = eps_theta = c (hoop term)");
    check(ok_b, "axisym: axial stretch -> eps_z = d");

    // Stiffness is symmetric (axisymmetric D and B^T D B are symmetric).
    const LinearElastic mat{10000.0, 0.3};
    const auto ke = katai::core::axisym::element_stiffness<E>(nodes, mat.axisymmetric_matrix());
    const double asym = (ke - ke.transpose()).cwiseAbs().maxCoeff() / ke.cwiseAbs().maxCoeff();
    check(asym < 1e-12, "axisym: element stiffness is symmetric");
    std::printf("  [%s] axisym kernels ok (stiffness asym = %.1e)\n", name, asym);
}

} // namespace

int main() {
    test_patch<Tri6Element>(kTri6Ref, "tri6");
    test_patch<Tri15Element>(kTri15Ref, "tri15");
    if (g_failures == 0) {
        std::printf("OK: axisymmetric element kernels (P1.6) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
