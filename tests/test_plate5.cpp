// 5-node quartic Timoshenko plate element verification (P2.4, tri15 parity). Same checks as
// test_plate but for the 5-node element (sits on a tri15 edge). A quartic element captures the
// cubic tip-load cantilever almost exactly, so very few elements suffice. Closed form:
//   delta = P L^3/(3 EI) + P L/(k G A')   (bending + Timoshenko shear)
// (See docs/references/structural-plate-formulation.md.)
#include <katai/fem/elements/plate.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

using katai::core::plate::PlateProps;
using katai::core::plate::NodeCoords5;

// Cantilever of length L (ne 5-node elements) at angle theta, clamped at node 0, tip load +
// moment. Returns tip nodal (ux, uy, phi) in global axes.
Eigen::Vector3d cantilever(double L, int ne, const PlateProps& p, double theta,
                           double fx_tip, double fy_tip, double m_tip) {
    const int nnodes = 4 * ne + 1;
    const int ndof = 3 * nnodes;
    const double c = std::cos(theta), s = std::sin(theta);
    const double h = L / (4 * ne);  // node spacing (quarter nodes)
    std::vector<double> X(nnodes), Y(nnodes);
    for (int i = 0; i < nnodes; ++i) { X[i] = c * (i * h); Y[i] = s * (i * h); }

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
    for (int e = 0; e < ne; ++e) {
        int gn[5];
        for (int k = 0; k < 5; ++k) gn[k] = 4 * e + k;  // natural order xi=-1,-0.5,0,0.5,1
        NodeCoords5 Xe;
        for (int k = 0; k < 5; ++k) { Xe(k, 0) = X[gn[k]]; Xe(k, 1) = Y[gn[k]]; }
        const auto Ke = katai::core::plate::stiffness5(Xe, p);
        for (int a = 0; a < 5; ++a)
            for (int ca = 0; ca < 3; ++ca)
                for (int b = 0; b < 5; ++b)
                    for (int cb = 0; cb < 3; ++cb)
                        K(3 * gn[a] + ca, 3 * gn[b] + cb) += Ke(3 * a + ca, 3 * b + cb);
    }

    Eigen::VectorXd f = Eigen::VectorXd::Zero(ndof);
    const int tip = nnodes - 1;
    f(3 * tip + 0) = fx_tip; f(3 * tip + 1) = fy_tip; f(3 * tip + 2) = m_tip;

    std::vector<char> isfixed(ndof, 0);
    isfixed[0] = isfixed[1] = isfixed[2] = 1;  // clamp node 0
    std::vector<int> freed;
    for (int d = 0; d < ndof; ++d) if (!isfixed[d]) freed.push_back(d);
    const int nf = static_cast<int>(freed.size());
    Eigen::MatrixXd Kff(nf, nf); Eigen::VectorXd ff(nf);
    for (int i = 0; i < nf; ++i) {
        ff(i) = f(freed[i]);
        for (int j = 0; j < nf; ++j) Kff(i, j) = K(freed[i], freed[j]);
    }
    const Eigen::VectorXd uf = Kff.ldlt().solve(ff);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof);
    for (int i = 0; i < nf; ++i) u(freed[i]) = uf(i);
    return {u(3 * tip + 0), u(3 * tip + 1), u(3 * tip + 2)};
}

void test_cantilever() {
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0, P = 100.0;
    const double db = P * L * L * L / (3.0 * p.EI), ds = P * L / p.shear_rigidity();
    const double delta = db + ds;

    // Quartic element captures the cubic cantilever almost exactly -> ne=2 is plenty.
    const Eigen::Vector3d u = cantilever(L, 2, p, 0.0, 0.0, -P, 0.0);
    std::printf("  cantilever (ne=2): tip_defl=%.6e  analytic=%.6e (bend=%.4e shear=%.4e)\n",
                -u(1), delta, db, ds);
    check(close(-u(1), delta, 1e-4), "tip deflection = PL^3/3EI + PL/kGA' (quartic, ne=2)");

    // Single element -> cubic bending is exact; expect tight agreement already.
    const Eigen::Vector3d u1 = cantilever(L, 1, p, 0.0, 0.0, -P, 0.0);
    std::printf("  cantilever (ne=1): tip_defl=%.6e  (rel err %.2e)\n",
                -u1(1), std::fabs(-u1(1) - delta) / delta);
    check(close(-u1(1), delta, 1e-6), "single quartic element exact for tip-load cantilever");

    // Rotation invariance at 30 deg.
    const double th = 30.0 * kPi / 180.0;
    const Eigen::Vector3d u2 = cantilever(L, 2, p, th, P * std::sin(th), -P * std::cos(th), 0.0);
    const double tdefl = -u2(0) * std::sin(th) + u2(1) * std::cos(th);
    check(close(std::fabs(tdefl), delta, 1e-4), "rotation invariance (30deg)");
}

void test_axial_moment() {
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0;
    const double N = 500.0;
    const Eigen::Vector3d ua = cantilever(L, 2, p, 0.0, N, 0.0, 0.0);
    check(close(ua(0), N * L / p.EA, 1e-8), "pure axial elongation = NL/EA (quartic)");
    const double M = 200.0;
    const Eigen::Vector3d um = cantilever(L, 2, p, 0.0, 0.0, 0.0, M);
    check(close(um(2), M * L / p.EI, 1e-6), "tip rotation = ML/EI");
    check(close(um(1), M * L * L / (2.0 * p.EI), 1e-6), "tip deflection = ML^2/2EI");

    // Recover M at the clamp via forces5: M should equal the applied tip moment (pure bending).
    // Build a 1-element beam to read forces directly.
    NodeCoords5 Xe; for (int k = 0; k < 5; ++k) { Xe(k, 0) = (L / 4.0) * k; Xe(k, 1) = 0.0; }
    // (forces5 exercised in the wall benchmark; here just confirm the displacement fields above.)
    (void)Xe;
}

void test_thin_no_locking() {
    const double E = 3.0e7, d = 0.02, nu = 0.15, L = 10.0, P = 1.0;  // d/L = 1/500
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double db = P * L * L * L / (3.0 * p.EI), ds = P * L / p.shear_rigidity();
    const Eigen::Vector3d u = cantilever(L, 2, p, 0.0, 0.0, -P, 0.0);
    std::printf("  thin (d/L=1/500): tip_defl=%.6e  PL^3/3EI=%.6e (shear frac %.2e)\n",
                -u(1), db, ds / db);
    check(close(-u(1), db + ds, 1e-3), "thin plate: no shear locking (quartic)");
}

} // namespace

int main() {
    test_cantilever();
    test_axial_moment();
    test_thin_no_locking();
    if (g_failures == 0) {
        std::printf("OK: 5-node quartic plate element verified (tri15 edge)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
