// Plate (Mindlin/Timoshenko beam) element verification (P2.4). Soil-free cantilever: one
// end clamped, a transverse point load P at the other -> the tip deflection against the
// analytic (Timoshenko: bending + shear):
//   delta = P L^3/(3 EI) + P L/(k G A')     (k G A' = plate.shear_rigidity)
// Also pure axial (NL/EA) and rotation invariance at a random orientation. Selective
// reduced shear integration means no shear locking in thin plates -> exact with few
// elements. (See docs/references/structural-plate-formulation.md.)
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
using katai::core::plate::NodeCoords;

// Assemble + solve a straight cantilever of length L (ne 3-node elements) oriented at angle
// theta, clamped at node 0, with a tip load (fx_tip, fy_tip) and tip moment m_tip. Returns
// the tip nodal displacement (ux, uy, phi) in GLOBAL axes.
Eigen::Vector3d cantilever(double L, int ne, const PlateProps& p, double theta,
                           double fx_tip, double fy_tip, double m_tip) {
    const int nnodes = 2 * ne + 1;
    const int ndof = 3 * nnodes;
    const double c = std::cos(theta), s = std::sin(theta);
    const double h = L / (2 * ne);  // node spacing along the beam
    std::vector<double> X(nnodes), Y(nnodes);
    for (int i = 0; i < nnodes; ++i) { X[i] = c * (i * h); Y[i] = s * (i * h); }

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
    for (int e = 0; e < ne; ++e) {
        const int gn[3] = {2 * e, 2 * e + 2, 2 * e + 1};  // xi=-1, +1, 0
        NodeCoords Xe;
        for (int k = 0; k < 3; ++k) { Xe(k, 0) = X[gn[k]]; Xe(k, 1) = Y[gn[k]]; }
        const auto Ke = katai::core::plate::stiffness(Xe, p);
        for (int a = 0; a < 3; ++a)
            for (int ca = 0; ca < 3; ++ca)
                for (int b = 0; b < 3; ++b)
                    for (int cb = 0; cb < 3; ++cb)
                        K(3 * gn[a] + ca, 3 * gn[b] + cb) += Ke(3 * a + ca, 3 * b + cb);
    }

    Eigen::VectorXd f = Eigen::VectorXd::Zero(ndof);
    const int tip = nnodes - 1;
    f(3 * tip + 0) = fx_tip; f(3 * tip + 1) = fy_tip; f(3 * tip + 2) = m_tip;

    // Clamp node 0 (ux, uy, phi) via large-penalty -> reduced solve by row/col removal.
    std::vector<int> fixed{0, 1, 2};
    std::vector<int> free;
    std::vector<char> isfixed(ndof, 0);
    for (int d : fixed) isfixed[d] = 1;
    for (int d = 0; d < ndof; ++d) if (!isfixed[d]) free.push_back(d);
    const int nf = static_cast<int>(free.size());
    Eigen::MatrixXd Kff(nf, nf); Eigen::VectorXd ff(nf);
    for (int i = 0; i < nf; ++i) {
        ff(i) = f(free[i]);
        for (int j = 0; j < nf; ++j) Kff(i, j) = K(free[i], free[j]);
    }
    const Eigen::VectorXd uf = Kff.ldlt().solve(ff);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof);
    for (int i = 0; i < nf; ++i) u(free[i]) = uf(i);
    return {u(3 * tip + 0), u(3 * tip + 1), u(3 * tip + 2)};
}

void test_cantilever_bending_shear() {
    // Solid wall: E=3e7 kPa, thickness d=0.5 m -> EA=E d, EI=E d^3/12, nu=0.15.
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0, P = 100.0;  // tip transverse load [kN/m]

    const double db = P * L * L * L / (3.0 * p.EI);          // bending
    const double ds = P * L / p.shear_rigidity();            // shear (Timoshenko)
    const double delta = db + ds;

    // Horizontal cantilever, tip load in -y (transverse). Tip transverse defl = -uy.
    const Eigen::Vector3d u = cantilever(L, 8, p, 0.0, 0.0, -P, 0.0);
    const double tip_defl = -u(1);
    std::printf("  cantilever: tip_defl=%.6e  analytic=%.6e (bend=%.4e shear=%.4e)\n",
                tip_defl, delta, db, ds);
    check(close(tip_defl, delta, 1e-3), "tip deflection = PL^3/3EI + PL/kGA' (Timoshenko)");

    // Rotation invariance: same beam at 30deg, transverse load perpendicular to the axis.
    const double th = 30.0 * kPi / 180.0;
    const Eigen::Vector3d u2 =
        cantilever(L, 8, p, th, P * std::sin(th), -P * std::cos(th), 0.0);
    const double tip_defl2 = -u2(0) * std::sin(th) + u2(1) * std::cos(th);  // transverse comp
    check(close(std::fabs(tip_defl2), delta, 1e-3), "rotation invariance (30deg cantilever)");
}

void test_pure_axial_and_moment() {
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0;

    // Pure axial: tip axial load N -> elongation NL/EA.
    const double N = 500.0;
    const Eigen::Vector3d ua = cantilever(L, 8, p, 0.0, N, 0.0, 0.0);
    check(close(ua(0), N * L / p.EA, 1e-6), "pure axial elongation = NL/EA");

    // Tip moment M -> tip rotation phi = ML/EI, tip deflection = ML^2/(2EI).
    const double M = 200.0;
    const Eigen::Vector3d um = cantilever(L, 8, p, 0.0, 0.0, 0.0, M);
    check(close(um(2), M * L / p.EI, 1e-3), "tip rotation = ML/EI (pure bending)");
    check(close(um(1), M * L * L / (2.0 * p.EI), 1e-3), "tip deflection = ML^2/2EI");
}

void test_thin_plate_no_locking() {
    // Very thin plate (d/L = 1/500): shear deformation negligible; selective reduced
    // integration must NOT lock (tip defl ~ PL^3/3EI, not ~0).
    const double E = 3.0e7, d = 0.02, nu = 0.15, L = 10.0, P = 1.0;
    PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double db = P * L * L * L / (3.0 * p.EI), ds = P * L / p.shear_rigidity();
    const Eigen::Vector3d u = cantilever(L, 8, p, 0.0, 0.0, -P, 0.0);
    std::printf("  thin plate (d/L=1/500): tip_defl=%.6e  PL^3/3EI=%.6e (shear frac %.2e)\n",
                -u(1), db, ds / db);
    check(close(-u(1), db + ds, 2e-3), "thin plate: no shear locking (defl ~ PL^3/3EI)");
}

} // namespace

int main() {
    test_cantilever_bending_shear();
    test_pure_axial_and_moment();
    test_thin_plate_no_locking();
    if (g_failures == 0) {
        std::printf("OK: plate (Mindlin/Timoshenko) element verified (cantilever, axial, "
                    "moment, no shear locking)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
