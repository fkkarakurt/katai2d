// A plate in axial symmetry is a shell of revolution: the element against the closed forms of
// shell and plate theory, with nothing else in the model.
//
// The plane strain plate (plate.hpp) has a meridional strain, curvature and shear; turned about the
// axis it gains a hoop strain u_r / r and a hoop curvature cos(beta) phi / r, and it is
// integrated per radian like the soil. Three structures that have a closed form, each assembled
// from these elements alone (3- and 5-noded) and loaded by a pressure:
//   (a) a cylinder of radius R, free at its ends, under an internal pressure p: the hoop force is
//       N_t = p R (statics alone) and the radial displacement u_r = p R^2 / (EA (1 - nu^2)) (with
//       N_s = 0 the hoop strain carries EA (1 - nu^2));
//   (b) a circular plate of radius a clamped at its edge under a uniform pressure q: the centre
//       deflection of a Mindlin plate, w = q a^4 / (64 D) + q a^2 / (4 k G A'), D = EI (the
//       plate's EI is the plane strain one: it already carries 1 / (1 - nu^2));
//   (c) the same plate simply supported: w = q a^4 (5 + nu) / (64 D (1 + nu)) + q a^2 / (4 k G A').
// And the plane strain limit: far from the axis a strip of the shell is the plane strain plate.
//
// verify: KV-STR-029
//   oracle:   closed_form
//   source:   Timoshenko, S.P. & Woinowsky-Krieger, S. (1959). Theory of Plates and Shells, 2nd ed. McGraw-Hill, New York -- the cylinder under internal pressure (ch. 15) and the circular plate, clamped and simply supported (ch. 3); Mindlin, R.D. (1951). Influence of rotatory inertia and shear on flexural motions of isotropic, elastic plates. J. Appl. Mech. 18, 31-38 -- the shear term
//   locator:  EA = 1e6 kN/m, EI = 1e3 kNm/m, nu = 0.3, k = 5/6; (a) R = 2 m, height 4 m, 16 elements, p = 10 kPa; (b), (c) a = 3 m, 24 elements graded towards the edge, q = 10 kPa; each with the 3- and the 5-noded element
//   quantity: (a) N_t at mid-height [kN/m] and u_r [m]; (b), (c) the centre deflection [m]
//   expected: (a) N_t = p R = 20 kN/m, u_r = p R^2 / (EA (1 - nu^2)); (b) q a^4 / (64 EI) + q a^2 / (4 k G A'); (c) q a^4 (5 + nu) / (64 EI (1 + nu)) + q a^2 / (4 k G A')
//   band:     1e-6 relative on (a) (an exact state of the element); 2e-3 on (b), (c) (quadratic and quartic interpolation of a quartic deflection, measured)
#include <katai/fem/elements/plate.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = katai::core::plate;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

// A meridian from (r0, z0) to (r1, z1) in ne elements of `order` (3 or 5) nodes, with the node
// positions spaced by `grade` (1 = uniform). Assembled, loaded by a pressure p normal to the
// meridian (positive along the element's n = (-s, c)), and solved with the DOFs `fixed` held.
struct Shell {
    int nn = 0;
    std::vector<double> r, z;
    Eigen::MatrixXd K;
    Eigen::VectorXd f, u;
    pl::PlateProps P;
    int order = 3, ne = 0;
};

Shell build(int order, int ne, double r0, double z0, double r1, double z1, double grade, double p) {
    Shell S;
    S.order = order; S.ne = ne;
    const int per = order - 1;
    S.nn = ne * per + 1;
    S.P.EA = 1.0e6; S.P.EI = 1.0e3; S.P.nu = 0.3; S.P.axisym = true;
    // node positions along the meridian, t in [0,1], graded towards t = 1 when grade < 1
    for (int i = 0; i < S.nn; ++i) {
        const double t = (double)i / (S.nn - 1);
        const double g = grade == 1.0 ? t : 1.0 - std::pow(1.0 - t, 1.0 / grade);
        S.r.push_back(r0 + g * (r1 - r0));
        S.z.push_back(z0 + g * (z1 - z0));
    }
    const int nd = 3 * S.nn;
    S.K = Eigen::MatrixXd::Zero(nd, nd);
    S.f = Eigen::VectorXd::Zero(nd);
    for (int e = 0; e < ne; ++e) {
        // element node order: ends first, then interior (plate.hpp: [A, B, middle] / 5-node xi order)
        std::vector<int> nodes;
        if (order == 3) nodes = {e * 2, e * 2 + 2, e * 2 + 1};
        else nodes = {e * 4, e * 4 + 1, e * 4 + 2, e * 4 + 3, e * 4 + 4};
        Eigen::MatrixXd Ke;
        if (order == 3) {
            pl::NodeCoords X;
            for (int k = 0; k < 3; ++k) { X(k, 0) = S.r[nodes[k]]; X(k, 1) = S.z[nodes[k]]; }
            Ke = pl::stiffness(X, S.P);
        } else {
            pl::NodeCoords5 X;
            for (int k = 0; k < 5; ++k) { X(k, 0) = S.r[nodes[k]]; X(k, 1) = S.z[nodes[k]]; }
            Ke = pl::stiffness5(X, S.P);
        }
        for (int a = 0; a < 3 * order; ++a)
            for (int b = 0; b < 3 * order; ++b)
                S.K(3 * nodes[a / 3] + a % 3, 3 * nodes[b / 3] + b % 3) += Ke(a, b);
        // consistent pressure load, integral of N_i p r ds, 6-point Gauss along the element (exact enough)
        const double xa = S.r[nodes[0]], ya = S.z[nodes[0]], xb = S.r[nodes[order == 3 ? 1 : 4]], yb = S.z[nodes[order == 3 ? 1 : 4]];
        const double L = std::hypot(xb - xa, yb - ya), c = (xb - xa) / L, s = (yb - ya) / L;
        const double gx[6] = {-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142};
        const double gw[6] = {0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924};
        for (int q = 0; q < 6; ++q) {
            const double xi = gx[q];
            Eigen::VectorXd N(order);
            if (order == 3) {
                const Eigen::Vector3d n3 = pl::shape(xi);
                N << n3(0), n3(1), n3(2);
            } else {
                // Lagrange on -1, -0.5, 0, 0.5, 1
                const double xs[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
                for (int i = 0; i < 5; ++i) {
                    double v = 1.0;
                    for (int j = 0; j < 5; ++j) if (j != i) v *= (xi - xs[j]) / (xs[i] - xs[j]);
                    N(i) = v;
                }
            }
            double rq = 0.0;
            for (int k = 0; k < order; ++k) rq += N(k) * S.r[nodes[k]];
            const double wq = gw[q] * 0.5 * L * rq * p;
            for (int k = 0; k < order; ++k) {
                S.f(3 * nodes[k] + 0) += wq * N(k) * (-s);
                S.f(3 * nodes[k] + 1) += wq * N(k) * c;
            }
        }
    }
    return S;
}

void solve(Shell& S, const std::vector<int>& fixed) {
    const int nd = 3 * S.nn;
    std::vector<int> fr;
    for (int d = 0; d < nd; ++d)
        if (std::find(fixed.begin(), fixed.end(), d) == fixed.end()) fr.push_back(d);
    Eigen::MatrixXd Kf((int)fr.size(), (int)fr.size());
    Eigen::VectorXd ff((int)fr.size());
    for (size_t i = 0; i < fr.size(); ++i) {
        ff((int)i) = S.f(fr[i]);
        for (size_t j = 0; j < fr.size(); ++j) Kf((int)i, (int)j) = S.K(fr[i], fr[j]);
    }
    const Eigen::VectorXd uf = Kf.fullPivLu().solve(ff);
    S.u = Eigen::VectorXd::Zero(nd);
    for (size_t i = 0; i < fr.size(); ++i) S.u(fr[i]) = uf((int)i);
}

}  // namespace

int main() {
    std::printf("== A plate in axial symmetry: a shell of revolution (KV-STR-029) ==\n");
    const double kGA = pl::PlateProps{1.0e6, 1.0e3, 0.3}.shear_rigidity();
    for (int order : {3, 5}) {
        {
            // (a) the cylinder: meridian vertical at r = R, pressure outward (n = (-s, c) with
            // s = 1, c = 0 points to -r, so p < 0 pushes outward)
            const double R = 2.0, H = 4.0, p = 10.0;
            Shell S = build(order, 16 / (order - 1) * 2, R, 0.0, R, H, 1.0, -p);
            solve(S, {1});   // u_z at the bottom: the rigid translation along the axis
            const int mid = S.nn / 2;
            const double ur = S.u(3 * mid);
            // the hoop force from the element at mid-height
            pl::PlateForces F{};
            if (order == 3) {
                pl::NodeCoords X; pl::Dof u;
                const int e = (mid - 1) / 2, nodes[3] = {e * 2, e * 2 + 2, e * 2 + 1};
                for (int k = 0; k < 3; ++k) { X(k, 0) = S.r[nodes[k]]; X(k, 1) = S.z[nodes[k]]; for (int c = 0; c < 3; ++c) u(3 * k + c) = S.u(3 * nodes[k] + c); }
                F = pl::forces(X, S.P, u, 0.0);
            } else {
                pl::NodeCoords5 X; pl::Dof5 u;
                const int e = (mid - 1) / 4, n0 = e * 4;
                for (int k = 0; k < 5; ++k) { X(k, 0) = S.r[n0 + k]; X(k, 1) = S.z[n0 + k]; for (int c = 0; c < 3; ++c) u(3 * k + c) = S.u(3 * (n0 + k) + c); }
                F = pl::forces5(X, S.P, u, 0.0);
            }
            const double want_u = p * R * R / (S.P.EA * (1.0 - 0.09));
            std::printf("      %d-noded cylinder: N_t %.6f kN/m (p R %.3f), u_r %.6e m (%.6e), N_s %.2e\n",
                        order, F.Nt, p * R, ur, want_u, F.N);
            char lbl[160];
            std::snprintf(lbl, sizeof lbl, "(a) %d-noded: a cylinder under internal pressure carries N_t = p R and expands p R^2 / (EA (1 - nu^2))", order);
            check(rel(F.Nt, p * R) < 1e-6 && rel(ur, want_u) < 1e-6 && std::fabs(F.N) < 1e-6 * p * R, lbl);
        }
        for (int clamped : {1, 0}) {
            // (b), (c) the circular plate: meridian along r from the axis to a, pressure down
            // (n = (-s, c) = (0, 1) is up, so p < 0 loads downward)
            const double a = 3.0, q = 10.0;
            Shell S = build(order, 24 / (order - 1) * 2, 0.0, 0.0, a, 0.0, 0.7, -q);
            std::vector<int> fixed = {0, 2};     // the axis: u_r = 0 and the slope phi = 0
            fixed.push_back(3 * (S.nn - 1) + 1); // the edge: w = 0
            if (clamped) fixed.push_back(3 * (S.nn - 1) + 2);   // and phi = 0
            // u_r is free elsewhere (the membrane is unloaded)
            solve(S, fixed);
            const double w0 = S.u(1);
            const double shear = q * a * a / (4.0 * kGA);
            const double want = clamped ? -(q * std::pow(a, 4) / (64.0 * S.P.EI) + shear)
                                        : -(q * std::pow(a, 4) * (5.0 + 0.3) / (64.0 * S.P.EI * 1.3) + shear);
            std::printf("      %d-noded circular plate, %s: centre w %.6e m, Mindlin %.6e m (%+.4f%%)\n",
                        order, clamped ? "clamped" : "simply supported", w0, want, 100.0 * (w0 - want) / want);
            char lbl[160];
            std::snprintf(lbl, sizeof lbl, "(%s) %d-noded: a %s circular plate deflects as the Mindlin plate", clamped ? "b" : "c", order,
                          clamped ? "clamped" : "simply supported");
            check(rel(w0, want) < 2e-3, lbl);
        }
    }
    {
        // The plane strain limit: a short straight piece of a shell far from the axis, bent by
        // its end moment, is the plane strain plate (the hoop terms are of order 1/r).
        pl::PlateProps P{1.0e6, 1.0e3, 0.3};
        pl::PlateProps A = P; A.axisym = true;
        pl::NodeCoords X;
        const double r0 = 1.0e6;
        X << r0, 0.0, r0 + 1.0, 0.0, r0 + 0.5, 0.0;
        const Eigen::MatrixXd Kp = pl::stiffness(X, P), Ka = pl::stiffness(X, A) / r0;
        std::printf("      far from the axis: |K_axisym / r - K_plane| / |K_plane| = %.2e\n", (Ka - Kp).norm() / Kp.norm());
        check((Ka - Kp).norm() / Kp.norm() < 1e-5, "far from the axis a strip of the shell is the plane strain plate (per radian / r)");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
