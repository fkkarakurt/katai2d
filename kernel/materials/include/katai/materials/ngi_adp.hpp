#pragma once
// NGI-ADP: an anisotropic undrained shear strength model for clay, in total stress.
//
// Grimstad, G., Andresen, L. & Jostad, H.P. (2012). NGI-ADP: Anisotropic shear strength model for
// clay. International Journal for Numerical and Analytical Methods in Geomechanics 36(4), 483-497;
// Andresen, L. & Jostad, H.P. (1999), Andresen (2002); Davis, E.H. & Christian, J.T. (1971).
// Bearing capacity of anisotropic cohesive soil. J. Soil Mech. Found. Div. ASCE 97(SM5), 753-769.
// Formulation and the corrections made to its printed form: docs/references/ngi-adp-formulation.md.
//
// The model, in its three-dimensional form (compression positive, y the vertical -- depositional
// -- direction, stresses [xx, yy, zz, xy], engineering shear strain):
//
//   modified deviator   s^_n  = dev( s_n - (1 - k) s0_n + k a_n ),  a = (Da/3, -2 Da/3, Da/3)
//                       s^_xy = c (t_xy - (1 - k) t0_xy),           c = (suA + suP) / (2 suDSS)
//                       Da = suA - suP
//   yield               f = sqrt(H(w) J2^) - k (suA + suP) / 2
//                       H(w) = cos^2( acos(1 - 2 a1 w) / 6 ),  w = 27/4 J3^^2 / J2^^3
//   plastic potential   g = sqrt(J2^)                      (no H: von Mises in the modified space)
//   hardening           k = 2 sqrt(gp / gpf) / (1 + gp / gpf),  k = 1 from gp >= gpf
//                       dgp = sqrt(H) sqrt( 2/3 sum (de_i - de_j)^2 + dg_xy^2 )   (plastic parts)
//   failure strain      gpf(th) elliptical between the plastic failure strains in triaxial
//                       compression (cos 2th = 1), simple shear (0) and triaxial extension (-1),
//                       cos 2th = sqrt(3)/2 s^_yy / sqrt(J2^)
//   elasticity          isotropic, G (the unloading-reloading modulus) and the undrained nu
//
// The strengths are those of plane strain: active suA (vertical compression), passive suP (vertical
// extension) and direct simple shear suDSS. a1 is fixed by the triaxial-compression strength
// suC,TX = r_tx suA (r_tx = 0.99): exact Tresca, a1 = 1, makes it equal to suA.
//
// THE FAILURE STRAINS ENTERED ARE TOTAL STRAINS, as a laboratory test reports them: the strain at
// which the test, run from the in-situ state, reaches its strength. The hardening reads plastic
// strains, so each is less the elastic part of its own test, the stress it climbs over the shear
// modulus: gpf_C = gf_C - (suC,TX - t0)/G, gpf_E = gf_E - (suE,TX + t0)/G, gpf_DSS = gf_DSS - suDSS/G,
// t0 = (s0_yy - s0_xx)/2 the point's initial shear stress. A test run from that state reaches its
// strength at the strain entered (KV-CST-024).
//
// THE INITIAL STATE IS THE POINT'S OWN: s0 is the stress the point carries when the model first
// acts on it (the K0 procedure's, or the phase before's), and at k = 0 the yield surface is that
// point. The in-situ shear stress follows from it rather than from a separate input.

#include <Eigen/Core>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>

namespace katai::core::adp {

struct Params {
    double G = 0.0;        // shear modulus [kPa]
    double nu = 0.495;     // undrained Poisson's ratio
    double suA = 0.0;      // plane strain active strength [kPa]
    double suP = 0.0;      // plane strain passive strength [kPa]
    double suDSS = 0.0;    // direct simple shear strength [kPa]
    double gfC = 0.0;      // total shear strain at failure, triaxial compression (3/2 eps1) [-]
    double gfE = 0.0;      // total shear strain at failure, triaxial extension [-]
    double gfDSS = 0.0;    // total shear strain at failure, direct simple shear [-]
    double r_tx = 0.99;    // suC,TX / suA (fixes the Tresca rounding a1)
};

using Vec4 = Eigen::Matrix<double, 4, 1>;   // [xx, yy, zz, xy]
using Mat4 = Eigen::Matrix<double, 4, 4>;

struct Derived {
    double S = 0.0, Da = 0.0, c = 1.0, a1 = 1.0, H1 = 0.75;
    double suC_tx = 0.0, suE_tx = 0.0;
};

inline Derived derive(const Params& P) {
    Derived d;
    d.S = P.suA + P.suP;
    d.Da = P.suA - P.suP;
    d.c = d.S / (2.0 * P.suDSS);
    // Triaxial: w = 1, J2^ = q^^2/3, so sqrt(H1/3) |q^| = S/2 and suC,TX = (Da + S x)/2 with
    // x = sqrt(3 / (4 H1)). The rounding is fixed by suC,TX = r_tx suA.
    const double x = (2.0 * P.r_tx * P.suA - d.Da) / d.S;
    d.H1 = std::clamp(0.75 / (x * x), 0.75, 1.0);
    const double ang = std::acos(std::clamp(std::sqrt(d.H1), -1.0, 1.0));
    d.a1 = std::clamp(0.5 * (1.0 - std::cos(6.0 * ang)), 0.0, 1.0);
    const double xs = std::sqrt(0.75 / d.H1);
    d.suC_tx = 0.5 * (d.Da + d.S * xs);
    d.suE_tx = 0.5 * (d.S * xs - d.Da);
    return d;
}

inline Mat4 elastic(const Params& P) {
    const double K = 2.0 * P.G * (1.0 + P.nu) / (3.0 * (1.0 - 2.0 * P.nu));
    const double lam = K - 2.0 * P.G / 3.0;
    Mat4 D = Mat4::Zero();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) D(i, j) = lam + (i == j ? 2.0 * P.G : 0.0);
    D(3, 3) = P.G;
    return D;
}

// Modified deviator s^ and its invariants.
struct Mod {
    Vec4 s = Vec4::Zero();
    double J2 = 0.0, J3 = 0.0, w = 0.0, H = 1.0;
};

inline Mod modified(const Vec4& sig, const Vec4& s0, double k, const Derived& d) {
    Mod m;
    Vec4 t;
    t(0) = sig(0) - (1.0 - k) * s0(0) + k * d.Da / 3.0;
    t(1) = sig(1) - (1.0 - k) * s0(1) - k * 2.0 * d.Da / 3.0;
    t(2) = sig(2) - (1.0 - k) * s0(2) + k * d.Da / 3.0;
    const double ph = (t(0) + t(1) + t(2)) / 3.0;
    m.s(0) = t(0) - ph;
    m.s(1) = t(1) - ph;
    m.s(2) = t(2) - ph;
    m.s(3) = d.c * (sig(3) - (1.0 - k) * s0(3));
    m.J2 = 0.5 * (m.s(0) * m.s(0) + m.s(1) * m.s(1) + m.s(2) * m.s(2)) + m.s(3) * m.s(3);
    m.J3 = m.s(0) * m.s(1) * m.s(2) - m.s(2) * m.s(3) * m.s(3);
    m.w = m.J2 > 1e-300 ? std::clamp(6.75 * m.J3 * m.J3 / (m.J2 * m.J2 * m.J2), 0.0, 1.0) : 0.0;
    m.H = std::pow(std::cos(std::acos(std::clamp(1.0 - 2.0 * d.a1 * m.w, -1.0, 1.0)) / 6.0), 2);
    return m;
}

// The plastic failure strain at the modified deviator's orientation (eq. 9-10, with cos^2 2th
// under the root: the form that returns gpf_C, gpf_DSS and gpf_E at cos 2th = 1, 0, -1).
inline double failure_strain(const Mod& m, double gC, double gE, double gD) {
    const double c2 = m.J2 > 1e-300 ? std::clamp(0.8660254037844386 * m.s(1) / std::sqrt(m.J2), -1.0, 1.0) : 1.0;
    const double RA = 0.5 * (gE - gC), RB = 0.5 * (gE + gC), RC = std::sqrt(gE * gC), RD = gD * RB / RC;
    const double num = RB * RD * std::sqrt((RD * RD - RC * RC) * c2 * c2 + RC * RC) - RD * RD * RA * c2;
    const double den = RB * RB - (RB * RB - RD * RD) * c2 * c2;
    return num / den;
}

inline double kappa_of(double gp, double gpf) {
    if (gp <= 0.0) return 0.0;
    if (gp >= gpf) return 1.0;
    const double r = gp / gpf;
    return 2.0 * std::sqrt(r) / (1.0 + r);
}

// The flow direction dg/dsigma (strain-like, engineering shear) and the plastic shear strain per
// unit multiplier, sqrt(H) sqrt(2/3 sum (n_i - n_j)^2 + n_xy^2).
inline void flow(const Mod& m, const Derived& d, Vec4& n, double& per_lambda) {
    const double r = std::sqrt(std::max(m.J2, 1e-300));
    n(0) = m.s(0) / (2.0 * r);
    n(1) = m.s(1) / (2.0 * r);
    n(2) = m.s(2) / (2.0 * r);
    n(3) = d.c * m.s(3) / r;
    const double a = n(0) - n(1), b = n(0) - n(2), e = n(1) - n(2);
    per_lambda = std::sqrt(m.H) * std::sqrt((2.0 / 3.0) * (a * a + b * b + e * e) + n(3) * n(3));
}

struct State {
    Vec4 sig = Vec4::Zero();   // compression positive
    Vec4 s0 = Vec4::Zero();    // the initial state the model translates from
    double gp = 0.0;           // plastic shear strain
};

// Plastic failure strains of the point, in the hardening's measure, from the total ones entered
// (header note). Kept above
// a small floor: an input whose elastic part reaches its total has no plastic room, and the
// input contract refuses it; the floor only keeps the arithmetic finite.
inline void plastic_failure_strains(const Params& P, const Derived& d, const Vec4& s0,
                                    double& gC, double& gE, double& gD) {
    const double t0 = 0.5 * (s0(1) - s0(0));
    const double floor = 1e-7;
    // In a triaxial test the hardening's measure is sqrt(4 H(1)/3) (eps1 - eps3) of the plastic
    // strain -- 1 at exact Tresca, 1.4% more at the model's rounding: the failure strain the test
    // is to reach is read in that measure.
    const double tx = std::sqrt(4.0 * d.H1 / 3.0);
    gC = std::max(tx * (P.gfC - (d.suC_tx - t0) / P.G), floor);
    gE = std::max(tx * (P.gfE - (d.suE_tx + t0) / P.G), floor);
    gD = std::max(P.gfDSS - P.suDSS / P.G, floor);
}

// One strain increment from a committed state: backward Euler on [sigma, dlambda, gp] (a Newton
// iteration with a finite-difference Jacobian; the hardening's square root at gp = 0 makes an
// analytical one awkward and buys nothing at six unknowns). Returns false if it did not converge.
// Given `T`, it also returns the consistent tangent of the increment: differentiating the residual at
// the answer with respect to the strain (only the trial stress depends on it) gives
// J dX = (D / S) d(eps) in its first four rows, so the tangent is one more Jacobian and a solve,
// where differencing the integration itself is four more integrations of the whole Newton iteration.
inline bool return_map(const Params& P, const Derived& d, const State& c, const Vec4& de, State& out,
                       Mat4* T = nullptr) {
    const Mat4 D = elastic(P);
    double gC, gE, gD;
    plastic_failure_strains(P, d, c.s0, gC, gE, gD);
    out = c;
    const Vec4 tr = c.sig + D * de;
    // k at a stress and a plastic strain: the failure strain is read at the orientation of the
    // modified deviator, which itself moves with k, so the two are solved together (the map's slope is
    // below 0.2, and three passes reach what more would). One pass read the orientation at a k that was not the answer:
    // measured on direct simple shear against the published curves, 5.7% and 4.1% low before the peak in
    // two clays (rms 2.3% and 1.5%) where triaxial compression and extension agreed to 0.1%; the peak,
    // where k = 1, does not depend on it.
    const auto kappa_at = [&](const Vec4& sig, double gp) {
        double k = kappa_of(gp, gC);
        for (int pass = 0; pass < 3; ++pass) {   // measured: three reach the curves, more do not move them
            const double next = kappa_of(gp, failure_strain(modified(sig, c.s0, k, d), gC, gE, gD));
            const bool settled = std::fabs(next - k) <= 1e-12;   // at the peak, or where the three strains agree
            k = next;
            if (settled) break;
        }
        return k;
    };
    {
        const double k0 = kappa_at(tr, c.gp);
        const Mod mk = modified(tr, c.s0, k0, d);
        if (std::sqrt(mk.H * mk.J2) - k0 * 0.5 * d.S <= 1e-10 * d.S) {   // elastic
            out.sig = tr;
            if (T) *T = D;
            return true;
        }
    }
    const double ss = d.S, gs = std::max(gD, 1e-6);
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;
    const auto residual = [&](const Vec6& X, Vec6& R) {
        const Vec4 sig = X.head<4>();
        const double dl = X(4), gp = X(5);
        const double k = kappa_at(sig, gp);
        const Mod m = modified(sig, c.s0, k, d);
        Vec4 n;
        double per;
        flow(m, d, n, per);
        R.head<4>() = (sig - tr + dl * (D * n)) / ss;
        R(4) = (std::sqrt(m.H * m.J2) - k * 0.5 * d.S) / ss;
        R(5) = (gp - c.gp - dl * per) / gs;
    };
    Vec6 X;
    X.head<4>() = tr;
    X(4) = 0.0;
    X(5) = c.gp;
    // A start off the singular corner: a small plastic step along the trial's direction.
    {
        const Mod m = modified(tr, c.s0, kappa_of(c.gp, gC), d);
        Vec4 n;
        double per;
        flow(m, d, n, per);
        const double f = std::sqrt(m.H * m.J2);
        const double dl0 = 0.5 * f / (2.0 * P.G);
        X(4) = dl0;
        X(5) = c.gp + dl0 * per;
        X.head<4>() = tr - dl0 * (D * n);
    }
    Vec6 R;
    const auto jacobian = [&]() {
        Mat6 J;
        for (int j = 0; j < 6; ++j) {
            Vec6 Xp = X;
            const double h = j < 4 ? 1e-7 * ss : 1e-7 * gs + 1e-12;
            Xp(j) += h;
            Vec6 Rp;
            residual(Xp, Rp);
            J.col(j) = (Rp - R) / h;
        }
        return J;
    };
    for (int it = 0; it < 60; ++it) {
        residual(X, R);
        if (R.norm() < 1e-12) {
            out.sig = X.head<4>();
            out.gp = X(5);
            if (X(4) < -1e-14) return false;
            if (T) {
                Eigen::Matrix<double, 6, 4> B = Eigen::Matrix<double, 6, 4>::Zero();
                B.topRows<4>() = D / ss;
                const Eigen::Matrix<double, 6, 4> dXde = jacobian().partialPivLu().solve(B);
                *T = dXde.topRows<4>();
                if (!T->allFinite()) *T = D;
            }
            return true;
        }
        const Mat6 J = jacobian();
        Vec6 dX = J.partialPivLu().solve(-R);
        if (!dX.allFinite()) return false;
        // Keep the multiplier and the plastic strain where they mean something.
        double a = 1.0;
        for (int halve = 0; halve < 30; ++halve) {
            Vec6 Xn = X + a * dX;
            if (Xn(4) < 0.0) Xn(4) = 0.5 * X(4);
            if (Xn(5) < c.gp) Xn(5) = 0.5 * (X(5) + c.gp);
            Vec6 Rn;
            residual(Xn, Rn);
            if (Rn.norm() < R.norm() || halve == 29) { X = Xn; break; }
            a *= 0.5;
        }
    }
    return false;
}

// The increment, cut into substeps when one step does not converge.
// A STRAIN TOO SMALL TO RETURN. The hardening has an infinite slope where it starts, kappa ~ 2 sqrt(gp / gf), and from the
// stress a point started at (s = s0) the flow direction is singular too, so the return of a very small increment may not
// converge, and cutting the increment only makes that worse: every level of the cut then costs its whole iteration. Measured on
// a shear column shaken in a nonlinear Dynamic phase (99 nodes): no result in 280 s, each stress point spending 15 ms of cuts
// on every evaluation. An increment below 1e-3 of the strength, as a stress, is taken elastically: that is the size of the
// error, the next larger increment returns to the surface, and a test from rest loses nothing it can read.
inline bool too_small_to_return(const Params& P, const Derived& d, const Vec4& de) {
    return 2.0 * P.G * de.lpNorm<Eigen::Infinity>() <= 1e-3 * d.S;
}

inline bool integrate(const Params& P, const State& c, const Vec4& de, State& out, int depth = 0) {
    const Derived d = derive(P);
    if (return_map(P, d, c, de, out)) return true;
    if (depth <= 1 && too_small_to_return(P, d, de)) {
        out = c;
        out.sig = c.sig + elastic(P) * de;
        return true;
    }
    if (depth > 10) return false;
    State mid;
    if (!integrate(P, c, 0.5 * de, mid, depth + 1)) return false;
    return integrate(P, mid, 0.5 * de, out, depth + 1);
}

// The consistent tangent, by differencing the integration itself: for an increment that had to be cut.
inline Mat4 tangent(const Params& P, const State& c, const Vec4& de, const State& out) {
    Mat4 T;
    const double h = 1e-8 + 1e-6 * de.lpNorm<Eigen::Infinity>();
    for (int j = 0; j < 4; ++j) {
        Vec4 dp = de;
        dp(j) += h;
        State o;
        if (!integrate(P, c, dp, o)) return elastic(P);
        T.col(j) = (o.sig - out.sig) / h;
    }
    return T;
}

// The increment and its consistent tangent.
inline bool integrate_tangent(const Params& P, const State& c, const Vec4& de, State& out, Mat4& T) {
    const Derived d = derive(P);
    if (return_map(P, d, c, de, out, &T)) return true;
    if (!integrate(P, c, de, out, 1)) return false;
    T = too_small_to_return(P, d, de) ? elastic(P) : tangent(P, c, de, out);
    return true;
}

}  // namespace katai::core::adp
