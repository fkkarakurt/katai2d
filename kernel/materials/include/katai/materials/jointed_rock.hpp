#pragma once
// JOINTED ROCK: a rock mass whose strength is that of its joints.
//
// The intact rock is elastic -- isotropic, or transversely isotropic about the normal to the first
// joint set, the stratification of a layered rock -- and it fails only by slip or opening on up to
// three sets of joints, each a family of parallel planes (smeared, "ubiquitous", over the whole
// material). In plane strain and in axial symmetry the planes contain the out-of-plane direction
// (the strike is normal to the model), so a plane is described by its dip alpha in the model plane,
// measured from the x axis, anticlockwise positive.
//
// On a plane with unit normal n and unit tangent s = (cos alpha, sin alpha), n = (-sin alpha,
// cos alpha), the traction is
//     sigma_n = n . sigma n        (tension positive)
//     tau     = s . sigma n
// and both are LINEAR in the stress. The strength of a plane is Coulomb's with a tension cut-off,
//     f_s+ =  tau + sigma_n tan(phi) - c <= 0
//     f_s- = -tau + sigma_n tan(phi) - c <= 0
//     f_t  =  sigma_n - sigma_t           <= 0      (sigma_t <= c / tan(phi), the apex)
// with the plastic potential g = |tau| + sigma_n tan(psi) for slip and sigma_n for opening.
//
// THE RETURN IS EXACT. Every surface is a plane in stress space and the elasticity is constant, so
// the backward Euler step needs no iteration: for a set J of active surfaces the multipliers solve
// the linear system  A dl = f(sigma_trial),  A_kj = a_k . D b_j, and the stress is
// sigma_trial - D B dl. Which set is active is found by enumerating every set of up to three
// surfaces (the in-plane stress has three components, so no more than three can be independent)
// and keeping the ones that satisfy the Kuhn-Tucker conditions (every dl >= 0, every surface <= 0
// afterwards). With associated flow exactly one does; with non-associated flow more than one can,
// and the one closest to the trial stress in the energy norm is taken. The consistent tangent of
// that set is closed form: D - D B A^-1 A^T D (non-symmetric when psi != phi).
//
// Stress and strain vectors are [xx, yy, zz, xy], tension positive, engineering shear strain; zz
// is the out-of-plane (plane strain) or hoop (axial symmetry) direction.
//
// Elasticity, in the axes of the first plane (s along it, n its normal, z out of plane; s and z
// span the plane of isotropy):
//     eps_ss = sig_ss / E1 - nu1 sig_zz / E1 - nu2 sig_nn / E2
//     eps_zz = sig_zz / E1 - nu1 sig_ss / E1 - nu2 sig_nn / E2
//     eps_nn = sig_nn / E2 - nu2 (sig_ss + sig_zz) / E2
//     gam_sn = tau_sn / G2
// (Lekhnitskii, S.G. (1963). Theory of Elasticity of an Anisotropic Elastic Body. Holden-Day, San
// Francisco), positive definite when E1, E2, G2 > 0, |nu1| < 1 and 1 - nu1 - 2 nu2^2 E1/E2 > 0.
// The strength model is the ubiquitous joint, or multilaminate, idealisation of Zienkiewicz, O.C.
// & Pande, G.N. (1977). Time-dependent multilaminate model of rocks -- a numerical study of
// deformation and failure of rock masses. Int. J. Numer. Anal. Meth. Geomech. 1(3), 219-247; the
// strength of a rock with one plane of weakness is Jaeger's (Jaeger, J.C. (1960). Shear failure of
// anisotropic rocks. Geological Magazine 97(1), 65-72).

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace katai::core::jr {

using Vec4 = Eigen::Vector4d;
using Mat4 = Eigen::Matrix4d;

struct Plane {
    double alpha = 0.0;    // dip from the x axis [rad], anticlockwise positive
    double c = 0.0;        // cohesion [kPa]
    double phi = 0.0;      // friction angle [rad]
    double psi = 0.0;      // dilatancy angle [rad]
    bool tcut = true;      // tension cut-off on this plane
    double sigt = 0.0;     // tensile strength [kPa] (>= 0, held at or below the apex c / tan phi)
};

struct Params {
    double E1 = 0.0, nu1 = 0.0;       // the plane of isotropy (the whole rock when isotropic)
    bool aniso = false;               // transversely isotropic about the normal to plane 0
    double E2 = 0.0, nu2 = 0.0, G2 = 0.0;
    int n = 1;                        // number of joint sets, 1..3
    Plane p[3];
};

// Whether the elastic constants describe a stable (positive definite) material; the reason if not.
inline const char* elastic_admissible(const Params& P) {
    if (!(P.E1 > 0.0)) return "E1 must be positive";
    if (!P.aniso) return (P.nu1 > -1.0 && P.nu1 < 0.5) ? nullptr : "nu must lie between -1 and 0.5";
    if (!(P.E2 > 0.0)) return "E2 must be positive";
    if (!(P.G2 > 0.0)) return "G2 must be positive";
    if (!(P.nu1 > -1.0 && P.nu1 < 1.0)) return "nu1 must lie between -1 and 1";
    if (!(1.0 - P.nu1 - 2.0 * P.nu2 * P.nu2 * P.E1 / P.E2 > 0.0))
        return "1 - nu1 - 2 nu2^2 E1/E2 must be positive (the elastic energy would not be)";
    return nullptr;
}

// The rotation taking global engineering strains [xx, yy, zz, xy] to the axes of a plane of dip
// alpha [ss, nn, zz, sn]; stresses transform with its transpose the other way (sigma_g = R^T sigma_l).
inline Mat4 strain_rotation(double alpha) {
    const double c = std::cos(alpha), s = std::sin(alpha);
    Mat4 R = Mat4::Zero();
    R(0, 0) = c * c;      R(0, 1) = s * s;      R(0, 3) = c * s;
    R(1, 0) = s * s;      R(1, 1) = c * c;      R(1, 3) = -c * s;
    R(2, 2) = 1.0;
    R(3, 0) = -2.0 * c * s; R(3, 1) = 2.0 * c * s; R(3, 3) = c * c - s * s;
    return R;
}

// The elastic operator, 4 x 4 on [xx, yy, zz, xy].
inline Mat4 elastic(const Params& P) {
    if (!P.aniso) {
        const double e = P.E1, v = P.nu1;
        const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
        Mat4 D = Mat4::Zero();
        D(0, 0) = D(1, 1) = D(2, 2) = f * (1.0 - v);
        D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(1, 2) = D(2, 1) = f * v;
        D(3, 3) = f * (1.0 - 2.0 * v) / 2.0;
        return D;
    }
    Mat4 C = Mat4::Zero();   // local [ss, nn, zz, sn]
    C(0, 0) = 1.0 / P.E1;      C(0, 1) = -P.nu2 / P.E2;   C(0, 2) = -P.nu1 / P.E1;
    C(1, 0) = -P.nu2 / P.E2;   C(1, 1) = 1.0 / P.E2;      C(1, 2) = -P.nu2 / P.E2;
    C(2, 0) = -P.nu1 / P.E1;   C(2, 1) = -P.nu2 / P.E2;   C(2, 2) = 1.0 / P.E1;
    C(3, 3) = 1.0 / P.G2;
    const Mat4 Dl = C.inverse();
    const Mat4 R = strain_rotation(P.p[0].alpha);
    return R.transpose() * Dl * R;
}

// The gradients of sigma_n and tau of a plane, as stress-vector duals ([xx, yy, zz, xy], so that
// a . sigma is the traction and a is also the plastic strain direction with engineering shear).
inline void traction_gradients(double alpha, Vec4& an, Vec4& at) {
    const double c = std::cos(alpha), s = std::sin(alpha);
    an << s * s, c * c, 0.0, -2.0 * c * s;     // n = (-s, c)
    at << -c * s, c * s, 0.0, c * c - s * s;   // s = (c, s)
}

inline double sigma_n(const Vec4& sig, double alpha) {
    Vec4 an, at;
    traction_gradients(alpha, an, at);
    return an.dot(sig);
}
inline double tau(const Vec4& sig, double alpha) {
    Vec4 an, at;
    traction_gradients(alpha, an, at);
    return at.dot(sig);
}

struct Surface {
    Vec4 a, b;      // f = a . sigma - r; plastic flow along b
    double r = 0.0;
    int plane = 0;
};

inline int surfaces(const Params& P, Surface out[9]) {
    int k = 0;
    for (int i = 0; i < std::clamp(P.n, 1, 3); ++i) {
        const Plane& q = P.p[i];
        Vec4 an, at;
        traction_gradients(q.alpha, an, at);
        const double tf = std::tan(q.phi), tp = std::tan(q.psi);
        out[k++] = Surface{at + tf * an, at + tp * an, q.c, i};
        out[k++] = Surface{-at + tf * an, -at + tp * an, q.c, i};
        // THE APEX OPENS THE JOINT, CUT-OFF OR NOT. Without a tension cut-off a joint carries
        // tension up to the apex of its Coulomb line, sigma_n = c / tan(phi); a trial beyond it
        // must come back to the apex, and that needs a normal plastic strain. Slip with psi < phi
        // has too little of one (with psi = 0 none at all: the two slip flows are +-a_t, and no
        // combination of them moves sigma_n), so without this surface such a trial had no
        // admissible active set and the point kept its elastic trial, far outside its strength --
        // measured on a slope on joints with c = 20 kPa, phi = 30, psi = 0 and no cut-off: a
        // "mechanism" at 43% of its self-weight. The joint opens at the apex instead.
        if (q.tcut) {
            double st = std::max(0.0, q.sigt);
            if (tf > 0.0) st = std::min(st, q.c / tf);
            out[k++] = Surface{an, an, st, i};
        } else if (tf > 0.0) {
            out[k++] = Surface{an, an, q.c / tf, i};
        }
    }
    return k;
}

struct Return {
    Vec4 sig = Vec4::Zero();
    Mat4 T = Mat4::Zero();      // the consistent tangent (D where elastic)
    bool plastic = false;
    bool ok = true;             // false: no admissible active set was found (the stress is the trial)
    int active[3] = {-1, -1, -1};
    int n_active = 0;
};

// The stress from a trial (elastic predictor) stress, with the operator D the predictor used.
inline Return return_map(const Params& P, const Mat4& D, const Vec4& trial) {
    Return R;
    R.sig = trial;
    R.T = D;
    Surface S[9];
    const int ns = surfaces(P, S);
    double scale = 1.0 + trial.cwiseAbs().maxCoeff();
    for (int k = 0; k < ns; ++k) scale = std::max(scale, std::fabs(S[k].r));
    const double tol = 1e-10 * scale;
    double ftr[9];
    bool any = false;
    for (int k = 0; k < ns; ++k) {
        ftr[k] = S[k].a.dot(trial) - S[k].r;
        any = any || ftr[k] > tol;
    }
    if (!any) return R;

    double best = -1.0;
    Vec4 best_sig = trial;
    int best_set[3] = {-1, -1, -1}, best_n = 0;
    Eigen::Matrix3d best_Ainv = Eigen::Matrix3d::Zero();
    const auto try_set = [&](const int* J, int m) {
        Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
        Eigen::Vector3d f = Eigen::Vector3d::Zero();
        for (int i = 0; i < m; ++i) {
            f(i) = ftr[J[i]];
            const Vec4 Db = D * S[J[i]].b;
            for (int j = 0; j < m; ++j) A(j, i) = S[J[j]].a.dot(Db);
        }
        const Eigen::MatrixXd Am = A.topLeftCorner(m, m);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(Am);
        if (lu.rank() < m) return;
        const Eigen::VectorXd dl = lu.solve(Eigen::VectorXd(f.head(m)));
        for (int i = 0; i < m; ++i)
            if (dl(i) < -1e-12 * (1.0 + dl.cwiseAbs().maxCoeff())) return;
        Vec4 bsum = Vec4::Zero();
        for (int i = 0; i < m; ++i) bsum += dl(i) * S[J[i]].b;
        const Vec4 dsig = D * bsum;
        const Vec4 sig = trial - dsig;
        for (int k = 0; k < ns; ++k)
            if (S[k].a.dot(sig) - S[k].r > tol) return;
        const double energy = bsum.dot(dsig);   // (trial - sigma) . C (trial - sigma)
        if (best < 0.0 || energy < best) {
            best = energy;
            best_sig = sig;
            best_n = m;
            for (int i = 0; i < 3; ++i) best_set[i] = i < m ? J[i] : -1;
            best_Ainv.setZero();
            best_Ainv.topLeftCorner(m, m) = Am.inverse();
        }
    };
    int J[3];
    // FIRST THE ACTIVE SET THE TRIAL ITSELF NAMES (Simo & Hughes 1998, Sec. 5.2): every violated
    // surface, a multiplier that comes out negative dropped, a surface violated afterwards added.
    // Where it ends in an admissible state that is the answer. This matters only for non-associated
    // flow, where more than one set can satisfy the Kuhn-Tucker conditions: measured on two joint
    // sets symmetric about the loading, slip on either set alone is admissible too and is the
    // closer in the energy norm, but which one a closest-set rule picked depended on round-off,
    // the symmetric answer was lost and the tangent jumped between the two. The enumeration below
    // remains for what this does not resolve (more than three surfaces violated, or a cycle).
    {
        int cur[9], m = 0;
        for (int k = 0; k < ns; ++k)
            if (ftr[k] > tol) cur[m++] = k;
        for (int guard = 0; guard < 12 && m > 0 && m <= 3; ++guard) {
            const double before = best;
            try_set(cur, m);
            if (best >= 0.0 && before < 0.0) break;   // admissible: done
            // solve again to see which multiplier is negative or which surface is violated
            Eigen::MatrixXd A(m, m);
            Eigen::VectorXd f(m);
            for (int i = 0; i < m; ++i) {
                f(i) = ftr[cur[i]];
                for (int j = 0; j < m; ++j) A(j, i) = S[cur[j]].a.dot(D * S[cur[i]].b);
            }
            Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
            if (lu.rank() < m) break;
            const Eigen::VectorXd dl = lu.solve(f);
            int worst = -1;
            for (int i = 0; i < m; ++i)
                if (dl(i) < 0.0 && (worst < 0 || dl(i) < dl(worst))) worst = i;
            if (worst >= 0) {
                for (int i = worst; i + 1 < m; ++i) cur[i] = cur[i + 1];
                --m;
                continue;
            }
            Vec4 bsum = Vec4::Zero();
            for (int i = 0; i < m; ++i) bsum += dl(i) * S[cur[i]].b;
            const Vec4 sig = trial - D * bsum;
            int add = -1;
            double fmax = tol;
            for (int k = 0; k < ns; ++k) {
                const double fk = S[k].a.dot(sig) - S[k].r;
                if (fk > fmax) { fmax = fk; add = k; }
            }
            if (add < 0 || m == 3) break;
            cur[m++] = add;
        }
    }
    const bool resolved = best >= 0.0;
    for (int i = 0; i < ns && !resolved; ++i) {
        J[0] = i;
        try_set(J, 1);
        for (int j = i + 1; j < ns; ++j) {
            J[1] = j;
            try_set(J, 2);
            for (int k = j + 1; k < ns; ++k) {
                J[2] = k;
                try_set(J, 3);
            }
        }
    }
    if (best < 0.0) {
        R.ok = false;
        return R;
    }
    R.plastic = true;
    R.sig = best_sig;
    R.n_active = best_n;
    for (int i = 0; i < 3; ++i) R.active[i] = best_set[i];
    // D - D B A^-1 A^T D, with the columns of B the flow directions and of A the gradients.
    Eigen::Matrix<double, 4, Eigen::Dynamic> DB(4, best_n), DA(4, best_n);
    for (int i = 0; i < best_n; ++i) {
        DB.col(i) = D * S[best_set[i]].b;
        DA.col(i) = D.transpose() * S[best_set[i]].a;
    }
    const Eigen::MatrixXd Ainv = best_Ainv.topLeftCorner(best_n, best_n);
    R.T = D - DB * Ainv * DA.transpose();
    return R;
}

}  // namespace katai::core::jr
