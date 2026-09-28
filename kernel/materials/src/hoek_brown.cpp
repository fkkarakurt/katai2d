#include <katai/materials/hoek_brown.hpp>

namespace katai::core::hoekbrown {

Return return_mapping(double s1t, double s2t, double s3t, const Params& P,
                      const Constants& C) {
    Return R{s1t, s2t, s3t, false, false, false};
    if (yield(s1t, s3t, P, C) <= 0.0) return R;   // elastic

    const LameConstants lame = lame_from(P.E, P.nu);
    const double lam = lame.lambda, mu = lame.mu;
    const double npsi = [&] {
        const double pm = psi_mobilised(s1t, P, C);
        const double sp = std::sin(pm);
        return (1.0 + sp) / (1.0 - sp);
    }();
    // D_e * m for m = [N_psi, 0, -1] (isotropic elasticity in principal space).
    const double d1 = (lam + 2.0 * mu) * npsi - lam;
    const double d2 = lam * (npsi - 1.0);
    const double d3 = lam * npsi - (lam + 2.0 * mu);

    // THE APEX IS A DOMAIN BOUNDARY, NOT A SPECIAL CASE. The criterion's bracket closes exactly at
    // s1 = sigma_t = s sigma_ci / m_b, so there is no surface on the tensile side of the
    // isotropic tension point. A trial that is already past it cannot be returned to a surface
    // that does not exist, and goes to the apex itself.
    if (s1t >= C.sigt) {
        R.s1 = R.s2 = R.s3 = C.sigt;
        R.plastic = true; R.apex = true;
        return R;
    }

    const auto f_of = [&](double L) {
        return yield(s1t - L * d1, s3t - L * d3, P, C);
    };
    // The corrector drives s1 down and s3 up, so f falls monotonically and without bound: a
    // doubling search finds a bracket on any input, and bisection cannot leave it. (Newton would
    // be faster and would step past the apex where the bracket's derivative vanishes; this is one
    // of the places where the slower method is the one that always answers.)
    double lo = 0.0, hi = 1.0;
    const double scale = std::fabs(s1t) + std::fabs(s3t) + std::fabs(C.sigc) + 1.0;
    hi = scale / std::max(1.0, d1 - d3);
    for (int i = 0; i < 200 && f_of(hi) > 0.0; ++i) hi *= 2.0;
    for (int i = 0; i < 80; ++i) {          // bisection: monotone in lambda, and it cannot escape
        const double mid = 0.5 * (lo + hi);
        (f_of(mid) > 0.0 ? lo : hi) = mid;
    }
    const double L = 0.5 * (lo + hi);
    R.s1 = s1t - L * d1;
    R.s2 = s2t - L * d2;
    R.s3 = s3t - L * d3;
    R.plastic = true;

    // THE ORDER IS PART OF THE ANSWER. The corrector moves s2 as well, so a state that started
    // ordered can leave the region where f_13 alone governs -- which is exactly where the second
    // yield function f_12 takes over. Clamping s2 back would report a stress the
    // criterion never admitted: measured against this tree's own Mohr-Coulomb return in the Tresca
    // degeneration, the clamp was 25 kPa out on a material whose cohesion is 50. So the edge is
    // returned to properly: BOTH surfaces active, two multipliers, two conditions.
    //
    //     sigma = sigma_tr - la * D m_a - lb * D m_b,   m_a = [N,0,-1],  m_b = [0,N,-1]
    //     f_13(sigma) = 0  and  f_12(sigma) = 0
    //
    // The directions are constant, so this is a 2x2 system in (la, lb) alone. It is solved by a
    // damped Newton on a numerical Jacobian -- two unknowns, smooth residuals, and the main-surface
    // solution as a starting point, which is inside the basin by construction.
    if (R.s2 > R.s1 + 1e-12 || R.s2 < R.s3 - 1e-12) {
        R.edge = true;
        const double e1 = lam * npsi - lam;                    // D m_b, component 1
        const double e2 = (lam + 2.0 * mu) * npsi - lam;       // component 2
        const double e3 = lam * npsi - (lam + 2.0 * mu);       // component 3
        double la = L, lb = 0.0;
        const auto residual = [&](double a, double b, double& r1, double& r2) {
            const double x1 = s1t - a * d1 - b * e1;
            const double x2 = s2t - a * d2 - b * e2;
            const double x3 = s3t - a * d3 - b * e3;
            r1 = yield(x1, x3, P, C);
            r2 = yield12(x2, x3, P, C);
        };
        double r1 = 0.0, r2 = 0.0;
        residual(la, lb, r1, r2);
        const double h = 1e-7 * (std::fabs(la) + std::fabs(lb) + 1e-12);
        for (int it = 0; it < 40; ++it) {
            if (std::fabs(r1) + std::fabs(r2) < 1e-10 * (std::fabs(C.sigc) + 1.0)) break;
            double a1, a2, b1, b2;
            residual(la + h, lb, a1, a2);
            residual(la, lb + h, b1, b2);
            const double J11 = (a1 - r1) / h, J12 = (b1 - r1) / h;
            const double J21 = (a2 - r2) / h, J22 = (b2 - r2) / h;
            const double det = J11 * J22 - J12 * J21;
            if (std::fabs(det) < 1e-300) break;
            const double da = (-r1 * J22 + r2 * J12) / det;
            const double db = (-r2 * J11 + r1 * J21) / det;
            double step = 1.0;
            for (int k = 0; k < 20; ++k) {                     // damping: never leave la, lb >= 0
                const double na = la + step * da, nb = lb + step * db;
                if (na >= 0.0 && nb >= 0.0) {
                    double q1, q2;
                    residual(na, nb, q1, q2);
                    if (std::fabs(q1) + std::fabs(q2) < std::fabs(r1) + std::fabs(r2)) {
                        la = na; lb = nb; r1 = q1; r2 = q2;
                        break;
                    }
                }
                step *= 0.5;
            }
        }
        R.s1 = s1t - la * d1 - lb * e1;
        R.s2 = s2t - la * d2 - lb * e2;
        R.s3 = s3t - la * d3 - lb * e3;
    }
    return R;
}

PlaneReturn plane_return(const PlaneStrainStress& trial, const Params& P,
                         const Constants& C) {
    enum { kA = 0, kB = 1, kZ = 2 };
    const double sxx = trial.in_plane(0), syy = trial.in_plane(1), sxy = trial.in_plane(2);
    const double mean = 0.5 * (sxx + syy), half = 0.5 * (sxx - syy);
    const double radius = std::sqrt(half * half + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = half / radius; sin2t = sxy / radius; }

    struct PV { double v; int src; };
    PV pv[3] = {{mean + radius, kA}, {mean - radius, kB}, {trial.zz, kZ}};
    std::sort(pv, pv + 3, [](const PV& x, const PV& y) { return x.v > y.v; });

    const Return R = return_mapping(pv[0].v, pv[1].v, pv[2].v, P, C);
    PlaneReturn out;
    out.plastic = R.plastic; out.apex = R.apex; out.edge = R.edge;
    if (!R.plastic) { out.stress = trial; return out; }

    const double ret[3] = {R.s1, R.s2, R.s3};
    double pa = 0.0, pb = 0.0, pz = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == kA) pa = ret[i];
        else if (pv[i].src == kB) pb = ret[i];
        else pz = ret[i];
    }
    const double m = 0.5 * (pa + pb), r = 0.5 * (pa - pb);
    out.stress.in_plane(0) = m + r * cos2t;
    out.stress.in_plane(1) = m - r * cos2t;
    out.stress.in_plane(2) = r * sin2t;
    out.stress.zz = pz;
    return out;
}

}  // namespace katai::core::hoekbrown
