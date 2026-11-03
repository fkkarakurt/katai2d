// Hardening Soil -- robust multi-surface integrator (hs_integrate: explicit substepping
// + Koiter active-set, Sloan/Potts&Gens). Replaces the fragile implicit both-active
// solver. Verified to (1) reproduce the shear hyperbola (cap off) in a drained triaxial,
// (2) give the cap K_iso in isotropic compression, and (3) integrate the Berlin Sand III
// oedometer (E=105 MPa) ROBUSTLY -- the high-stiffness case where the previous implicit
// solver diverged. (See docs/references/hardening-soil-formulation.md.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::HsIntegrated;
using katai::core::hs_integrate;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

// One drained-triaxial step (sigma3 const): solve lateral strain via the tangent.
void triaxial_step(const HardeningSoilParams& p, double sigma3, double de1,
                   Eigen::Vector3d& sig, double& gp, double& pp, double& eps1) {
    double dlat = 0.0; HsIntegrated r;
    for (int it = 0; it < 25; ++it) {
        r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, dlat, dlat));
        const double g = r.stress(2) - sigma3;
        if (std::fabs(g) <= 1e-9 * (1.0 + std::fabs(sigma3))) break;
        const double dgdl = r.tangent(2, 1) + r.tangent(2, 2);
        dlat -= g / dgdl;
    }
    sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += de1;
}

void test_triaxial_hyperbola() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 35.0 * kPi / 180.0; p.dilatancy = 0.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.cap_beta = 0.0;  // cap off -> pure shear
    const double sigma3 = 100.0;
    const double Ei = p.Ei(sigma3), qa = p.q_asymptote(sigma3), qf = p.q_failure(sigma3);
    auto hyper = [&](double q) { return (1.0 / Ei) * q / (1.0 - q / qa); };

    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0, pp = 0, eps1 = 0, maxrel = 0; int n = 0; double frac = 0.25;
    for (int s = 0; s < 6000; ++s) {
        triaxial_step(p, sigma3, 5e-6, sig, gp, pp, eps1);
        const double q = sig(0) - sig(2);
        if (q >= frac * qf && frac < 0.76) {
            maxrel = std::fmax(maxrel, std::fabs(eps1 - hyper(q)) / hyper(q));
            ++n; frac += 0.25;
        }
    }
    std::printf("  triaxial (substep): checks=%d  max rel err vs hyperbola=%.3e\n", n, maxrel);
    check(n == 3 && maxrel < 1.5e-3, "robust integrator reproduces the shear hyperbola");
}

void test_isotropic_cap() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.m = 0.5; p.p_ref = 100.0;
    p.friction = 30.0 * kPi / 180.0; p.nu_ur = 0.2;
    p.cap_alpha = 0.6; p.cap_beta = 2.0e-3;
    const double p0 = 100.0, e = 1e-6;
    Eigen::Vector3d sig(p0, p0, p0);
    const HsIntegrated r = hs_integrate(p, sig, 0.0, p0, Eigen::Vector3d(e, e, e));
    const double Ke = p.bulk(p0), Hcap = p.cap_hardening_modulus(r.pp);
    const double K_iso = Ke * Hcap / (Ke + Hcap);
    const double dp = (r.stress(0) + r.stress(1) + r.stress(2)) / 3.0 - p0;
    std::printf("  isotropic (substep): K_iso=%.1f  return secant=%.1f\n", K_iso, dp / (3 * e));
    check(close(dp / (3 * e), K_iso, 5e-3), "robust integrator gives cap K_iso");
}

void test_berlin_sand_oedometer_robust() {
    // Berlin Sand III: the HIGH-stiffness case that broke
    // the previous implicit both-active solver. Here it must integrate robustly.
    HardeningSoilParams p;
    p.p_ref = 100.0; p.E50_ref = 105.0e3; p.Eur_ref = 315.0e3; p.Eoed_ref = 105.0e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38.0 * kPi / 180.0;
    p.dilatancy = 6.0 * kPi / 180.0; p.cohesion = 1.0; p.Rf = 0.9;
    p.cap_alpha = 1.0; p.cap_beta = 1.0e-3;  // any reasonable cap

    Eigen::Vector3d sig(5.0, 5.0, 5.0);
    double gp = 0, pp = 5.0, eps1 = 0;
    bool finite = true, monotone = true; double prev_gp = 0, K0 = 0, Eoed = 0, s1b = 5.0;
    for (int i = 0; i < 2000; ++i) {
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(5e-5, 0, 0));
        if (!r.stress.allFinite()) { finite = false; break; }
        if (r.gamma_p < prev_gp - 1e-10) monotone = false;
        prev_gp = r.gamma_p;
        if (s1b < p.p_ref && r.stress(0) >= p.p_ref) {
            Eoed = (r.stress(0) - s1b) / 5e-5; K0 = r.stress(2) / r.stress(0);
        }
        s1b = r.stress(0);
        sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += 5e-5;
        if (sig(0) > 3.0 * p.p_ref) break;
    }
    std::printf("  Berlin Sand oedometer: finite=%d monotone=%d  K0=%.3f Eoed=%.0f\n",
                finite, monotone, K0, Eoed);
    check(finite && monotone, "Berlin Sand III oedometer integrates robustly (no divergence)");
    check(K0 > 0.2 && K0 < 0.8 && Eoed > 0.0, "Berlin Sand oedometer gives physical K0, Eoed");
}

// THE IMPLICIT RETURN OF THE STIFFNESS FLOOR (hs_integrate): a point of a fill with a cohesion,
// pulled into tension across the cut-off. The material and the increments are three integrations
// of the first phase of an embankment on soft soil, as that run handed them over. Checked:
// (1) a point pulled from rest -- whose walk ran to its ceiling -- is returned in one step, and
//     agrees with the limit of the walk to 1e-4: uniform walks of 20 000 and 40 000 substeps
//     (replayed plans: no error control, no implicit return), extrapolated to an infinite count as
//     the first-order error the walk has there (R falling as 1/n) extrapolates;
// (2) its tangent is the derivative of the return, by central differences of replays of it;
// (3) an increment whose answer depends on its path -- a cut-off plane met and left again in the
//     first 0.3% of it -- carries the one step's error, 4 Pa in a stress of 90 Pa (under 1e-3 of
//     the 5 kPa the tolerance is read against), and the return is CONSISTENT: chained over 1024
//     equal parts of the increment it converges onto the walk's limit.
void test_floor_implicit_return() {
    HardeningSoilParams p;
    p.E50_ref = 25.0e3; p.Eur_ref = 75.0e3; p.Eoed_ref = 25.0e3; p.m = 0.5; p.p_ref = 100.0;
    p.cohesion = 1.0; p.friction = 30.0 * kPi / 180.0; p.dilatancy = 0.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.cap_alpha = 1.2184094641781964; p.cap_beta = 0.0019959451779286543; p.K0nc = 0.5;
    p.sigma_t = 0.0;
    auto uniform_walk = [&](const Eigen::Vector3d& s, double gp, double pp, const Eigen::Vector3d& de) {
        auto at = [&](int n) {
            katai::core::HsSubstepPlan plan;
            plan.n = n;
            plan.dT = 1.0 / n;
            return hs_integrate(p, s, gp, pp, de, 0.0, nullptr, &plan);
        };
        const HsIntegrated a = at(20000), b = at(40000);
        HsIntegrated lim = b;
        lim.stress = 2.0 * b.stress - a.stress;
        lim.gamma_p = 2.0 * b.gamma_p - a.gamma_p;
        lim.pp = 2.0 * b.pp - a.pp;
        return lim;
    };
    // the measure two successive walks are compared in: the stress over its 5 kPa floor, the two
    // hardening variables relative to themselves
    auto differ = [](const HsIntegrated& a, const HsIntegrated& b) {
        const double ds = (a.stress - b.stress).norm() / std::max(b.stress.norm(), 5.0);
        const double dg = std::fabs(a.gamma_p - b.gamma_p) / std::max(std::fabs(b.gamma_p), 1e-12);
        const double dp = std::fabs(a.pp - b.pp) / std::max(std::fabs(b.pp), 1e-12);
        return std::max(ds, std::max(dg, dp));
    };
    {
        const Eigen::Vector3d s0 = Eigen::Vector3d::Zero();
        const Eigen::Vector3d de(0.000106361899975711, 3.3881317890172014e-21, -0.00027398110327038483);
        katai::core::HsSubstepPlan plan;
        const HsIntegrated r = hs_integrate(p, s0, 0.0, 0.0, de, 0.0, &plan);
        const HsIntegrated w = uniform_walk(s0, 0.0, 0.0, de);
        std::printf("  floor, pulled from rest: implicit nsub=%d s=(%.6f %.2e %.2e) gp=%.6e pp=%.6f\n"
                    "                           walk limit   s=(%.6f %.2e %.2e) gp=%.6e pp=%.6f\n",
                    r.nsub, r.stress(0), r.stress(1), r.stress(2), r.gamma_p, r.pp, w.stress(0), w.stress(1),
                    w.stress(2), w.gamma_p, w.pp);
        check(r.nsub == 1 && plan.implicit, "a point pulled from rest across the cut-off is returned implicitly");
        std::printf("  floor, pulled from rest: implicit vs walk limit %.2e\n", differ(r, w));
        check(differ(r, w) < 1e-4, "the implicit return agrees with the limit of the walk");
        // (2) the tangent, by central differences of replays of the implicit return
        double worst = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double h = 1e-7 * de.norm();
            Eigen::Vector3d dp = de, dm = de;
            dp(j) += h;
            dm(j) -= h;
            const HsIntegrated rp = hs_integrate(p, s0, 0.0, 0.0, dp, 0.0, nullptr, &plan);
            const HsIntegrated rm = hs_integrate(p, s0, 0.0, 0.0, dm, 0.0, nullptr, &plan);
            const Eigen::Vector3d col = (rp.stress - rm.stress) / (2.0 * h);
            worst = std::max(worst, (col - r.tangent.col(j)).norm() / std::max(r.tangent.norm(), 1.0));
        }
        std::printf("  floor implicit tangent vs central differences: %.2e\n", worst);
        check(worst < 1e-5, "the implicit return's tangent is its derivative");
    }
    {
        const Eigen::Vector3d s0(0.030517498771983995, 0.0009718762049324102, 0.021941964441669203);
        const double gp0 = 0.00012872319607270638, pp0 = 1.0494678337865455;
        const Eigen::Vector3d de(1.6718255103433506e-05, -5.082197683525802e-21, -0.0002446996066676391);
        katai::core::HsSubstepPlan plan;
        const HsIntegrated r = hs_integrate(p, s0, gp0, pp0, de, 0.0, &plan);
        const HsIntegrated w = uniform_walk(s0, gp0, pp0, de);
        // the increment in 1024 equal parts, each integrated as hs_integrate integrates it
        constexpr int kParts = 1024;
        HsIntegrated c = r;
        c.stress = s0; c.gamma_p = gp0; c.pp = pp0;
        for (int k = 0; k < kParts; ++k) c = hs_integrate(p, c.stress, c.gamma_p, c.pp, de / kParts);
        std::printf("  floor, path-dependent: nsub=%d s=(%.6f %.6f %.2e) | in %d parts (%.6f %.6f) | walk limit "
                    "s=(%.6f %.6f %.2e)\n  one step %.2e, %d parts %.2e from the limit\n",
                    r.nsub, r.stress(0), r.stress(1), r.stress(2), kParts, c.stress(0), c.stress(1), w.stress(0),
                    w.stress(1), w.stress(2), differ(r, w), kParts, differ(c, w));
        check(r.nsub == 1 && plan.implicit, "the path-dependent increment is returned implicitly too");
        check(differ(r, w) < 1e-3, "one implicit step is within 1e-3 of the walk's limit");
        check(differ(c, w) < 0.1 * differ(r, w), "the implicit return converges onto the walk's limit");
    }
    // (3) soft clay at the surface of a slope under fill, pulled sideways: both lateral stresses end on
    // their cut-off planes with the shear face and the cap (the edge). Read with one face, the minor's
    // plane takes a negative multiplier; the edge's other face gives it, and the set holds.
    {
        HardeningSoilParams q;
        q.E50_ref = 5600.0; q.Eur_ref = 20000.0; q.Eoed_ref = 5000.0; q.m = 1.0; q.p_ref = 100.0;
        q.cohesion = 10.0; q.friction = 25.0 * kPi / 180.0; q.dilatancy = 0.0; q.Rf = 0.9; q.nu_ur = 0.2;
        q.cap_alpha = 0.8981985941422174; q.cap_beta = 0.01045656517386395; q.K0nc = 0.5773817382593005;
        q.sigma_t = 0.0;
        const Eigen::Vector3d s0(1.648555476431281, 6.576204169810876e-06, 0.00034263449314631683);
        const double gp0 = 0.000676145769572623, pp0 = 1.9162882618444304;
        const Eigen::Vector3d de(0.0008218384100099963, -2.710505431213761e-20, -0.0008222232803236482);
        katai::core::HsSubstepPlan plan;
        const HsIntegrated r = hs_integrate(q, s0, gp0, pp0, de, 0.0, &plan);
        auto at = [&](int n) {
            katai::core::HsSubstepPlan u;
            u.n = n;
            u.dT = 1.0 / n;
            return hs_integrate(q, s0, gp0, pp0, de, 0.0, nullptr, &u);
        };
        const HsIntegrated a = at(20000), b = at(40000);
        HsIntegrated w = b;
        w.stress = 2.0 * b.stress - a.stress;
        w.gamma_p = 2.0 * b.gamma_p - a.gamma_p;
        w.pp = 2.0 * b.pp - a.pp;
        std::printf("  floor, two-plane edge: implicit nsub=%d s=(%.6f %.2e %.2e) gp=%.6e pp=%.6f\n"
                    "                         walk limit   s=(%.6f %.2e %.2e) gp=%.6e pp=%.6f\n",
                    r.nsub, r.stress(0), r.stress(1), r.stress(2), r.gamma_p, r.pp, w.stress(0), w.stress(1),
                    w.stress(2), w.gamma_p, w.pp);
        check(r.nsub == 1 && plan.implicit, "a point on the two-plane edge is returned implicitly");
        std::printf("  floor, two-plane edge: implicit vs walk limit %.2e\n", differ(r, w));
        check(differ(r, w) < 1e-4, "the edge's implicit return agrees with the limit of the walk");
        double worst = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double h = 1e-7 * de.norm();
            Eigen::Vector3d dp = de, dm = de;
            dp(j) += h;
            dm(j) -= h;
            const HsIntegrated rp = hs_integrate(q, s0, gp0, pp0, dp, 0.0, nullptr, &plan);
            const HsIntegrated rm = hs_integrate(q, s0, gp0, pp0, dm, 0.0, nullptr, &plan);
            const Eigen::Vector3d col = (rp.stress - rm.stress) / (2.0 * h);
            worst = std::max(worst, (col - r.tangent.col(j)).norm() / std::max(r.tangent.norm(), 1.0));
        }
        std::printf("  floor, two-plane edge: tangent vs central differences %.2e\n", worst);
        check(worst < 1e-5, "the edge's tangent is its derivative");
    }
    // (4) fill at the toe of an embankment, at the compression corner (sigma_b = sigma_c, off the cut-off),
    // pulled sideways: the face (a, c) alone turns the ordering over; with the face (a, b) the two
    // stay equal, as the walk keeps them.
    {
        HardeningSoilParams q;
        q.E50_ref = 15000.0; q.Eur_ref = 45000.0; q.Eoed_ref = 15000.0; q.m = 0.5; q.p_ref = 100.0;
        q.cohesion = 3.0; q.friction = 30.0 * kPi / 180.0; q.dilatancy = 0.0; q.Rf = 0.9; q.nu_ur = 0.2;
        q.cap_alpha = 1.2738860785916148; q.cap_beta = 0.0034323938176975257; q.K0nc = 0.5;
        q.sigma_t = 0.0;
        const Eigen::Vector3d s0(4.340165321101571, 0.04527568332118562, 0.04538116750829531);
        const double gp0 = 0.0008105400300646346, pp0 = 3.6808739598287943;
        const Eigen::Vector3d de(0.0002530522252955523, 0.0, -0.00018579331663896972);
        katai::core::HsSubstepPlan plan;
        const HsIntegrated r = hs_integrate(q, s0, gp0, pp0, de, 0.0, &plan);
        auto at = [&](int n) {
            katai::core::HsSubstepPlan u;
            u.n = n;
            u.dT = 1.0 / n;
            return hs_integrate(q, s0, gp0, pp0, de, 0.0, nullptr, &u);
        };
        const HsIntegrated a = at(20000), b = at(40000);
        HsIntegrated w = b;
        w.stress = 2.0 * b.stress - a.stress;
        w.gamma_p = 2.0 * b.gamma_p - a.gamma_p;
        w.pp = 2.0 * b.pp - a.pp;
        std::printf("  floor, compression corner: implicit nsub=%d s=(%.6f %.6f %.6f) gp=%.6e pp=%.6f\n"
                    "                             walk limit   s=(%.6f %.6f %.6f) gp=%.6e pp=%.6f\n",
                    r.nsub, r.stress(0), r.stress(1), r.stress(2), r.gamma_p, r.pp, w.stress(0), w.stress(1),
                    w.stress(2), w.gamma_p, w.pp);
        check(r.nsub == 1 && plan.implicit, "a point at the compression corner is returned implicitly");
        check(std::fabs(r.stress(1) - r.stress(2)) < 1e-9 * (1.0 + r.stress.norm()),
              "the two stresses of the corner stay equal");
        std::printf("  floor, compression corner: implicit vs walk limit %.2e\n", differ(r, w));
        check(differ(r, w) < 1e-3, "the corner's implicit return agrees with the limit of the walk");
        double worst = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double h = 1e-7 * de.norm();
            Eigen::Vector3d dp = de, dm = de;
            dp(j) += h;
            dm(j) -= h;
            const HsIntegrated rp = hs_integrate(q, s0, gp0, pp0, dp, 0.0, nullptr, &plan);
            const HsIntegrated rm = hs_integrate(q, s0, gp0, pp0, dm, 0.0, nullptr, &plan);
            const Eigen::Vector3d col = (rp.stress - rm.stress) / (2.0 * h);
            worst = std::max(worst, (col - r.tangent.col(j)).norm() / std::max(r.tangent.norm(), 1.0));
        }
        std::printf("  floor, compression corner: tangent vs central differences %.2e\n", worst);
        check(worst < 1e-5, "the corner's tangent is its derivative");
    }
}

// A POINT WITHOUT COHESION THAT ENDS AT THE APEX: sand (c = 0, phi 32, psi 2, cut-off at zero) in the
// stiffness floor, on the failure plateau, given a strain whose plastic part the apex can take --
// compression in one direction with more than (1 + sin psi) / (1 - sin psi) times as much extension
// in the others. The answer is the apex, its hardening untouched: the walk slid there along the
// cut-off with the last substep's elastic part left over, 0.0348 / 0.0035 / 0.0004 kPa at 1 000 /
// 10 000 / 40 000 substeps, and ran every call to the ceiling (the ground heaved beside a wall pushed
// into the sand: one to five seconds a call). With the extension too small for the apex the point
// is integrated as before, and is not at the apex.
void test_apex_without_cohesion() {
    HardeningSoilParams p;
    p.E50_ref = 40.0e3; p.Eur_ref = 120.0e3; p.Eoed_ref = 40.0e3; p.m = 0.5; p.p_ref = 100.0;
    p.cohesion = 0.0; p.friction = 32.0 * kPi / 180.0; p.dilatancy = 2.0 * kPi / 180.0; p.Rf = 0.9;
    p.nu_ur = 0.2; p.cap_alpha = 1.324586561619468; p.cap_beta = 0.0012449161950475467; p.K0nc = 0.47;
    p.sigma_t = 0.0;
    const Eigen::Vector3d s0(5.059234277866875, 5.330585043705342, 1.9092184612089125);
    const double gp0 = 0.024436522062272314, pp0 = 68.47648619243508;
    const Eigen::Vector3d de(0.0043973899224799395, 1.0842021724855044e-19, -0.005432015299931425);
    const HsIntegrated r = hs_integrate(p, s0, gp0, pp0, de);
    std::printf("  apex, c = 0: nsub=%d s=(%.3e %.3e %.3e) gp=%.6e pp=%.6f\n", r.nsub, r.stress(0), r.stress(1),
                r.stress(2), r.gamma_p, r.pp);
    check(r.nsub <= 1 && r.saturated == 0, "a point that ends at the apex is returned there at once");
    check(r.stress.cwiseAbs().maxCoeff() == 0.0, "its stress is the apex");
    check(r.gamma_p == gp0 && std::fabs(r.pp - pp0) <= 1e-9 * pp0, "its hardening is untouched");
    // the extension a quarter of the compression: not the apex's
    const Eigen::Vector3d de2(0.0043973899224799395, 0.0, -0.25 * 0.0043973899224799395);
    const HsIntegrated r2 = hs_integrate(p, s0, gp0, pp0, de2);
    std::printf("  apex, c = 0, mostly compressed: nsub=%d s=(%.3f %.3f %.3f)\n", r2.nsub, r2.stress(0),
                r2.stress(1), r2.stress(2));
    check(r2.stress.cwiseAbs().maxCoeff() > 1.0, "a strain the apex cannot take is integrated, off the apex");
}

} // namespace

int main() {
    test_triaxial_hyperbola();
    test_isotropic_cap();
    test_berlin_sand_oedometer_robust();
    test_floor_implicit_return();
    test_apex_without_cohesion();
    if (g_failures == 0) {
        std::printf("OK: robust multi-surface HS integrator (substepping + Koiter) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
