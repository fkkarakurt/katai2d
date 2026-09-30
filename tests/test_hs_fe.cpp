// Hardening Soil -- P2.3d-3a: shear-hardening HS wired into the FE material point
// (integrate_point, plane strain). Verifies the principal-space wrapper: elastic
// predictor with the stress-dependent Eur, spectral decomposition of the trial stress,
// compression-positive sign conversion, shear return, and coaxial reconstruction.
//  (a) an elastic step reproduces the Eur predictor exactly,
//  (b) FRAME INDIFFERENCE: rotating the committed stress + strain increment rotates the
//      result identically (objectivity -> the eigen-decomposition/reconstruction is right),
//  (c) at large shear strain the principal stress difference is bounded by the MC qf.
// (See docs/references/hardening-soil-formulation.md and material_model.hpp hs_forward.)
#include <katai/materials/material_model.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

using katai::core::GaussState;
using katai::core::HardeningSoilParams;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::integrate_point;

namespace {

constexpr double kPi = 3.14159265358979323846;

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

MaterialModel make_hs() {
    MaterialModel m;
    m.type = MaterialType::HardeningSoil;
    m.hs.E50_ref = 3.0e4; m.hs.Eur_ref = 9.0e4; m.hs.Eoed_ref = 3.0e4;
    m.hs.m = 0.5; m.hs.p_ref = 100.0; m.hs.friction = 35.0 * kPi / 180.0;
    m.hs.dilatancy = 0.0; m.hs.Rf = 0.9; m.hs.nu_ur = 0.2; m.hs.cap_beta = 0.0;
    return m;
}

// Rotate an in-plane symmetric tensor by theta. Stress: [sxx,syy,sxy]; strain uses the
// engineering shear gamma = 2 eps_xy (factor handled by the caller's convention).
Eigen::Vector3d rotate_sym(const Eigen::Vector3d& v, double th, bool strain) {
    const double c = std::cos(th), s = std::sin(th);
    const double xy = strain ? 0.5 * v(2) : v(2);  // tensor shear
    const double xx = v(0), yy = v(1);
    const double xxr = xx * c * c + yy * s * s + 2 * xy * s * c;
    const double yyr = xx * s * s + yy * c * c - 2 * xy * s * c;
    const double xyr = (yy - xx) * s * c + xy * (c * c - s * s);
    return Eigen::Vector3d(xxr, yyr, strain ? 2.0 * xyr : xyr);
}

// Principal stress difference (tension-positive) and the minor principal (comp-positive).
void principals(const GaussState& s, double& q, double& sigma3_comp) {
    const double sxx = s.stress(0), syy = s.stress(1), sxy = s.stress(2);
    const double mean = 0.5 * (sxx + syy);
    const double radius = std::sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    const double a[3] = {mean + radius, mean - radius, s.stress_zz};
    const double smax = std::max(std::max(a[0], a[1]), a[2]);
    const double smin = std::min(std::min(a[0], a[1]), a[2]);
    q = smax - smin;
    sigma3_comp = -smax;  // most tensile tension-positive principal = comp-positive sigma3
}

void test_volumetric_then_shear() {
    // HS has no elastic deviatoric range (the shear surface passes through the origin):
    // any deviatoric loading is plastic. A volumetric (isotropic-stress) increment is
    // the only elastic one, but plane-strain isotropic STRAIN still produces a small
    // deviator, so even that hardens slightly -- this is correct HS behaviour. Verify
    // hardening grows monotonically under continued shear and stays bounded by qf(sigma3).
    const MaterialModel m = make_hs();
    GaussState s;
    s.stress = Eigen::Vector3d(-100, -100, 0);
    s.stress_zz = -100;
    double prev_gp = 0.0;
    bool monotone = true;
    for (int i = 0; i < 50; ++i) {
        GaussState tr; Eigen::Matrix3d tan;
        integrate_point(m, s, Eigen::Vector3d(0.0, -1e-4, 0.0), tr, tan);
        if (tr.gamma_p < prev_gp - 1e-14) monotone = false;
        prev_gp = tr.gamma_p;
        s = tr;
        double q, s3; principals(s, q, s3);
        if (q > m.hs.q_failure(s3) * (1.0 + 1e-3)) monotone = false;  // bound by qf(sigma3)
    }
    check(monotone && prev_gp > 0.0, "shear loading hardens monotonically, bounded by qf(sigma3)");
}

void test_frame_indifference() {
    const MaterialModel m = make_hs();
    // Anisotropic committed stress + a shear-inducing strain increment that yields.
    GaussState comm;
    comm.stress = Eigen::Vector3d(-160, -100, -20);
    comm.stress_zz = -110;
    comm.gamma_p = 0.0;
    const Eigen::Vector3d de(-1.5e-3, 4.0e-4, 8.0e-4);

    GaussState trA; Eigen::Matrix3d tanA;
    integrate_point(m, comm, de, trA, tanA);
    check(trA.gamma_p > 0.0, "frame test is at a plastic (yielding) state");

    const double th = 30.0 * kPi / 180.0;
    GaussState commR = comm;
    commR.stress = rotate_sym(comm.stress, th, false);
    const Eigen::Vector3d deR = rotate_sym(de, th, true);
    GaussState trB; Eigen::Matrix3d tanB;
    integrate_point(m, commR, deR, trB, tanB);

    // trB must equal the rotation of trA (objectivity).
    const Eigen::Vector3d expected = rotate_sym(trA.stress, th, false);
    const double err = (trB.stress - expected).cwiseAbs().maxCoeff() /
                       trA.stress.cwiseAbs().maxCoeff();
    std::printf("  frame indifference: rel err=%.3e  gamma_p=%.4e\n", err, trA.gamma_p);
    check(err < 1e-8, "HS is frame indifferent (rotation commutes with the update)");
    check(close(trB.stress_zz, trA.stress_zz, 1e-9), "sigma_zz unchanged by in-plane rotation");
    check(close(trB.gamma_p, trA.gamma_p, 1e-9), "hardening is frame independent");
}

void test_confined_loading_admissible() {
    // Confined compression (eps_yy down, eps_xx=0): sigma3 grows, the deviator hardens.
    // At every state the deviator must satisfy q <= qf(sigma3_current) (MC admissibility),
    // hardening must be monotone, and the stress must stay finite. (The shear<->MC-failure
    // coordination for paths that drive sigma3 toward tension is a later refinement,
    // alongside the cap, P2.3d-3b.)
    const MaterialModel m = make_hs();
    GaussState s;
    s.stress = Eigen::Vector3d(-100, -100, 0);
    s.stress_zz = -100;
    double prev_gp = 0.0;
    bool ok = true;
    double last_ratio = 0.0;
    for (int i = 0; i < 2000; ++i) {
        GaussState tr; Eigen::Matrix3d tan;
        integrate_point(m, s, Eigen::Vector3d(0.0, -5e-5, 0.0), tr, tan);
        if (tr.gamma_p < prev_gp - 1e-12) ok = false;
        prev_gp = tr.gamma_p;
        s = tr;
        if (!s.stress.allFinite()) ok = false;
        double q, s3; principals(s, q, s3);
        last_ratio = q / m.hs.q_failure(s3);
        if (q > m.hs.q_failure(s3) * (1.0 + 1e-3)) ok = false;
    }
    std::printf("  confined loading: gamma_p=%.4e  final q/qf(sigma3)=%.3f\n",
                prev_gp, last_ratio);
    check(ok && prev_gp > 0.0, "confined loading: monotone hardening, q<=qf(sigma3), finite");
}

// (d) THE TENSION CUT-OFF IS IN THE TANGENT. A shallow point pulled into tension is capped after
// the model's own return; the tangent handed to Newton has to be the derivative of the capped
// stress, not of the one before the cap. Until 2026-09-29 it was the latter: a central difference
// of the whole update disagreed with it at EVERY capped point, by up to 1e5 relative -- the row
// the cap pins to sigma_t still carried the full elastic stiffness. One and two capped
// principals, each well inside its active set so the difference does not straddle a switch.
// With c > 0: a cohesionless point cannot carry tension under Mohr-Coulomb in the first place, so
// the model's own return already stops it at the apex and the cut-off has nothing left to do.
void test_tension_cutoff_tangent() {
    MaterialModel m = make_hs();
    m.hs.cohesion = 10.0;
    m.tension_cutoff = true;
    m.tensile_strength = 0.0;
    struct Case { Eigen::Vector3d de; int ncap; const char* what; };
    const Case cases[] = {
        {Eigen::Vector3d(4e-4, -1e-4, 0.5e-4), 1, "one principal capped: tangent = d(capped)/d(eps)"},
        {Eigen::Vector3d(2e-4, 2e-4, 0.2e-4), 2, "two principals capped: tangent = d(capped)/d(eps)"},
    };
    for (const Case& c : cases) {
        GaussState g;
        g.stress = Eigen::Vector3d(-3.0, -4.0, 0.3);   // tension positive: a shallow point
        g.stress_zz = -8.0;
        const Eigen::Vector3d sc(8.0, 4.0, 3.0);
        g.pp = katai::core::hs_initial_pp(m.hs, sc);
        g.gamma_p = katai::core::hs_initial_gamma_p(m.hs, sc);
        GaussState t;
        Eigen::Matrix3d D;
        bool pl = false;
        katai::core::hs_forward(m, g, c.de, t, &D, nullptr, &pl);
        int ncap = 0;
        for (double s : {0.5 * (t.stress(0) + t.stress(1)) +
                             std::hypot(0.5 * (t.stress(0) - t.stress(1)), t.stress(2)),
                         0.5 * (t.stress(0) + t.stress(1)) -
                             std::hypot(0.5 * (t.stress(0) - t.stress(1)), t.stress(2)),
                         t.stress_zz})
            ncap += std::fabs(s) < 1e-9 ? 1 : 0;
        Eigen::Matrix3d F;
        const double h = 1e-7 * c.de.norm();
        for (int k = 0; k < 3; ++k) {
            Eigen::Vector3d e = Eigen::Vector3d::Zero();
            e(k) = h;
            GaussState tp, tm;
            katai::core::hs_forward(m, g, c.de + e, tp);
            katai::core::hs_forward(m, g, c.de - e, tm);
            F.col(k) = (tp.stress - tm.stress) / (2.0 * h);
        }
        // Relative to the elastic stiffness where the capped block is zero: with both in-plane
        // principals pinned to sigma_t the in-plane tangent IS zero, and so is the difference.
        const double err = (D - F).norm() / std::max(F.norm(), m.hs.Eur(m.hs.p_limit()));
        std::printf("  tension cut-off tangent: %d capped, |D - FD|/|FD| = %.2e\n", ncap, err);
        check(pl && ncap == c.ncap && err < 1e-3, c.what);
    }
}

// (e) NO STRENGTH BELOW THE STIFFNESS FLOOR. p_limit (= 0.1 p_ref) keeps E_i and E_ur off zero; it
// is not part of the failure criterion. A cohesionless point sheared at a minor stress of 2 kPa
// must stop at the Mohr-Coulomb deviator of its own stress, q_f = 2 sigma3 sin(phi)/(1 - sin(phi))
// = 5.38 kPa at phi = 35 -- not at q_f(p_limit) = 26.9 kPa, which is where it stopped until
// 2026-09-29: an apparent cohesion of up to 7 kPa on the unsafe side, at exactly the depth a
// footing's bearing capacity is decided. The vertical strain is driven and the horizontal strain
// chosen at every step so that sigma_xx stays at -2 kPa.
void test_strength_is_not_floored() {
    const MaterialModel m = make_hs();
    const double s3 = 2.0;
    GaussState s;
    s.stress = Eigen::Vector3d(-s3, -s3, 0.0);
    s.stress_zz = -s3;
    double qmax = 0.0;
    bool finite = true;
    for (int i = 0; i < 4000; ++i) {
        // Pick eps_xx by two secant corrections so that sigma_xx stays at -s3.
        double exx = 0.0;
        GaussState tr; Eigen::Matrix3d tan;
        for (int it = 0; it < 30; ++it) {
            integrate_point(m, s, Eigen::Vector3d(exx, -2e-6, 0.0), tr, tan);
            const double r = tr.stress(0) + s3;
            if (std::fabs(r) < 1e-10) break;
            exx -= r / tan(0, 0);
        }
        s = tr;
        if (!s.stress.allFinite()) finite = false;
        double q, sm; principals(s, q, sm);
        qmax = std::max(qmax, q);
    }
    const double qf = m.hs.q_failure(s3);
    std::printf("  shallow cohesionless shear: q_max = %.3f kPa, q_f(%.0f kPa) = %.3f, q_f(p_limit) = %.3f\n",
                qmax, s3, qf, m.hs.q_failure(m.hs.p_limit()));
    check(finite && qmax <= qf * (1.0 + 1e-3) + 1e-9 && qmax >= 0.97 * qf,
          "a cohesionless point below p_limit fails at the Mohr-Coulomb deviator of its own stress");
}

// (f) THE MOBILISED DILATANCY RULE OF HARDENING SOIL. Rowe's law only above a mobilised friction of
// 3/4 sin(phi); zero below it; a non-positive psi taken as it is above it; HSsmall not affected
// (its Li & Dafalias branch is pinned in test_hssmall). phi = psi = 41 degrees puts phi_cv at 0,
// where the rule without the threshold dilated from the first increment of shear.
void test_mobilised_dilatancy_rule() {
    auto rule = [](double phi_deg, double psi_deg) {
        HardeningSoilParams h;
        h.friction = phi_deg * kPi / 180.0;
        h.dilatancy = psi_deg * kPi / 180.0;
        return katai::core::detail::hs_dilatancy(h);
    };
    const double s41 = std::sin(41.0 * kPi / 180.0);
    const auto d41 = rule(41.0, 41.0);
    const bool below = d41(0.74 * s41) == 0.0;
    const bool above = std::fabs(d41(0.76 * s41) - 0.76 * s41) < 1e-12;   // phi_cv = 0: Rowe = sin phi_m
    const auto dneg = rule(35.0, -5.0);
    const double s35 = std::sin(35.0 * kPi / 180.0);
    const bool neg_below = dneg(0.7 * s35) == 0.0;
    const bool neg_above = std::fabs(dneg(0.9 * s35) - std::sin(-5.0 * kPi / 180.0)) < 1e-15;
    const auto dsmall = rule(35.0, 5.0);   // phi_cv above the threshold: the threshold changes nothing
    const double scv = dsmall.sin_cs;
    const bool small_ok = dsmall(0.99 * scv) == 0.0 &&
                          std::fabs(dsmall(0.5 * (scv + s35)) -
                                    (0.5 * (scv + s35) - scv) / (1.0 - 0.5 * (scv + s35) * scv)) < 1e-15;
    std::printf("  dilatancy rule: phi=psi=41 below/above 3/4 sin(phi): %.3f / %.3f; psi=-5 above: %.4f\n",
                d41(0.74 * s41), d41(0.76 * s41), dneg(0.9 * s35));
    check(below && above, "HS: psi_m = 0 below 3/4 sin(phi), Rowe above it");
    check(neg_below && neg_above, "HS: a negative psi is psi_m above the threshold, zero below it");
    check(small_ok, "HS: where phi_cv lies above the threshold the rule is Rowe's, unchanged");
}

} // namespace

int main() {
    test_volumetric_then_shear();
    test_frame_indifference();
    test_confined_loading_admissible();
    test_tension_cutoff_tangent();
    test_strength_is_not_floored();
    test_mobilised_dilatancy_rule();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil FE material point (shear) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
