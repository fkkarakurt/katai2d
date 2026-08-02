// Verification of the Mohr-Coulomb constitutive kernels (P1.2a).
//
// These tests exercise the building blocks of the return-mapping algorithm in
// isolation: the plane-strain elastic predictor (including the out-of-plane
// reaction sigma_zz), the spectral decomposition of the stress tensor, and the
// yield function. Each is checked against closed-form expectations so that any
// sign or convention error is caught before the return mapping is layered on.
#include <katai/materials/material_model.hpp>
#include <katai/materials/mohr_coulomb.hpp>

#include <cmath>
#include <cstdio>

using katai::core::elastic_predictor;
using katai::core::GaussState;
using katai::core::integrate_point;
using katai::core::lame_from;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::LameConstants;
using katai::core::mc_return_mapping;
using katai::core::mc_yield;
using katai::core::McReturn;
using katai::core::MohrCoulombParams;
using katai::core::PlaneStrainStress;
using katai::core::principal_stresses;
using katai::core::PrincipalStresses;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

constexpr double kE = 10000.0;
constexpr double kNu = 0.3;

// The elastic predictor must reproduce the plane-strain Hooke response: an
// in-plane axial strain induces both an in-plane Poisson stress and an
// out-of-plane reaction d_sigma_zz = lambda (d_eps_xx + d_eps_yy).
void test_elastic_predictor() {
    const LameConstants lame = lame_from(kE, kNu);
    const double lam = lame.lambda, mu = lame.mu;

    // (a) Uniaxial strain in x: only d_eps_xx is non-zero.
    {
        const double exx = 1.0e-3;
        const auto s =
            elastic_predictor({}, Eigen::Vector3d(exx, 0.0, 0.0), lame);
        check(close(s.in_plane(0), (lam + 2.0 * mu) * exx), "uniaxial sxx");
        check(close(s.in_plane(1), lam * exx), "uniaxial syy (Poisson)");
        check(close(s.in_plane(2), 0.0), "uniaxial sxy = 0");
        check(close(s.zz, lam * exx), "uniaxial szz reaction");
    }
    // (b) Pure shear: trace is zero, so sigma_zz is unaffected.
    {
        const double gxy = 2.0e-3;
        const auto s =
            elastic_predictor({}, Eigen::Vector3d(0.0, 0.0, gxy), lame);
        check(close(s.in_plane(2), mu * gxy), "shear sxy = mu*gamma");
        check(close(s.in_plane(0), 0.0) && close(s.in_plane(1), 0.0),
              "shear normal stresses = 0");
        check(close(s.zz, 0.0), "shear szz = 0");
    }
    // (c) Incremental superposition from a non-zero base state.
    {
        PlaneStrainStress base;
        base.in_plane = Eigen::Vector3d(-40.0, -20.0, 5.0);
        base.zz = -18.0;
        const auto s =
            elastic_predictor(base, Eigen::Vector3d(0.0, 0.0, 0.0), lame);
        check(close(s.in_plane(0), -40.0) && close(s.zz, -18.0),
              "zero increment leaves state unchanged");
    }
}

// Spectral decomposition: sigma_zz is one principal; the other two follow from
// the in-plane Mohr circle. Validate ordering and rotational invariance.
void test_principal_stresses() {
    PlaneStrainStress s;
    s.in_plane = Eigen::Vector3d(-100.0, -50.0, 30.0);
    s.zz = -60.0;
    const auto p = principal_stresses(s);

    const double radius = std::sqrt(25.0 * 25.0 + 30.0 * 30.0);  // 39.0512...
    const double in_a = -75.0 + radius;   // -35.9488
    const double in_b = -75.0 - radius;   // -114.0512
    check(close(p.s1, in_a), "sigma_1 = larger in-plane principal");
    check(close(p.s2, -60.0), "sigma_2 = sigma_zz (intermediate)");
    check(close(p.s3, in_b), "sigma_3 = smaller in-plane principal");
    check(p.s1 >= p.s2 && p.s2 >= p.s3, "principals sorted descending");

    // Rotational invariance: rotating the in-plane axes must not change the
    // principal magnitudes.
    const double th = 20.0 * 3.14159265358979323846 / 180.0;
    const double c2 = std::cos(2.0 * th), s2 = std::sin(2.0 * th);
    const double mean = -75.0, hdiff = -25.0, sxy = 30.0;
    PlaneStrainStress r;
    r.in_plane = Eigen::Vector3d(mean + hdiff * c2 + sxy * s2,
                                 mean - hdiff * c2 - sxy * s2,
                                 -hdiff * s2 + sxy * c2);
    r.zz = -60.0;
    const auto pr = principal_stresses(r);
    check(close(pr.s1, p.s1) && close(pr.s2, p.s2) && close(pr.s3, p.s3),
          "principal stresses invariant under in-plane rotation");
}

// Yield function: zero on a state placed exactly on the failure surface,
// strictly negative for hydrostatic compression (no deviatoric stress).
void test_yield_function() {
    const double phi = 30.0 * 3.14159265358979323846 / 180.0;
    const double c = 10.0;
    const double sphi = std::sin(phi), cphi = std::cos(phi);

    // Place sigma_3 = -100 and solve f = 0 for sigma_1; sigma_2 in between.
    const double s3 = -100.0;
    const double s1 = (s3 * (1.0 - sphi) + 2.0 * c * cphi) / (1.0 + sphi);
    PlaneStrainStress on_surface;
    on_surface.in_plane = Eigen::Vector3d(s1, s3, 0.0);  // axis-aligned
    on_surface.zz = -60.0;                               // intermediate
    check(close(mc_yield(principal_stresses(on_surface), c, phi), 0.0, 1e-9),
          "yield = 0 on the failure surface");

    PlaneStrainStress hydro;
    hydro.in_plane = Eigen::Vector3d(-50.0, -50.0, 0.0);
    hydro.zz = -50.0;
    check(mc_yield(principal_stresses(hydro), c, phi) < 0.0,
          "hydrostatic compression is admissible (f < 0)");
}

double deg(double d) { return d * 3.14159265358979323846 / 180.0; }

double yield_of(const PlaneStrainStress& s, double c, double phi) {
    return mc_yield(principal_stresses(s), c, phi);
}

// Return mapping: the defining property is that an inadmissible trial stress is
// projected exactly onto the yield surface (f = 0), with the correct active
// region (face / edge / apex) and preserved principal directions (coaxiality).
void test_return_mapping() {
    const double c = 10.0, phi = deg(30.0), psi = deg(10.0);
    const MohrCoulombParams p{20000.0, 0.3, c, phi, psi};
    const double surf_tol = 1.0e-7;

    // (a) Admissible trial is left untouched (no plastic correction).
    {
        PlaneStrainStress el;
        el.in_plane = Eigen::Vector3d(-5.0, -8.0, 1.0);
        el.zz = -6.0;
        const McReturn r = mc_return_mapping(el, p);
        check(!r.plastic, "admissible trial stays elastic");
        check(close(r.stress.in_plane(0), el.in_plane(0)) &&
                  close(r.stress.zz, el.zz),
              "elastic step leaves stress unchanged");
    }
    // (b) Smooth face: three distinct principals, f returns to 0.
    {
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(50.0, -50.0, 0.0);
        s.zz = 0.0;
        const McReturn r = mc_return_mapping(s, p);
        const PrincipalStresses pr = principal_stresses(r.stress);
        check(r.plastic, "smooth: plastic correction applied");
        check(std::fabs(yield_of(r.stress, c, phi)) < surf_tol,
              "smooth: returned stress on yield surface (f=0)");
        check(pr.s1 > pr.s2 && pr.s2 > pr.s3, "smooth: distinct principals");
    }
    // (c) Right edge (triaxial compression): sigma_2 = sigma_3 after return.
    {
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(100.0, -20.0, 0.0);
        s.zz = -20.0;
        const McReturn r = mc_return_mapping(s, p);
        const PrincipalStresses pr = principal_stresses(r.stress);
        check(std::fabs(yield_of(r.stress, c, phi)) < surf_tol,
              "right edge: f=0");
        check(close(pr.s2, pr.s3, 1.0e-6), "right edge: sigma_2 = sigma_3");
    }
    // (d) Left edge (triaxial extension): sigma_1 = sigma_2 after return.
    {
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(20.0, -100.0, 0.0);
        s.zz = 20.0;
        const McReturn r = mc_return_mapping(s, p);
        const PrincipalStresses pr = principal_stresses(r.stress);
        check(std::fabs(yield_of(r.stress, c, phi)) < surf_tol,
              "left edge: f=0");
        check(close(pr.s1, pr.s2, 1.0e-6), "left edge: sigma_1 = sigma_2");
    }
    // (e) Apex: hydrostatic tension beyond the tip returns to p = c*cot(phi).
    {
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(50.0, 50.0, 0.0);
        s.zz = 50.0;
        const McReturn r = mc_return_mapping(s, p);
        const PrincipalStresses pr = principal_stresses(r.stress);
        const double apex = c * std::cos(phi) / std::sin(phi);  // c*cot(phi)
        check(close(pr.s1, apex, 1.0e-6) && close(pr.s3, apex, 1.0e-6),
              "apex: all principals = c*cot(phi)");
    }
    // (f) Coaxiality: a sheared trial returns along the same principal frame.
    {
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(30.0, -30.0, 40.0);
        s.zz = 0.0;
        const PrincipalStresses tr = principal_stresses(s);
        const McReturn r = mc_return_mapping(s, p);
        const PrincipalStresses re = principal_stresses(r.stress);
        check(std::fabs(yield_of(r.stress, c, phi)) < surf_tol,
              "sheared: f=0");
        check(close(re.cos2t, tr.cos2t, 1.0e-6) &&
                  close(re.sin2t, tr.sin2t, 1.0e-6),
              "return is coaxial with the trial stress");
    }
    // (g) Non-dilatant (psi = 0) path still returns admissibly.
    {
        const MohrCoulombParams p0{20000.0, 0.3, c, phi, 0.0};
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(60.0, -40.0, 10.0);
        s.zz = -10.0;
        const McReturn r = mc_return_mapping(s, p0);
        check(std::fabs(yield_of(r.stress, c, phi)) < surf_tol,
              "psi=0: returned stress on yield surface");
    }
    // (h) Tresca limit (phi = psi = 0): the cone degenerates to a prism with no
    //     apex; the criterion reduces to sigma_1 - sigma_3 = 2c and hydrostatic
    //     states are always admissible. This is the constitutive basis of the
    //     Prandtl Nc = 2 + pi bearing-capacity benchmark (P1.3).
    {
        const MohrCoulombParams pt{20000.0, 0.3, c, 0.0, 0.0};
        PlaneStrainStress s;
        s.in_plane = Eigen::Vector3d(30.0, -30.0, 0.0);
        s.zz = 0.0;
        const McReturn r = mc_return_mapping(s, pt);
        const PrincipalStresses pr = principal_stresses(r.stress);
        check(std::fabs((pr.s1 - pr.s3) - 2.0 * c) < 1e-6,
              "phi=0: returns to Tresca surface (sigma_1 - sigma_3 = 2c)");
        PlaneStrainStress h;
        h.in_plane = Eigen::Vector3d(50.0, 50.0, 0.0);
        h.zz = 50.0;
        check(!mc_return_mapping(h, pt).plastic,
              "phi=0: hydrostatic state is admissible");
    }
}

// Consistent (algorithmic) tangent. Because the return map is piecewise linear,
// within an active region the tangent must be the exact local Jacobian of the
// (already verified) stress update. We also check its symmetry structure:
// associated flow (psi = phi) yields a symmetric tangent, non-associated does not.
void test_consistent_tangent() {
    const double c = 10.0, phi = deg(30.0);
    GaussState committed;  // virgin state (zero stress)

    // (a) An elastic step returns the elastic operator exactly.
    {
        const MaterialModel m{MaterialType::MohrCoulomb, 20000.0, 0.3,
                              c, phi, deg(10.0)};
        const Eigen::Matrix3d De = m.elastic_plane_strain();
        GaussState trial;
        Eigen::Matrix3d T;
        integrate_point(m, committed, Eigen::Vector3d(1e-6, -1e-7, 0.0), trial, T);
        check((T - De).cwiseAbs().maxCoeff() < 1e-6 * De.cwiseAbs().maxCoeff(),
              "elastic step: tangent = D_e");
    }
    // (b) Plastic step: D_T reproduces the update for a within-region perturbation,
    //     i.e. sigma(eps+delta) - sigma(eps) == D_T * delta (exact, piecewise linear).
    {
        const MaterialModel m{MaterialType::MohrCoulomb, 20000.0, 0.3,
                              c, phi, deg(10.0)};
        const Eigen::Vector3d de(2e-3, -2e-3, 1e-3);  // yields (smooth region)
        GaussState trial;
        Eigen::Matrix3d T;
        integrate_point(m, committed, de, trial, T);

        const Eigen::Vector3d delta(3e-7, -1e-7, 2e-7);
        GaussState trial_p;
        Eigen::Matrix3d Tp;
        integrate_point(m, committed, de + delta, trial_p, Tp);
        const Eigen::Vector3d predicted = T * delta;
        const Eigen::Vector3d actual = trial_p.stress - trial.stress;
        const double scale = trial.stress.cwiseAbs().maxCoeff() + 1.0;
        check((predicted - actual).cwiseAbs().maxCoeff() < 1e-6 * scale,
              "plastic step: D_T is the consistent (local Jacobian) tangent");
    }
    // (c) Symmetry structure of the elastoplastic tangent.
    {
        auto tangent_for = [&](double psi) {
            const MaterialModel m{MaterialType::MohrCoulomb, 20000.0, 0.3,
                                  c, phi, psi};
            GaussState trial;
            Eigen::Matrix3d T;
            integrate_point(m, committed, Eigen::Vector3d(2e-3, -2e-3, 5e-4),
                            trial, T);
            return T;
        };
        const Eigen::Matrix3d Ta = tangent_for(phi);        // associated
        const Eigen::Matrix3d Tn = tangent_for(deg(5.0));   // non-associated
        const double asym_a =
            (Ta - Ta.transpose()).cwiseAbs().maxCoeff() / Ta.cwiseAbs().maxCoeff();
        const double asym_n =
            (Tn - Tn.transpose()).cwiseAbs().maxCoeff() / Tn.cwiseAbs().maxCoeff();
        check(asym_a < 1e-4, "associated flow (psi=phi): symmetric tangent");
        check(asym_n > 1e-3, "non-associated flow: asymmetric tangent");
    }
}

// Direct validation of the closed-form consistent tangent against a central finite
// difference of the stress update, D_fd(:,j) = [sigma(eps+h e_j) - sigma(eps-h e_j)]
// / 2h. The points carry shear (so the in-plane frame rotates -> exercises the spin
// term) and sit interior to each plastic region (so the FD stencil does not straddle
// a region boundary). This is the independent check that the analytic D_T -- with
// its sigma_zz coupling and per-region Jacobian -- equals the true derivative.
void test_tangent_vs_finite_difference() {
    const double c = 10.0, phi = deg(30.0);
    GaussState committed;

    auto fd_tangent = [&](const MaterialModel& m, const Eigen::Vector3d& de) {
        constexpr double h = 1e-6;
        Eigen::Matrix3d fd;
        for (int j = 0; j < 3; ++j) {
            Eigen::Vector3d ep = de, em = de;
            ep(j) += h;
            em(j) -= h;
            GaussState tp, tm;
            Eigen::Matrix3d dummy;
            integrate_point(m, committed, ep, tp, dummy);
            integrate_point(m, committed, em, tm, dummy);
            fd.col(j) = (tp.stress - tm.stress) / (2.0 * h);
        }
        return fd;
    };

    struct Case { const char* name; Eigen::Vector3d de; };
    // Strains chosen via the elastic predictor (from a virgin state) to reproduce
    // the trial states verified in test_return_mapping. The smooth case carries
    // shear (rotating frame -> spin term); the edge cases are axis-aligned (the
    // sigma_2 = sigma_3 / sigma_1 = sigma_2 coincidence is between an in-plane and
    // the out-of-plane principal, exercising the sigma_zz coupling and the edge
    // Jacobian). All sit interior to their region, so the FD stencil stays put.
    const Case cases[] = {
        {"smooth (with shear)", Eigen::Vector3d(3e-3, -1e-3, 1e-3)},
        {"right edge (triaxial compression)", Eigen::Vector3d(3.03e-3, -4.77e-3, 0.0)},
        {"left edge (triaxial extension)", Eigen::Vector3d(4.77e-3, -3.03e-3, 0.0)},
    };
    for (double psi : {phi, deg(5.0)}) {  // associated and non-associated
        const MaterialModel m{MaterialType::MohrCoulomb, 20000.0, 0.3, c, phi, psi};
        for (const Case& cs : cases) {
            GaussState trial;
            Eigen::Matrix3d T;
            integrate_point(m, committed, cs.de, trial, T);
            const Eigen::Matrix3d fd = fd_tangent(m, cs.de);
            const double err = (T - fd).cwiseAbs().maxCoeff();
            const double scale = T.cwiseAbs().maxCoeff() + 1.0;
            check(err < 1e-3 * scale,
                  "analytic tangent matches central finite difference");
            if (err >= 1e-3 * scale)
                std::fprintf(stderr, "  [%s, psi=%.1fdeg] rel err = %.2e\n",
                             cs.name, psi / deg(1.0), err / scale);
        }
    }
}

// ===================== Tension cut-off (Rankine cap) ==========================
// Formulation doc sec 7: three associated Rankine planes (PLAXIS MMM Eq 3-11)
// reduce, with sorted principals, to sigma_1 <= sigma_t; the combined MC+cap
// return adds the tension face / TT edge / T apex, the MC-tension line and its
// two vertices -- all affine, all with constant exact region Jacobians.

// Closed forms reachable by hand (doc 7g) + the bit-identical guarantees.
void test_tension_cutoff_closed_forms() {
    const LameConstants lame = lame_from(kE * 2.0, kNu);  // E = 20000, nu = 0.3
    const double lam = lame.lambda, G = lame.mu;
    const double st = 5.0;
    const MohrCoulombParams mp{kE * 2.0, kNu, 10.0, deg(30.0), 0.0, true, st};

    // (a) Uniaxial stretch: tension face. After yield sigma_xx = sigma_t and the
    //     lateral stresses sit at lam*sigma_t/(lam+2G); further stretching leaves
    //     BOTH frozen (associated flow removes exactly the new elastic increment).
    {
        const auto tr = elastic_predictor({}, Eigen::Vector3d(1.0e-3, 0.0, 0.0), lame);
        const McReturn r = mc_return_mapping(tr, mp);
        const double lat = lam * st / (lam + 2.0 * G);
        check(r.plastic && r.tension, "uniaxial stretch: tension face active");
        check(close(r.stress.in_plane(0), st), "uniaxial: sigma_xx = sigma_t");
        check(close(r.stress.in_plane(1), lat), "uniaxial: sigma_yy = lam*st/(lam+2G)");
        check(close(r.stress.zz, lat), "uniaxial: sigma_zz = lam*st/(lam+2G)");
        check(close(r.stress.in_plane(2), 0.0), "uniaxial: shear stays zero");
        const auto tr2 =
            elastic_predictor(r.stress, Eigen::Vector3d(1.0e-3, 0.0, 0.0), lame);
        const McReturn r2 = mc_return_mapping(tr2, mp);
        check(close(r2.stress.in_plane(0), st) && close(r2.stress.in_plane(1), lat) &&
                  close(r2.stress.zz, lat),
              "uniaxial continuation: plateau frozen at the closed form");
    }
    // (b) Hydrostatic tensile trial -> tension apex (st, st, st) exactly (the MC
    //     apex c*cot(phi) = 17.3 lies OUTSIDE the cap, so the old apex return
    //     would have been silently wrong here).
    {
        PlaneStrainStress tr;
        tr.in_plane << 20.0, 20.0, 0.0;
        tr.zz = 20.0;
        const McReturn r = mc_return_mapping(tr, mp);
        check(r.tension, "hydrostatic: tension region engaged");
        check(close(r.stress.in_plane(0), st) && close(r.stress.in_plane(1), st) &&
                  close(r.stress.zz, st) && close(r.stress.in_plane(2), 0.0),
              "hydrostatic trial returns to the tension apex (st, st, st)");
    }
    // (c) Equal-biaxial stretch: tension-tension edge (also the isotropic in-plane
    //     trial, radius = 0). r3 follows the TT-edge closed form.
    {
        const auto tr =
            elastic_predictor({}, Eigen::Vector3d(1.0e-3, 1.0e-3, 0.0), lame);
        const McReturn r = mc_return_mapping(tr, mp);
        const double S =
            (tr.in_plane(0) + tr.in_plane(1) - 2.0 * st) / (2.0 * lam + 2.0 * G);
        check(r.tension, "biaxial: tension region engaged");
        check(close(r.stress.in_plane(0), st) && close(r.stress.in_plane(1), st),
              "biaxial stretch: both in-plane stresses capped at sigma_t");
        check(close(r.stress.zz, tr.zz - lam * S),
              "biaxial: sigma_zz per the TT-edge closed form");
    }
    // (d) Cap inactive: a compressive MC-violating trial is BIT-IDENTICAL with the
    //     cut-off on or off (the pure-MC path is untouched).
    {
        PlaneStrainStress tr;
        tr.in_plane << -80.0, -260.0, 30.0;
        tr.zz = -120.0;
        MohrCoulombParams off = mp;
        off.tension_cutoff = false;
        const McReturn a = mc_return_mapping(tr, mp);
        const McReturn b = mc_return_mapping(tr, off);
        check(a.plastic && !a.tension, "compressive: MC plastic, cap inactive");
        check((a.stress.in_plane - b.stress.in_plane).cwiseAbs().maxCoeff() == 0.0 &&
                  a.stress.zz == b.stress.zz &&
                  (a.tangent - b.tangent).cwiseAbs().maxCoeff() == 0.0,
              "compressive: cut-off on == off, bit-identical incl. tangent");
    }
    // (e) Cut-off off: a tension-only trial passes through untouched (the legacy
    //     behaviour); on: the same trial is clipped to sigma_t.
    {
        PlaneStrainStress tr;
        tr.in_plane << 12.0, 2.0, 0.0;
        tr.zz = 6.0;
        MohrCoulombParams off = mp;
        off.tension_cutoff = false;
        const McReturn b = mc_return_mapping(tr, off);
        check(!b.plastic &&
                  (b.stress.in_plane - tr.in_plane).cwiseAbs().maxCoeff() == 0.0,
              "cut-off off: tension carried unclipped (legacy path)");
        const McReturn a = mc_return_mapping(tr, mp);
        check(a.tension && close(principal_stresses(a.stress).s1, st),
              "cut-off on: same trial clipped to sigma_t");
    }
}

// Region inversion: place an admissible target ON each tension region, push the
// trial out along a positive combination of the region's elastic flow directions
// D_e m (its return cone), and require the mapping to land back on the target
// exactly. Verifies the active-set selection AND the closed-form multipliers at
// once. The 4x4 algorithmic Jacobian is then validated against a central finite
// difference of the full map sigma_tr -> sigma_ret (trials are interior to their
// cones, so the stencil stays in-region).
void test_tension_region_inversion() {
    const double E = 20000.0, nu = 0.3, c = 10.0, phi = deg(30.0), psi = deg(10.0);
    const LameConstants lame = lame_from(E, nu);
    const double lam = lame.lambda, G = lame.mu;
    const double st = 5.0;
    const double sphi = std::sin(phi), spsi = std::sin(psi);
    const double sm = ((1.0 + sphi) * st - 2.0 * c * std::cos(phi)) / (1.0 - sphi);
    const MohrCoulombParams mp{E, nu, c, phi, psi, true, st};

    auto De = [&](const Eigen::Vector3d& v) {
        const double tr = v.sum();
        return Eigen::Vector3d(lam * tr + 2.0 * G * v(0), lam * tr + 2.0 * G * v(1),
                               lam * tr + 2.0 * G * v(2));
    };
    const Eigen::Vector3d e1(1, 0, 0), e2(0, 1, 0), e3(0, 0, 1);
    const Eigen::Vector3d m13(1.0 + spsi, 0.0, -(1.0 - spsi));   // MC smooth flow
    const Eigen::Vector3d m12(1.0 + spsi, -(1.0 - spsi), 0.0);   // MC (1,2) face
    const Eigen::Vector3d m23(0.0, 1.0 + spsi, -(1.0 - spsi));   // MC (2,3) face

    struct Case {
        const char* name;
        Eigen::Vector3d target;  // sorted principals on the region
        Eigen::Vector3d push;    // trial = target + push, push in the return cone
    };
    const double mid = 0.5 * (st + sm);
    const Case cases[] = {
        {"tension face", Eigen::Vector3d(st, 0.0, -10.0), De(e1) * 1.0e-3},
        {"tension-tension edge", Eigen::Vector3d(st, st, -10.0),
         De(e1) * 1.0e-3 + De(e2) * 6.0e-4},
        {"tension apex", Eigen::Vector3d(st, st, st),
         De(e1) * 1.2e-3 + De(e2) * 8.0e-4 + De(e3) * 5.0e-4},
        {"MC-tension line", Eigen::Vector3d(st, mid, sm),
         De(m13) * 6.0e-4 + De(e1) * 4.0e-4},
        {"vertex V_L", Eigen::Vector3d(st, st, sm),
         De(e1) * 5.0e-4 + De(e2) * 3.0e-4 + De(m13) * 4.0e-4 + De(m23) * 2.0e-4},
        {"vertex V_R", Eigen::Vector3d(st, sm, sm),
         De(e1) * 4.0e-4 + De(m13) * 5.0e-4 + De(m12) * 3.0e-4},
    };
    for (const Case& cs : cases) {
        const Eigen::Vector3d t = cs.target + cs.push;
        // Arrange (t1, t3) in-plane WITH shear (rotated frame) and t2 out-of-plane:
        // the spectral bookkeeping and the coaxial reconstruction are exercised.
        const double c2 = std::cos(0.3), s2r = std::sin(0.3);
        PlaneStrainStress tr;
        const double mean = 0.5 * (t(0) + t(2)), rad = 0.5 * (t(0) - t(2));
        tr.in_plane << mean + rad * c2, mean - rad * c2, rad * s2r;
        tr.zz = t(1);
        const McReturn ret = mc_return_mapping(tr, mp);
        const PrincipalStresses p = principal_stresses(ret.stress);
        const double tol = 1.0e-8 * (t.cwiseAbs().maxCoeff() + 1.0);
        check(ret.plastic && ret.tension, cs.name);
        const bool hit = std::fabs(p.s1 - cs.target(0)) < tol &&
                         std::fabs(p.s2 - cs.target(1)) < tol &&
                         std::fabs(p.s3 - cs.target(2)) < tol;
        check(hit, "region inversion lands on the constructed target");
        if (!hit)
            std::fprintf(stderr, "  [%s] got (%.6f, %.6f, %.6f) want (%.6f, %.6f, %.6f)\n",
                         cs.name, p.s1, p.s2, p.s3, cs.target(0), cs.target(1),
                         cs.target(2));

        Eigen::Matrix4d fd;  // central FD of sigma_tr -> sigma_ret, Voigt [xx,yy,xy,zz]
        const double h = 1.0e-5;
        for (int j = 0; j < 4; ++j) {
            PlaneStrainStress tp = tr, tm = tr;
            if (j < 3) {
                tp.in_plane(j) += h;
                tm.in_plane(j) -= h;
            } else {
                tp.zz += h;
                tm.zz -= h;
            }
            const McReturn rp = mc_return_mapping(tp, mp);
            const McReturn rm = mc_return_mapping(tm, mp);
            for (int i = 0; i < 3; ++i)
                fd(i, j) = (rp.stress.in_plane(i) - rm.stress.in_plane(i)) / (2.0 * h);
            fd(3, j) = (rp.stress.zz - rm.stress.zz) / (2.0 * h);
        }
        const double jerr = (ret.algo_jacobian - fd).cwiseAbs().maxCoeff();
        check(jerr < 1.0e-4, "algorithmic Jacobian matches the central FD in-region");
        if (jerr >= 1.0e-4)
            std::fprintf(stderr, "  [%s] J err = %.3e\n", cs.name, jerr);
    }
}

// Deterministic pseudo-random sweep: whatever the trial, the returned stress must
// satisfy BOTH surfaces and be a fixed point of the mapping (idempotency). With
// the cut-off off the same population must carry s1 > sigma_t somewhere -- the
// witness that the cap actually binds and the sweep covers tension states.
void test_tension_feasibility_sweep() {
    const double E = 20000.0, nu = 0.3, c = 10.0, phi = deg(30.0);
    unsigned long long lcg = 0x243F6A8885A308D3ull;  // fixed seed: reproducible
    auto urand = [&]() {
        lcg = lcg * 6364136223846793005ull + 1442695040888963407ull;
        return (double)(lcg >> 11) / 9007199254740992.0;
    };
    int bad_cap = 0, bad_mc = 0, bad_idem = 0, witness = 0, tension_hits = 0;
    for (double st : {0.0, 5.0}) {
        for (double psi : {0.0, deg(10.0)}) {
            const MohrCoulombParams on{E, nu, c, phi, psi, true, st};
            const MohrCoulombParams off{E, nu, c, phi, psi, false, st};
            for (int k = 0; k < 2000; ++k) {
                PlaneStrainStress tr;
                tr.in_plane << 300.0 * urand() - 200.0, 300.0 * urand() - 200.0,
                    120.0 * urand() - 60.0;
                tr.zz = 300.0 * urand() - 200.0;
                const McReturn r = mc_return_mapping(tr, on);
                const PrincipalStresses p = principal_stresses(r.stress);
                const double tol =
                    1.0e-6 * (std::fabs(p.s1) + std::fabs(p.s3) + c + st + 1.0);
                if (p.s1 > st + tol) ++bad_cap;
                if (mc_yield(p, c, phi) > tol) ++bad_mc;
                if (r.tension) ++tension_hits;
                const McReturn r2 = mc_return_mapping(r.stress, on);
                if ((r2.stress.in_plane - r.stress.in_plane).cwiseAbs().maxCoeff() +
                        std::fabs(r2.stress.zz - r.stress.zz) >
                    tol)
                    ++bad_idem;
                if (principal_stresses(mc_return_mapping(tr, off).stress).s1 >
                    st + tol)
                    ++witness;
            }
        }
    }
    check(bad_cap == 0, "sweep(8000): returned s1 <= sigma_t everywhere");
    check(bad_mc == 0, "sweep(8000): returned f_mc <= 0 everywhere");
    check(bad_idem == 0, "sweep(8000): the mapping is idempotent");
    check(tension_hits > 500, "sweep: tension regions genuinely exercised");
    check(witness > 500, "sweep: cut-off off carries tension (the cap binds)");
    if (bad_cap + bad_mc + bad_idem)
        std::fprintf(stderr, "  bad_cap=%d bad_mc=%d bad_idem=%d\n", bad_cap, bad_mc,
                     bad_idem);
}

// Through integrate_point (both kinematics): the MaterialModel fields reach the
// return mapping, the tension-region tangent is the local Jacobian of the update,
// and the axisymmetric branch caps the hoop stress as well.
void test_tension_integrate_point() {
    MaterialModel m;
    m.type = MaterialType::MohrCoulomb;
    m.youngs_modulus = 20000.0;
    m.poisson_ratio = 0.3;
    m.cohesion = 10.0;
    m.friction_angle = deg(30.0);
    m.dilatancy_angle = 0.0;
    m.tension_cutoff = true;
    m.tensile_strength = 5.0;
    const LameConstants lame = lame_from(m.youngs_modulus, m.poisson_ratio);
    GaussState committed;

    // (a) plane strain: uniaxial stretch -> tension face; D_T is the local Jacobian.
    {
        GaussState trial;
        Eigen::Matrix3d T;
        const Eigen::Vector3d de(1.0e-3, 0.0, 0.0);
        integrate_point(m, committed, de, trial, T);
        const double lat = lame.lambda * 5.0 / (lame.lambda + 2.0 * lame.mu);
        check(close(trial.stress(0), 5.0), "integrate_point: sigma_xx capped");
        check(close(trial.stress(1), lat) && close(trial.stress_zz, lat),
              "integrate_point: face closed-form laterals");
        const Eigen::Vector3d delta(3e-7, -1e-7, 2e-7);
        GaussState tp;
        Eigen::Matrix3d Tp;
        integrate_point(m, committed, de + delta, tp, Tp);
        const Eigen::Vector3d pred = T * delta, act = tp.stress - trial.stress;
        check((pred - act).cwiseAbs().maxCoeff() <
                  1e-6 * (trial.stress.cwiseAbs().maxCoeff() + 1.0),
              "integrate_point: tension-region tangent is consistent");
    }
    // (b) axisymmetric: hydrostatic stretch -> tension apex in r, z AND hoop.
    {
        GaussState trial;
        Eigen::Matrix4d T;
        Eigen::Vector4d de;
        de << 1.0e-3, 1.0e-3, 0.0, 1.0e-3;
        katai::core::integrate_point_axisym(m, committed, de, trial, T);
        check(close(trial.stress(0), 5.0) && close(trial.stress(1), 5.0) &&
                  close(trial.stress_zz, 5.0) && close(trial.stress(2), 0.0),
              "axisym: tension apex caps r, z and hoop stresses");
    }
}

} // namespace

int main() {
    test_elastic_predictor();
    test_principal_stresses();
    test_yield_function();
    test_return_mapping();
    test_consistent_tangent();
    test_tangent_vs_finite_difference();
    test_tension_cutoff_closed_forms();
    test_tension_region_inversion();
    test_tension_feasibility_sweep();
    test_tension_integrate_point();
    if (g_failures == 0) {
        std::printf("OK: Mohr-Coulomb kernels (P1.2a) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
