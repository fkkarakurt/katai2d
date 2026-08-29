// The Hoek-Brown criterion at the material point (PLAXIS 2D Material Models Manual §4; the 2002
// edition of Hoek, Carranza-Torres & Corkum). This tree had no rock model at all: rock had to be
// entered as a Mohr-Coulomb fit, which is a straight line drawn through a curve -- it can be made
// to match over a narrow band of confining stress and is wrong outside it in both directions.
//
// FOUR CLOSED FORMS AND ONE IDENTITY. The first three come from the manual's own equations and
// check that the criterion was transcribed rather than remembered; the fourth checks the return
// mapping lands ON the surface; and the fifth is the one that would catch the error nothing else
// would -- the ordering.
//
//  (1) the uni-axial compressive strength of the rock mass, sigma_c = -|sigma_ci| s^a (Eq 4-5):
//      the state (s1 = 0, s3 = sigma_c) must be exactly on the surface;
//  (2) the tensile strength, sigma_t = s |sigma_ci| / m_b (Eq 4-6): the isotropic tension point
//      must be exactly on it too -- it is the apex, where the criterion's own bracket closes;
//  (3) the envelope itself (Eq 4-1) at several confining stresses;
//  (4) a return from an inadmissible trial must satisfy f = 0 to round-off;
//  (5) THE ORDERING. The manual writes the criterion with sigma'_1 <= sigma'_2 <= sigma'_3, its
//      sigma'_1 being the most COMPRESSIVE; this solver sorts the other way. A model that got
//      that backwards would run, converge, and be wrong in the direction that matters -- so the
//      envelope is checked against the manual's own numbers rather than against itself.
//
// AND THE IDENTITY. With m_b -> 0 and s = 1 the bracket becomes constant and the criterion
// degenerates to s3 - s1 + |sigma_ci| = 0, which is Tresca: a friction-free material of cohesion
// |sigma_ci| / 2. That case this tree already solves, by a completely different return mapping
// that has been verified since P1.2b. The two must agree to round-off on any trial state, and if
// they do, the whole apparatus -- decomposition, flow direction, corrector, ordering -- is right
// on at least one line through it.
//
// verify: KV-CST-014
//   oracle:   closed_form
//   source:   PLAXIS 2D Material Models Manual (2025.1) §4, the generalised Hoek-Brown criterion of Hoek, Carranza-Torres & Corkum (2002): yield Eq 4-1 and 4-8, the rock-mass constants m_b (Eq 4-2), s (Eq 4-3) and a (Eq 4-4), the rock-mass compressive strength Eq 4-5 and tensile strength Eq 4-6, the non-associated flow of Eq 4-10/4-12 and the mobilised dilatancy of Eq 4-13/4-14. The Tresca degeneration is checked against this tree's own verified Mohr-Coulomb return (mohr_coulomb.cpp, P1.2b)
//   locator:  f = s3 - s1 + |sigma_ci| (m_b (-s1/|sigma_ci|) + s)^a with s1 >= s2 >= s3 tension-positive (the manual's sigma'_3 and sigma'_1 respectively); m_b = m_i exp((GSI-100)/(28-14D)); s = exp((GSI-100)/(9-3D)); a = 1/2 + (exp(-GSI/15) - exp(-20/3))/6; sigma_c = -|sigma_ci| s^a; sigma_t = s |sigma_ci| / m_b
//   quantity: the yield function at the closed-form uni-axial and isotropic-tensile states [kPa]; the failure sigma_1 at four confining stresses against Eq 4-1 [kPa]; the yield function after a return from an inadmissible trial [kPa]; and the returned principal stresses of the Tresca degeneration against this tree's Mohr-Coulomb return [kPa]; and the tensile limit the criterion admits with and without the optional cut-off of sec 4.3.7 [kPa]
//   expected: zero at the two closed-form states, the Eq 4-1 envelope at the four confining stresses, zero after every return, and the Mohr-Coulomb answer in the degeneration
//   band:     1e-9 relative on the constants, the two strengths and the envelope -- they are the same arithmetic evaluated two ways, so anything larger would be a transcription error rather than a discretisation, and all of them came out at 0.000e+00. 1e-6 kPa absolute on the yield function after a return (nine trials, all on the surface and still ordered) and on the apex. For the Tresca degeneration the band is on the LIMIT rather than on a number: at m_b = 1e-6 the two returns differ by 1.803e-04 kPa and at m_b = 1e-10 by 1.803e-08, a factor of 1e4 for a factor of 1e4 -- the difference is the criterion linearising towards Tresca, not the corrector, and a floor there would have meant the opposite. The edge return is what made that possible: clamping the intermediate stress instead left 2.5e+01 kPa on a material of cohesion 50. The tension cut-off is checked on the ADMISSIBILITY of an isotropic tensile state at 0.99 and 1.01 times the cap rather than on the constant, because a cap stored and never read would pass the constant; the uncapped rock admits both, which is what says the cut-off is doing the binding
#include <katai/materials/hoek_brown.hpp>
#include <katai/materials/mohr_coulomb.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::MohrCoulombParams;
using katai::core::PlaneStrainStress;
namespace hb = katai::core::hoekbrown;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// A weak, disturbed rock mass: 50 MPa intact rock, GSI 50, D 0.5 -- the kind of mass a tunnel
// through weathered ground is driven in, and far enough from intact rock that m_b, s and a all
// differ visibly from their intact values.
hb::Params rock() {
    hb::Params P;
    P.E = 3.0e6; P.nu = 0.25;
    P.sigci = 50000.0;   // kPa
    P.mi = 10.0; P.gsi = 50.0; P.D = 0.5;
    P.psi = 0.0; P.sig_psi = 0.0;
    return P;
}

void test_constants_and_closed_forms() {
    std::printf("-- (1) the rock-mass constants and the two strengths the manual gives in closed form --\n");
    const hb::Params P = rock();
    const hb::Constants C = hb::constants_of(P);
    // Eq 4-2/4-3/4-4, evaluated here independently of the header's own arithmetic.
    const double mb_ref = 10.0 * std::exp((50.0 - 100.0) / (28.0 - 14.0 * 0.5));
    const double s_ref = std::exp((50.0 - 100.0) / (9.0 - 3.0 * 0.5));
    const double a_ref = 0.5 + (std::exp(-50.0 / 15.0) - std::exp(-20.0 / 3.0)) / 6.0;
    std::printf("   m_b = %.6f (ref %.6f), s = %.6e (ref %.6e), a = %.6f (ref %.6f)\n",
                C.mb, mb_ref, C.s, s_ref, C.a, a_ref);
    check(std::fabs(C.mb / mb_ref - 1.0) < 1e-9, "m_b = m_i exp((GSI-100)/(28-14D))");
    check(std::fabs(C.s / s_ref - 1.0) < 1e-9, "s = exp((GSI-100)/(9-3D))");
    check(std::fabs(C.a / a_ref - 1.0) < 1e-9, "a = 1/2 + (exp(-GSI/15) - exp(-20/3))/6");

    const double sigc_ref = -50000.0 * std::pow(s_ref, a_ref);
    const double sigt_ref = s_ref * 50000.0 / mb_ref;
    std::printf("   rock-mass sigma_c = %.4f kPa (ref %.4f),  sigma_t = %.6f kPa (ref %.6f)\n",
                C.sigc, sigc_ref, C.sigt, sigt_ref);
    check(std::fabs(C.sigc / sigc_ref - 1.0) < 1e-9, "sigma_c = -|sigma_ci| s^a  (Eq 4-5)");
    check(std::fabs(C.sigt / sigt_ref - 1.0) < 1e-9, "sigma_t = s |sigma_ci| / m_b  (Eq 4-6)");

    // The two states the closed forms name must lie exactly ON the surface.
    const double f_uni = hb::yield(0.0, C.sigc, P, C);
    const double f_apex = hb::yield(C.sigt, C.sigt, P, C);
    std::printf("   f at the uni-axial state = %.3e kPa;  f at the isotropic tensile apex = %.3e kPa\n",
                f_uni, f_apex);
    check(std::fabs(f_uni) < 1e-9 * std::fabs(C.sigc), "uni-axial compression sits on the surface");
    check(std::fabs(f_apex) < 1e-9 * std::fabs(P.sigci), "and so does the tensile apex");
}

void test_envelope() {
    std::printf("\n-- (2) the envelope itself, against Eq 4-1 at four confining stresses --\n");
    const hb::Params P = rock();
    const hb::Constants C = hb::constants_of(P);
    std::printf("   sigma3 [kPa]   sigma1 (Eq 4-1)      f there\n");
    int n = 0;
    for (double s3conf : {-500.0, -2000.0, -8000.0, -20000.0}) {
        // The manual's Eq 4-1 with ITS ordering: sigma'_1 (most compressive) from sigma'_3 (least).
        // In this file's ordering the confining stress is s1 and the failure stress is s3.
        const double s1 = s3conf;
        const double bracket = C.mb * (-s1 / std::fabs(P.sigci)) + C.s;
        const double s3fail = s1 - std::fabs(P.sigci) * std::pow(bracket, C.a);
        const double f = hb::yield(s1, s3fail, P, C);
        std::printf("   %10.1f   %16.4f   %10.3e\n", s1, s3fail, f);
        check(std::fabs(f) < 1e-9 * std::fabs(s3fail), "the envelope point is on the surface");
        // And a point inside it is admissible, a point beyond it is not.
        // s3fail is compressive (negative), so 0.9x is LESS loaded and 1.1x is more. This pair is
        // what caught the sign of the yield function: the envelope point alone reads f = 0 with
        // either sign, and only a question that has a side can tell them apart.
        check(hb::yield(s1, 0.9 * s3fail, P, C) < 0.0, "a less loaded state is admissible");
        check(hb::yield(s1, 1.1 * s3fail, P, C) > 0.0, "a more loaded one is not");
        ++n;
    }
    check(n == 4, "checked four confining stresses");
}

void test_return_lands_on_the_surface() {
    std::printf("\n-- (3) a return from an inadmissible trial lands ON the surface --\n");
    const hb::Params P = rock();
    const hb::Constants C = hb::constants_of(P);
    int n = 0;
    for (double over : {1.2, 2.0, 5.0}) {
        for (double s1 : {-1000.0, -5000.0, 0.0}) {
            const double bracket = C.mb * (-s1 / std::fabs(P.sigci)) + C.s;
            const double s3fail = s1 - std::fabs(P.sigci) * std::pow(bracket, C.a);
            const double s3 = over * s3fail;                 // beyond the surface
            const double s2 = 0.5 * (s1 + s3);
            const hb::Return R = hb::return_mapping(s1, s2, s3, P, C);
            const double f = hb::yield(R.s1, R.s3, P, C);
            check(R.plastic, "the trial was recognised as inadmissible");
            check(std::fabs(f) < 1e-6, "the returned state satisfies f = 0");
            check(R.s1 >= R.s2 - 1e-9 && R.s2 >= R.s3 - 1e-9,
                  "and the returned principals are still ordered");
            ++n;
        }
    }
    std::printf("   %d returns, all on the surface and ordered\n", n);
    // A state pulled into tension beyond the apex returns TO the apex, which is the isotropic
    // tension point -- the criterion has no surface past it.
    const hb::Return A = hb::return_mapping(2.0 * C.sigt, 2.0 * C.sigt, 2.0 * C.sigt, P, C);
    std::printf("   pulled past the apex: returned to (%.6f, %.6f, %.6f), sigma_t = %.6f\n",
                A.s1, A.s2, A.s3, C.sigt);
    check(A.apex && std::fabs(A.s1 - C.sigt) < 1e-9 * std::fabs(C.sigt) &&
              std::fabs(A.s3 - C.sigt) < 1e-9 * std::fabs(C.sigt),
          "a state pulled past the apex returns to the isotropic tension point");
}

void test_tresca_identity() {
    std::printf("\n-- (4) the identity: with m_b -> 0 and s = 1 this IS Tresca, which is already verified --\n");
    // The degeneration is a LIMIT, not an equality: m_b = 0 exactly would divide by zero in the
    // transformed stresses, so the test approaches it and shows the approach. At m_b = 1e-6 the
    // criterion still differs from Tresca by its own linearisation, O(m_b sigma / sigma_ci); make
    // m_b ten thousand times smaller and the difference must fall by the same factor. If it did
    // not -- if it flattened out at some floor -- that floor would be the return mapping's error
    // and not the approximation's, and this test would be reporting the wrong thing.
    double diff_at[2] = {0.0, 0.0};
    const double mis[2] = {1.0e-6, 1.0e-10};
    for (int m = 0; m < 2; ++m) {
        hb::Params P = rock();
        P.gsi = 100.0; P.D = 0.0;     // s = 1, a = 1/2
        P.mi = mis[m];                // m_b -> 0: the bracket stops depending on the stress
        P.sigci = 100.0;              // so the criterion is s1 - s3 - 100 = 0, i.e. c = 50, phi = 0
        P.psi = 0.0;
        const hb::Constants C = hb::constants_of(P);

        MohrCoulombParams mc;
        mc.youngs_modulus = P.E; mc.poisson_ratio = P.nu;
        mc.cohesion = 0.5 * P.sigci; mc.friction_angle = 0.0; mc.dilatancy_angle = 0.0;

        const double trials[][4] = {{-300.0, -100.0, 0.0, -200.0},
                                    {-500.0, -50.0, 40.0, -150.0},
                                    {-80.0, -400.0, -30.0, -240.0},
                                    {50.0, -260.0, 10.0, -100.0}};
        int edges = 0;
        for (const auto& t : trials) {
            PlaneStrainStress st;
            st.in_plane << t[0], t[1], t[2];
            st.zz = t[3];
            const auto pr = katai::core::principal_stresses(st);
            const hb::Return R = hb::return_mapping(pr.s1, pr.s2, pr.s3, P, C);
            const katai::core::McReturn M = katai::core::mc_return_mapping(st, mc);
            const auto mp = katai::core::principal_stresses(M.stress);
            const double d = std::max({std::fabs(R.s1 - mp.s1), std::fabs(R.s2 - mp.s2),
                                       std::fabs(R.s3 - mp.s3)});
            diff_at[m] = std::fmax(diff_at[m], d);
            if (R.edge) ++edges;
            if (m == 0) check(R.plastic == M.plastic, "both models agree on whether the step was plastic");
        }
        std::printf("   m_i = %-8.0e -> worst |difference| vs the Mohr-Coulomb return = %.3e kPa"
                    "   (%d of 4 trials returned on the EDGE)\n", mis[m], diff_at[m], edges);
    }
    const double ratio = diff_at[0] / std::fmax(diff_at[1], 1e-300);
    std::printf("   the difference fell by %.0fx when m_b fell by 1e4 -- so it IS the "
                "approximation, not the return\n", ratio);
    check(diff_at[0] < 1e-3, "at m_b = 1e-6 the two returns agree to the size of the linearisation");
    check(diff_at[1] < 1e-6, "and at 1e-10 they agree to round-off");
    check(ratio > 3.0e3 && ratio < 3.0e4,
          "the difference scales LINEARLY with m_b, which is what says it is the approximation "
          "and not an error in the corrector");
}

void test_mobilised_dilatancy() {
    std::printf("\n-- (5) the mobilised dilatancy follows Eq 4-13/4-14 --\n");
    hb::Params P = rock();
    P.psi = 0.2; P.sig_psi = 4000.0;
    const hb::Constants C = hb::constants_of(P);
    check(std::fabs(hb::psi_mobilised(0.0, P, C) - P.psi) < 1e-12,
          "at sigma'_3 = 0 the mobilised dilatancy is the input value");
    check(std::fabs(hb::psi_mobilised(-P.sig_psi, P, C)) < 1e-12,
          "at -sigma_psi it has died out");
    check(std::fabs(hb::psi_mobilised(-0.5 * P.sig_psi, P, C) - 0.5 * P.psi) < 1e-12,
          "and it decreases linearly between the two");
    check(hb::psi_mobilised(-2.0 * P.sig_psi, P, C) == 0.0, "beyond it, it stays zero");
    check(hb::psi_mobilised(0.5 * C.sigt, P, C) > P.psi,
          "in the tensile range it is raised, so plastic expansion stays possible (Eq 4-14)");
}

}  // namespace

// (6) The optional tension cut-off of sec 4.3.7: it caps the criterion's OWN sigma_t and can only
// lower it. Checked on the thing the cap is for -- the isotropic tensile state the criterion admits
// -- rather than on the constant, because a cap that were stored and not read would pass the latter.
void test_tension_cutoff() {
    std::printf("\n-- (6) the optional tension cut-off caps sigma_t, downwards only (sec 4.3.7) --\n");
    const hb::Params P0 = rock();
    const double sigt0 = hb::constants_of(P0).sigt;
    std::printf("   the criterion's own sigma_t = %.4f kPa\n", sigt0);

    // A cap BELOW it binds: the isotropic tensile state just past the cap must be inadmissible,
    // and the one just inside it admissible. That is what a tensile strength MEANS.
    hb::Params Pc = rock();
    Pc.tension_cutoff = true; Pc.sigt_user = 0.4 * sigt0;
    const hb::Constants Cc = hb::constants_of(Pc);
    check(std::fabs(Cc.sigt - 0.4 * sigt0) < 1e-9 * sigt0,
          "a cut-off below sigma_t becomes the tensile limit");
    const double in = 0.99 * Cc.sigt, out = 1.01 * Cc.sigt;
    check(hb::yield(in, in, Pc, Cc) < 0.0, "an isotropic tension just inside the cap is admissible");
    check(hb::yield(out, out, Pc, Cc) > 0.0, "and just outside it is not");
    // The same two states against the UNCAPPED rock: both admissible, so the cap is what moved
    // the answer and not the states.
    const hb::Constants C0 = hb::constants_of(P0);
    check(hb::yield(out, out, P0, C0) < 0.0,
          "and the uncapped rock admits both, so it is the cut-off that bound");

    // A cap ABOVE sigma_t claims a strength the criterion has not got, and is ignored.
    hb::Params Pa = rock();
    Pa.tension_cutoff = true; Pa.sigt_user = 5.0 * sigt0;
    check(std::fabs(hb::constants_of(Pa).sigt - sigt0) < 1e-12 * sigt0,
          "a cut-off above sigma_t does nothing: the criterion's own limit governs");
    // And off is off.
    hb::Params Pn = rock();
    Pn.tension_cutoff = false; Pn.sigt_user = 0.1 * sigt0;
    check(std::fabs(hb::constants_of(Pn).sigt - sigt0) < 1e-12 * sigt0,
          "with the cut-off off, the entered value is not read at all");
}

int main() {
    std::printf("Hoek-Brown (rock behaviour) at the material point\n\n");
    test_constants_and_closed_forms();
    test_envelope();
    test_return_lands_on_the_surface();
    test_tresca_identity();
    test_mobilised_dilatancy();
    test_tension_cutoff();
    if (g_failures == 0) {
        std::printf("\nOK: the criterion is the manual's, and its degeneration is this tree's own Tresca\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
