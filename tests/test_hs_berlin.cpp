// Hardening Soil -- the Berlin Sand III parameter set in a drained triaxial test. The robust
// two-surface substepping integrator (hs_integrate, cap+shear) drives a drained triaxial at
// sigma3 = 200 kPa with the cap CALIBRATED to (K0^NC, Eoed_ref); we verify the full response:
//   - deviatoric: qf = (c cot phi + sigma3) 2 sin phi/(1-sin phi) = 646 kPa,
//     reached as a perfectly-plastic plateau; secant at qf/2 == E50 (HS hyperbola);
//   - volumetric: initial CONTRACTION (cap, p rising) then net
//     DILATION (shear, psi=6 deg) -- eps_v turns negative and reaches ~ -0.005 at axial 5%.
// The volumetric dilation needs (a) the Rowe dilatancy with the correct sign
// (eps_v^p/eps_q^p = -sin psi_m) and (b) the perfectly-plastic Mohr-Coulomb flow on the qf
// plateau (the hardening shear surface cannot grow past qf -> it must hand over to the MC
// failure surface with continued dilatant flow). (See hardening-soil-formulation.md sec 4-5.)
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

HardeningSoilParams berlin_sand_iii() {
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 105e3; p.Eur_ref = 315e3; p.Eoed_ref = 105e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38 * kPi / 180;
    p.dilatancy = 6 * kPi / 180; p.cohesion = 1.0; p.Rf = 0.9;
    return p;
}

void test_berlin_triaxial() {
    HardeningSoilParams p = berlin_sand_iii();
    // Calibrate the cap to (K0^NC=0.38, Eoed_ref) -- the standard HS workflow; the cap is
    // then ON for the fully-coupled cap+shear triaxial (volumetric response).
    katai::core::hs_calibrate_cap(p, 0.38);
    const double sigma3 = 200.0;
    const double qf = p.q_failure(sigma3);
    const double Ei = p.Ei(sigma3), qa = p.q_asymptote(sigma3), E50 = p.E50(sigma3);
    std::printf("  Berlin Sand III triaxial (sigma3=%.0f): qf=%.1f E50=%.0f "
                "[cap alpha=%.2f beta=%.2e]\n", sigma3, qf, E50, p.cap_alpha, p.cap_beta);
    check(close(qf, 646.0, 0.01), "qf = 646 kPa ((c cot phi + sigma3) 2 sin phi/(1-sin phi))");

    // Drained triaxial: sigma3 const, axial strain driven; pp initialised to the isotropic
    // consolidation stress sigma3 (NC). Solve lateral strain each step for sigma3 = const.
    // Parameterised by the integration tolerance, because the distance from the plateau is a
    // function of it -- see the plateau check below.
    double q50_err = -1.0, epsv_max_contr = 0.0, epsv_final = 0.0, q_final = 0.0;
    const auto walk = [&](double stol) {
        Eigen::Vector3d sig(sigma3, sigma3, sigma3);
        double gp = 0, pp = sigma3, eps1 = 0, epsv = 0;
        q50_err = -1.0; epsv_max_contr = 0.0; epsv_final = 0.0; q_final = 0.0;
        for (int s = 0; s < 12000; ++s) {
            double dlat = 0.0; HsIntegrated r;
            for (int it = 0; it < 30; ++it) {
                r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(5e-6, dlat, dlat), stol);
                const double g = r.stress(2) - sigma3;
                if (std::fabs(g) <= 1e-8 * (1 + sigma3)) break;
                dlat -= g / (r.tangent(2, 1) + r.tangent(2, 2));
            }
            epsv += 5e-6 + 2.0 * dlat;  // axial + 2*lateral (compression +)
            sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += 5e-6;
            epsv_max_contr = std::fmax(epsv_max_contr, epsv);  // peak contraction (positive)
            const double q = sig(0) - sig(2);
            if (q50_err < 0 && q >= 0.5 * qf) q50_err = std::fabs((q / eps1) - E50) / E50;
            q_final = q; epsv_final = epsv;
            if (eps1 > 0.05) break;
        }
        return q_final;
    };
    walk(0.0);   // the shipped integration tolerance
    // Printed to three decimals, not one: what this line is read for is HOW CLOSE the deviator
    // gets to its plateau, and the two returns of the failure bound that this tree once carried
    // were recorded as differing by 0.16% of q_f. One decimal cannot show a difference that size,
    // which is how a recorded difference can quietly stop reproducing -- as that one had.
    std::printf("  reached q=%.3f (qf=%.3f, %.3f%% of it) secant@50%%-vs-E50 err=%.2e | "
                "eps_v: peak-contr=%+.6f  final(@5%%)=%+.6f\n",
                q_final, qf, 100.0 * q_final / qf, q50_err, epsv_max_contr, epsv_final);
    // Deviatoric (shear) agreement -- unchanged from the shear-dominated validation.
    const double q_shipped = q_final;
    check(close(q_shipped, qf, 5e-3), "q reaches the qf plateau (perfect plasticity at failure)");
    check(q50_err >= 0 && q50_err < 0.05, "secant at qf/2 = E50 (HS hyperbola)");
    // Volumetric: initial contraction then net dilation.
    check(epsv_max_contr > 2e-4 && epsv_max_contr < 3e-3,
          "initial volumetric contraction (cap), small peak ~0.1%");
    check(epsv_final < -2.5e-3,
          "net dilation at 5% axial strain (psi=6 deg), eps_v < -0.0025");
    check(epsv_final > -8e-3, "dilation magnitude physical (not runaway)");

    // WHY THE PLATEAU IS NOT REACHED EXACTLY, asserted rather than tolerated. The check above is
    // a 0.5% band and the shortfall at the shipped integration tolerance is 0.31%, so the band
    // cannot see the shortfall at all -- and a band that cannot see the thing it is standing in
    // front of is not a guard. The shortfall has a cause: the Mohr-Coulomb failure bound is a
    // one-sided clamp on sigma1, so it RATCHETS on whatever error the substepping leaves in
    // sigma3, and the size of that error is STOL. It is therefore not a defect in the bound but a
    // declared, controllable quantity, and the falsifiable form of that claim is that TIGHTENING
    // the integration must move the deviator TOWARDS the plateau. Measured with
    // study_hs_integration (5% / 20% axial strain): 1e-3 gives 98.671% / 96.177% -- falling with
    // strain, an artefact that imitates softening -- 1e-5 (shipped) 99.686% / 99.686%, 1e-7
    // 99.949% / 99.984%, 1e-9 99.993% / 100.011%.
    //
    // Until 2026-08-24 this shortfall was recorded as the bound's one-sidedness alone, curable by
    // returning along the consistent direction De.n_s (a seam, KATAI_HS_MCPROJ). Re-measured, that
    // return leaves the deviator at exactly the SAME 99.686% and still costs the strip footing, so
    // it was removed. The direction below is what the record now rests on.
    const double q_tight = walk(1e-8);
    std::printf("  plateau vs integration tolerance: shipped %.3f%% of qf -> stol 1e-8 %.3f%%\n",
                100.0 * q_shipped / qf, 100.0 * q_tight / qf);
    // Absolute distance on both sides, not the signed shortfall: at 1e-8 the deviator lands
    // 0.03% ABOVE q_f (the drift loop's own convergence), and a signed comparison would call an
    // arbitrarily large overshoot an improvement.
    check(std::fabs(qf - q_tight) < std::fabs(qf - q_shipped),
          "tightening the integration moves the deviator TOWARDS the failure plateau, which is "
          "what makes the shortfall integration error rather than a defect in the bound");
    check(close(q_tight, qf, 1e-3),
          "...and it gets within 0.1% of it, so the plateau is reached and not merely approached");
}

} // namespace

int main() {
    test_berlin_triaxial();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil verified on Berlin Sand III (drained triaxial: "
                    "qf, E50, AND volumetric contraction->dilation)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
