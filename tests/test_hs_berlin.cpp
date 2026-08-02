// Hardening Soil -- validation against the documented, PLAXIS-validated Berlin Sand III
// drained triaxial test (Rocscience/PLAXIS 2014, Table 15.1 + Fig 15.4; cross-checked with
// the PLAXIS 2D Material Models Manual sec 6, Fig 15.4). This is the concrete "as accurate
// as PLAXIS" proof. The robust two-surface substepping integrator (hs_integrate, cap+shear)
// drives a drained triaxial at sigma3 = 200 kPa with the cap CALIBRATED to (K0^NC, Eoed_ref);
// we verify the full Fig 15.4 response:
//   - deviatoric: qf = (c cot phi + sigma3) 2 sin phi/(1-sin phi) = 646 kPa (plateau ~640),
//     reached as a perfectly-plastic plateau; secant at qf/2 == E50 (HS hyperbola);
//   - volumetric (Fig 15.4 lower-left): initial CONTRACTION (cap, p rising) then net
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
    // then ON for the fully-coupled cap+shear triaxial (Fig 15.4 volumetric).
    katai::core::hs_calibrate_cap(p, 0.38);
    const double sigma3 = 200.0;
    const double qf = p.q_failure(sigma3);
    const double Ei = p.Ei(sigma3), qa = p.q_asymptote(sigma3), E50 = p.E50(sigma3);
    std::printf("  Berlin Sand III triaxial (sigma3=%.0f): qf=%.1f E50=%.0f "
                "[cap alpha=%.2f beta=%.2e]\n", sigma3, qf, E50, p.cap_alpha, p.cap_beta);
    check(close(qf, 646.0, 0.01), "qf = 646 kPa (matches PLAXIS Fig 15.4 plateau ~640)");

    // Drained triaxial: sigma3 const, axial strain driven; pp initialised to the isotropic
    // consolidation stress sigma3 (NC). Solve lateral strain each step for sigma3 = const.
    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0, pp = sigma3, eps1 = 0, epsv = 0;
    double q50_err = -1.0, epsv_max_contr = 0.0, epsv_final = 0.0, q_final = 0.0;
    for (int s = 0; s < 12000; ++s) {
        double dlat = 0.0; HsIntegrated r;
        for (int it = 0; it < 30; ++it) {
            r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(5e-6, dlat, dlat));
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
    std::printf("  reached q=%.1f (qf=%.1f) secant@50%%-vs-E50 err=%.2e | "
                "eps_v: peak-contr=%+.5f  final(@5%%)=%+.5f\n",
                q_final, qf, q50_err, epsv_max_contr, epsv_final);
    // Deviatoric (shear) parity -- unchanged from the shear-dominated validation.
    check(close(q_final, qf, 5e-3), "q reaches the qf plateau (perfect plasticity at failure)");
    check(q50_err >= 0 && q50_err < 0.05, "secant at qf/2 = E50 (HS hyperbola, PLAXIS parity)");
    // Volumetric (Fig 15.4): initial contraction then net dilation.
    check(epsv_max_contr > 2e-4 && epsv_max_contr < 3e-3,
          "initial volumetric contraction (cap), small peak ~0.1% (Fig 15.4)");
    check(epsv_final < -2.5e-3,
          "net dilation at 5% axial strain (psi=6 deg), eps_v < -0.0025 (Fig 15.4 ~ -0.006)");
    check(epsv_final > -8e-3, "dilation magnitude physical (not runaway)");
}

} // namespace

int main() {
    test_berlin_triaxial();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil validated against Berlin Sand III (PLAXIS Fig 15.4: "
                    "qf, E50, AND volumetric contraction->dilation)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
