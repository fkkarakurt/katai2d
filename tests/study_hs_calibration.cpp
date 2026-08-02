// HS cap calibration study (NOT a CI test; EXCLUDE_FROM_ALL). Diagnoses the K0^NC
// reproduction against the AUTHORITATIVE PLAXIS Material Models Manual closed forms:
//   - cap surface f_c = qtilde^2/M^2 + p^2 - pp^2          (MMM Eq 6-26)
//   - M <-> K0^NC via Brinkgreve (1994)                    (MMM Eq 10-13, Eq 6-29)
//       M = 3*sqrt[ (1-K0)^2/(1+2K0)^2
//                   + (1-K0)(1-2nu)(L-1) / ((1+2K0)(1-2nu)L - (1-K0)(1+nu)) ],  L = lambda*/kappa*
//     approx: M ~= 3.0 - 2.8 K0^NC                          (MMM Eq 10-14)
//   - for HS, lambda*/kappa* = Ks/Kc                        (Soft Soil = compression/swelling slope)
//       Ks/Kc ~= (Eur_ref/Eoed_ref) * K0 / ((1+2K0)(1-2nu)) (MMM Eq 6-30)
// Goal: show that the von Mises cap with alpha = M(Brinkgreve) reproduces K0^NC in a
// CAP-ONLY oedometer (proving von Mises == PLAXIS qtilde cap on the symmetric oedometer
// path), and quantify how much the shear mechanism (also active in HS) deviates K0 --
// the documented PLAXIS HS K0^NC limitation (MMM sec 6.4.3: K0^NC has a valid range; out
// of range values are rejected and replaced by the nearest possible value).
//   cmake --build ... --target study_hs_calibration
#include <katai/materials/hardening_soil_plastic.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::hs_integrate;
using katai::core::HsIntegrated;

namespace {
constexpr double kPi = 3.14159265358979323846;

// MMM Eq 6-30: Ks/Kc for HS (= lambda*/kappa*).
double ks_over_kc(const HardeningSoilParams& p, double K0) {
    return (p.Eur_ref / p.Eoed_ref) * K0 / ((1.0 + 2.0 * K0) * (1.0 - 2.0 * p.nu_ur));
}

// MMM Eq 10-13 (Brinkgreve 1994): exact M <-> K0^NC, with L = lambda*/kappa*.
double brinkgreve_M(double K0, double nu, double L) {
    const double a = (1.0 - K0) * (1.0 - K0) / ((1.0 + 2.0 * K0) * (1.0 + 2.0 * K0));
    const double num = (1.0 - K0) * (1.0 - 2.0 * nu) * (L - 1.0);
    const double den = (1.0 + 2.0 * K0) * (1.0 - 2.0 * nu) * L - (1.0 - K0) * (1.0 + nu);
    return 3.0 * std::sqrt(a + num / den);
}

// FULL-model confined (oedometer) probe (shear+cap, robust hs_integrate). Strain-driven
// de1, de2=de3=0. Reports the lateral stresses sig2,sig3 separately (the qtilde cap can
// split them) and Eoed at sigma1 ~ p_ref.
void full_oedometer(const HardeningSoilParams& p, double& s1, double& s2, double& s3,
                    double& Eoed) {
    const double pr = p.p_ref, p0 = 0.02 * pr;
    Eigen::Vector3d sig(p0, p0, p0);
    double gp = 0.0, pp = p0, s1b = p0;
    const double de1 = 1.0e-4;
    s1 = s2 = s3 = Eoed = 0.0;
    for (int i = 0; i < 6000; ++i) {
        const HsIntegrated r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, 0, 0));
        if (s1b < pr && r.stress(0) >= pr) {
            Eoed = (r.stress(0) - s1b) / de1;
            s1 = r.stress(0); s2 = r.stress(1); s3 = r.stress(2);
        }
        s1b = r.stress(0); sig = r.stress; gp = r.gamma_p; pp = r.pp;
        if (sig(0) > 1.25 * pr) break;
    }
}

void diagnose(const char* name, HardeningSoilParams p, double K0nc) {
    const double L = ks_over_kc(p, K0nc);
    const double M_exact = brinkgreve_M(K0nc, p.nu_ur, L);
    const double M_approx = 3.0 - 2.8 * K0nc;
    std::printf("\n=== %s (target K0nc=%.4f) ===\n", name, K0nc);
    std::printf("  Ks/Kc (Eq 6-30) = %.4f ; M_exact (Eq 10-13) = %.4f ; M_approx (Eq 10-14) = %.4f\n",
                L, M_exact, M_approx);

    const double Ks_ref = p.Eur_ref / (3.0 * (1.0 - 2.0 * p.nu_ur));
    const double beta_plaxis = p.p_ref * (L - 1.0) / Ks_ref;
    const double Kp = katai::core::hs_cap_Kp(p, K0nc);
    std::printf("  Ks_ref=%.0f beta_plaxis=%.3e  beta_Kp(pref/Kp)=%.3e\n",
                Ks_ref, beta_plaxis, p.p_ref / Kp);
    // von Mises cap (= axisymmetric qtilde), M from Brinkgreve, over a beta sweep.
    std::printf("  --- von Mises cap ; M=M_approx=%.4f --- (Eoed=Eoed_ref <-> K0 floor)\n", M_approx);
    for (double b : {beta_plaxis, p.p_ref / Kp, 2.0 * beta_plaxis, 5.0 * beta_plaxis}) {
        HardeningSoilParams q = p; q.cap_alpha = M_approx; q.cap_beta = b;
        double s1, s2, s3, eo; full_oedometer(q, s1, s2, s3, eo);
        std::printf("     beta=%.3e -> K0=%.4f  Eoed=%.0f (t %.0f)\n",
                    b, 0.5 * (s2 + s3) / s1, eo, p.Eoed_ref);
    }
}
} // namespace

// Berlin Sand III drained triaxial at sigma3=200 WITH cap on -- prints q, axial stress and
// VOLUMETRIC strain vs axial strain to compare against PLAXIS Fig 15.4 (contraction then
// dilation, psi=6deg). Cap shape M from Brinkgreve (PLAXIS's own cap), pp initialised to the
// isotropic consolidation stress sigma3.
void berlin_triaxial_capon(double M, double beta) {
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 105e3; p.Eur_ref = 315e3; p.Eoed_ref = 105e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38 * kPi / 180;
    p.dilatancy = 6 * kPi / 180; p.cohesion = 1.0; p.Rf = 0.9;
    p.cap_alpha = M; p.cap_beta = beta;
    const double sigma3 = 200.0, qf = p.q_failure(sigma3);
    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0, pp = sigma3, eps1 = 0, epsv = 0;
    std::printf("  Berlin triaxial cap-on (M=%.3f beta=%.3e) qf=%.1f:\n", M, beta, qf);
    std::printf("     eps1      q        axial_sig   eps_v\n");
    double next = 0.0;
    for (int s = 0; s < 12000; ++s) {
        double dlat = 0.0; HsIntegrated r;
        for (int it = 0; it < 30; ++it) {
            r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(5e-6, dlat, dlat));
            const double g = r.stress(2) - sigma3;
            if (std::fabs(g) <= 1e-8 * (1 + sigma3)) break;
            dlat -= g / (r.tangent(2, 1) + r.tangent(2, 2));
        }
        epsv += 5e-6 + 2.0 * dlat; sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += 5e-6;
        if (eps1 >= next - 1e-12) {
            std::printf("     %.4f   %7.1f   %7.1f   % .5f\n",
                        eps1, sig(0) - sig(2), sig(0), epsv);
            next += 0.005;
        }
        if (eps1 > 0.05) break;
    }
}

int main() {
    HardeningSoilParams sand;
    sand.p_ref = 100; sand.E50_ref = 105e3; sand.Eur_ref = 315e3; sand.Eoed_ref = 105e3;
    sand.m = 0.55; sand.nu_ur = 0.2; sand.friction = 38 * kPi / 180;
    sand.dilatancy = 6 * kPi / 180; sand.cohesion = 1.0; sand.Rf = 0.9;
    diagnose("Berlin Sand III", sand, 0.38);
    // Berlin triaxial volumetric vs Fig 15.4: sweep cap aspect M (with calibrated beta from
    // hs_calibrate_cap so Eoed matches), find which gives net dilation (psi=6deg).
    {
        HardeningSoilParams cal = sand;
        katai::core::hs_calibrate_cap(cal, 0.38);
        std::printf("\n  [calibrated cap_alpha=%.3f cap_beta=%.3e]\n", cal.cap_alpha, cal.cap_beta);
        for (double M : {cal.cap_alpha, 5.0, 20.0}) {
            berlin_triaxial_capon(M, cal.cap_beta);
        }
    }

    HardeningSoilParams clay;
    clay.p_ref = 100; clay.E50_ref = 1.2e4; clay.Eur_ref = 6e4; clay.Eoed_ref = 8e3;
    clay.m = 1.0; clay.nu_ur = 0.2; clay.friction = 24 * kPi / 180;
    clay.dilatancy = 0; clay.cohesion = 5; clay.Rf = 0.9;
    diagnose("Clay (m=1)", clay, 1.0 - std::sin(24 * kPi / 180));
    return 0;
}
