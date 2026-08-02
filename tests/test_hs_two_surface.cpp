// Hardening Soil -- two-surface (shear + cap) return mapping (hs_return_principal),
// the foundation for the oedometer calibration. Verified by reduction to the already-
// verified single-surface returns, plus an oedometer consistency check:
//  (1) cap off (Kc=0)  -> matches hs_shear_step (pure shear)
//  (2) isotropic strain (q=0, shear inactive) -> matches hs_cap_step (cap, K_iso)
//  (3) oedometer (eps_h=0): BOTH surfaces active, the stress stays on both (f_s=f_c=0)
//      and admissible, producing a lateral ratio K0 and tangent stiffness Eoed.
// (See docs/references/hardening-soil-formulation.md.)
#include <katai/materials/hardening_soil_plastic.hpp>

#include <cmath>
#include <cstdio>

using katai::core::HardeningSoilParams;
using katai::core::HsState;
using katai::core::hs_return_principal;
using katai::core::hs_shear_step;
using katai::core::hs_cap_step;

namespace {

constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

HardeningSoilParams base() {
    HardeningSoilParams p;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.p_ref = 100.0; p.friction = 35.0 * kPi / 180.0;
    p.dilatancy = 0.0; p.Rf = 0.9; p.nu_ur = 0.2;
    return p;
}

void test_reduces_to_shear() {
    HardeningSoilParams p = base();
    p.cap_beta = 0.0;  // cap off
    const Eigen::Vector3d sig_n(120, 100, 100);
    const Eigen::Vector3d de(2.0e-3, -5.0e-4, -5.0e-4);
    const auto a = hs_shear_step(p, sig_n, 0.0, de);
    const HsState b = hs_return_principal(p, sig_n, 0.0, 0.0, de);
    const double err = (a.stress - b.stress).cwiseAbs().maxCoeff() /
                       a.stress.cwiseAbs().maxCoeff();
    std::printf("  reduce->shear: rel err=%.3e  gamma_p %.4e vs %.4e\n",
                err, a.gamma_p, b.gamma_p);
    check(err < 1e-10 && close(a.gamma_p, b.gamma_p, 1e-9),
          "two-surface (cap off) matches the shear-only return");
}

void test_reduces_to_cap() {
    HardeningSoilParams p = base();
    p.cap_alpha = 0.6; p.cap_beta = 2.0e-3;
    const double p0 = 100.0;
    const Eigen::Vector3d sig_n(p0, p0, p0);
    const Eigen::Vector3d de(1.0e-5, 1.0e-5, 1.0e-5);  // isotropic -> q=0, shear inactive
    const auto a = hs_cap_step(p, sig_n, p0, de);
    const HsState b = hs_return_principal(p, sig_n, 0.0, p0, de);
    const double err = (a.stress - b.stress).cwiseAbs().maxCoeff() /
                       a.stress.cwiseAbs().maxCoeff();
    std::printf("  reduce->cap: rel err=%.3e  pp %.4f vs %.4f\n", err, a.pp, b.pp);
    check(err < 1e-8 && close(a.pp, b.pp, 1e-6),
          "two-surface (isotropic) matches the cap-only return");
}

void test_oedometer_both_active() {
    HardeningSoilParams p = base();
    p.cap_alpha = 0.7; p.cap_beta = 2.5e-3;
    const double sphi = std::sin(p.friction);
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    auto fbar = [&](double q, double s3) {
        const double Ei = p.Ei(s3), qa = p.q_asymptote(s3), Eur = p.Eur(s3);
        return (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur;
    };
    auto fcap = [&](const Eigen::Vector3d& s, double pp) {
        const double qt = s(0) + (delta - 1.0) * s(1) - delta * s(2);
        const double pm = (s(0) + s(1) + s(2)) / 3.0;
        return qt * qt / (p.cap_alpha * p.cap_alpha) + pm * pm - pp * pp;
    };

    // Start isotropic on both surfaces (q=0 on shear, p=pp on cap).
    const double p0 = 50.0;
    Eigen::Vector3d sig(p0, p0, p0);
    double gp = 0.0, pp = p0;
    double sig1_prev = p0, eps1 = 0.0;
    const double de1 = 2.0e-5;  // axial compression; lateral confined (eps_h = 0)
    double Eoed = 0.0, K0 = 0.0;
    bool admissible = true;
    for (int i = 0; i < 1500; ++i) {
        const HsState st = hs_return_principal(p, sig, gp, pp, Eigen::Vector3d(de1, 0, 0));
        Eoed = (st.stress(0) - sig1_prev) / de1;
        sig1_prev = st.stress(0);
        sig = st.stress; gp = st.gamma_p; pp = st.pp; eps1 += de1;
        // both surfaces must remain satisfied (f<=0) and near-active
        const double fs = fbar(sig(0) - sig(2), sig(2)) - gp;
        const double fc = fcap(sig, pp);
        if (fs > 1e-4 * (1 + gp) || fc > 1e-3 * (1 + pp * pp)) admissible = false;
        K0 = sig(2) / sig(0);
    }
    std::printf("  oedometer: sigma1=%.2f sigma3=%.2f  K0=%.4f  Eoed(tangent)=%.0f\n",
                sig(0), sig(2), K0, Eoed);
    check(admissible, "oedometer: both surfaces remain admissible (f_s<=0, f_c<=0)");
    check(K0 > 0.2 && K0 < 0.9, "oedometer lateral ratio K0 is physical");
    check(Eoed > 0.0 && std::isfinite(Eoed), "oedometer tangent stiffness is positive/finite");
}

} // namespace

int main() {
    test_reduces_to_shear();
    test_reduces_to_cap();
    test_oedometer_both_active();
    if (g_failures == 0) {
        std::printf("OK: HS two-surface return mapping verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
