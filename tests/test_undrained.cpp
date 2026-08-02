// Undrained (A) constitutive law verified against classical soil mechanics
// (Skempton's B parameter) and the undrained bulk modulus identity. Effective
// stiffness (E', nu') plus the pore-fluid bulk stiffness Kw/n = (3(nu_u-nu')/
// ((1-2nu_u)(1+nu'))) K' must give Ku = K' + Kw/n = E'(1+nu_u)/(3(1+nu')(1-2nu_u))
// and Skempton B = (Kw/n)/Ku. The undrained stiffness D_u = D' + (Kw/n) m m^T adds
// only volumetric (bulk) stiffness, leaving the shear modulus G' unchanged.
// (See docs/references/effective-stress-formulation.md.)
#include <katai/materials/material_model.hpp>

#include <cmath>
#include <cstdio>

using katai::core::MaterialModel;
using katai::core::MaterialType;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

void test_undrained_constitutive() {
    constexpr double E = 1.0e4, nu = 0.3, nu_u = 0.495;
    const MaterialModel m{MaterialType::LinearElastic, E, nu, 0, 0, 0};

    const double K_eff = E / (3.0 * (1.0 - 2.0 * nu));       // K'
    const double G = E / (2.0 * (1.0 + nu));                 // G' (unchanged)
    const double kwn = m.kw_over_n(nu_u);

    // (1) Kw/n must give the correct undrained bulk modulus.
    const double Ku = K_eff + kwn;
    const double Ku_exact = E * (1.0 + nu_u) / (3.0 * (1.0 + nu) * (1.0 - 2.0 * nu_u));
    check(close(Ku, Ku_exact), "Kw/n gives the correct undrained bulk modulus Ku");

    // (2) Skempton's B = (Kw/n)/Ku; for nu_u=0.495, nu'=0.3 -> Kw/n = 45 K', B = 45/46.
    const double B = kwn / Ku;
    check(close(kwn, 45.0 * K_eff, 1e-6), "Kw/n = 45 K' for nu_u=0.495, nu'=0.3");
    check(close(B, 45.0 / 46.0, 1e-6), "Skempton B = 45/46 ~ 0.978 (near-saturated)");
    std::printf("  undrained: Kw/n=%.1f*K'  Ku/K'=%.2f  Skempton B=%.4f\n",
                kwn / K_eff, Ku / K_eff, B);

    // (3) D_u = D' + (Kw/n) m m^T: symmetric; shear unchanged (G'); the volumetric
    //     stiffness rises to Ku. Pure shear strain produces NO mean-stress change.
    const Eigen::Matrix3d Du = m.undrained_plane_strain(nu_u);
    const Eigen::Matrix3d Dp = m.elastic_plane_strain();
    check((Du - Du.transpose()).cwiseAbs().maxCoeff() < 1e-9 * Du.cwiseAbs().maxCoeff(),
          "undrained D_u is symmetric");
    check(close(Du(2, 2), G) && close(Du(2, 2), Dp(2, 2)),
          "undrained: shear modulus G' unchanged");
    // Mean stress under unit volumetric strain [1,1,0]: (D_u row0 + row1)/2 dotted
    // with the bulk should reflect Ku (plane strain: sigma_xx = (K + 4G/3) e_x + ...).
    const Eigen::Vector3d ev(1.0, 1.0, 0.0);  // isotropic in-plane strain
    const Eigen::Vector3d s_u = Du * ev, s_p = Dp * ev;
    const double mean_u = 0.5 * (s_u(0) + s_u(1)), mean_p = 0.5 * (s_p(0) + s_p(1));
    check(close(mean_u - mean_p, 2.0 * kwn),
          "undrained: volumetric response increased by 2*(Kw/n) (water)");
}

} // namespace

int main() {
    test_undrained_constitutive();
    if (g_failures == 0) {
        std::printf("OK: undrained (A) constitutive law (Skempton B) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
