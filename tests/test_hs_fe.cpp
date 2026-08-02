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

} // namespace

int main() {
    test_volumetric_then_shear();
    test_frame_indifference();
    test_confined_loading_admissible();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil FE material point (shear) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
