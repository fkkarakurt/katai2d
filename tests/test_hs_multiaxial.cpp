// Hardening Soil away from the triaxial compression corner.
//
// Every earlier HS test walks the triaxial COMPRESSION corner (sigma2 = sigma3) or the oedometer,
// and the model had been written for exactly that corner: one shear flow direction (1, R, R) and
// a von Mises cap. On that corner both are right. Off it they are not, and nothing measured it:
// an undrained triaxial EXTENSION of a normally consolidated clay came out 2 to 2.6 times too
// strong, with the mean effective stress DOUBLING on a psi = 0 path where it cannot rise at all.
// The stress path under an excavation floor is that path.
//
// The model this checks is the one the literature defines (Schanz 1998, Ch. 4; Schanz, Vermeer &
// Bonnier 1999; Benz 2007, Eqns 7.15-7.23):
//   * shear yield functions f_ij = fbar(sigma_i - sigma_j) - gamma_p on the Mohr-Coulomb planes,
//     with the corners of the hexagon treated as the meeting of two planes;
//   * plastic potentials g_ij = (sigma_i - sigma_j)/2 - (sigma_i + sigma_j)/2 sin(psi_m), so that
//     d eps_v^p = sin(psi_m) d gamma_p  (Schanz 1998, Eqn 4.40) on every plane and at every corner;
//   * the cap measured with q~ = sigma1 + (delta - 1) sigma2 - delta sigma3,
//     delta = (3 + sin phi)/(3 - sin phi)  (Schanz 1998, Eqns 4.71-4.74), which reduces to q on the compression
//     corner -- so the oedometer and every compression result is untouched by construction.
//
// Compression positive inside this file, as in the integrator.
//
// verify: KV-CST-018
//   oracle:   closed_form
//   source:   Schanz, T. (1998). Zur Modellierung des mechanischen Verhaltens von Reibungsmaterialien. Mitteilung 45, Institut fuer Geotechnik, Universitaet Stuttgart, Ch. 4; Benz, T. (2007). Small-strain stiffness of soils and its numerical consequences. Mitteilung 55, Institut fuer Geotechnik, Universitaet Stuttgart, Sec. 7.2
//   locator:  Schanz (1998) Eqns 4.34-4.40 / Benz (2007) Eqns 7.15-7.16 and 7.22: shear yield f_ij = fbar(sigma_i - sigma_j) - gamma_p, plastic potential g_ij = (sigma_i - sigma_j)/2 - (sigma_i + sigma_j)/2 sin psi_m, d gamma_p = d lambda (dg/ds1 - dg/ds2 - dg/ds3) = d lambda per active plane, hence d eps_v^p = sin psi_m d gamma_p on a face and at either corner, and no plastic strain along the intermediate stress on a face
//   quantity: the plastic strain increment d eps^p = d eps - C_e d sigma of one drained increment from a state on the shear surface (cap far away), at the compression corner, the extension corner and a general state [-]
//   expected: d eps_v^p / d gamma_p = -sin psi_m (Rowe, phi = 35, psi = 10 at 80% of q_f), d eps_2^p = 0 on the face, sigma2 = sigma3 (resp. sigma1 = sigma2) preserved at the corner
//   band:     2e-3 on the ratio and 1e-3 of the largest component on the intermediate strain, as asserted below -- measured exact to five decimals (-0.09633 against -0.09633) in all three places, and 1e-20 on the intermediate strain
//
// verify: KV-CST-019
//   oracle:   published_benchmark
//   source:   Borgh, R. (2018). Sensitivity analysis of numerical models for diaphragm walls with cross-walls in soft soils. Master's thesis, Department of Architecture and Civil Engineering, Chalmers University of Technology, Gothenburg
//   locator:  Sec. 3.5.2 Figs 3.18-3.19 (undrained triaxial element tests of the HSsmall clay sets, compression and extension) with the parameter sets of App. D (Clays 4-7: E50, Eoed, Eur, m = 1, nu_ur = 0.15, p_ref = 100, c', phi', psi = 0, gamma_0.7, G0, K0nc, OCR / POP) and the in-situ K0 of Table 3.8; start states and end points read off Fig. 3.18
//   quantity: mean effective stress p' and deviator q at 3% axial strain after K0 consolidation, in undrained compression and in undrained extension [kPa]
//   expected: (p', q) compression / extension: Clay4 (80.4, 106.4) / (81.0, 75.0); Clay5 (117.7, 144.7) / (118.6, 104.7); Clay6 (154.7, 190.6) / (156.6, 138.7); Clay7 (252.4, 309.8) / (255.3, 223.0)
//   band:     10% on each, as asserted below -- measured q -4.9 to -7.7% in compression and -1.9 to -4.7% in extension, p' +0.9 to +4.2%; the compression deviator is systematically low and declared as open: the published paths reach the Mohr-Coulomb line by 3% while this model is at 0.89-0.92 of it, which points at the cap calibration (alpha, beta) and the initial gamma_p convention, not at the flow rule this test was written for. Before the multi-surface formulation the extension deviator was 2.0 to 2.6 times the reference with p' doubling
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>

#include <cmath>
#include <cstdio>

using namespace katai::core;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what);
    if (!ok) ++g_failures;
}

struct Clay {
    const char* name;
    double E50, Eoed, Eur, c, phi, g07, G0, k0nc, ocr, pop;
    double k0, sv;          // consolidation K0 and sigma_v' of the reference test
    double pc_ref, qc_ref;  // (p', q) at 3% axial strain in compression (reference curves)
    double pe_ref, qe_ref;  // and in extension
};

MaterialModel build(const Clay& k) {
    MaterialParams p;
    p.E50_ref = k.E50; p.Eoed_ref = k.Eoed; p.Eur_ref = k.Eur; p.m = 1.0; p.p_ref = 100;
    p.Rf = 0.9; p.nu_ur = 0.15; p.c = k.c; p.phi_rad = k.phi * kPi / 180; p.psi_rad = 0;
    p.G0_ref = k.G0; p.gamma07 = k.g07; p.k0nc_auto = false; p.k0nc = k.k0nc;
    p.tension_cutoff = true; p.tensile_strength = 0;
    p.E = k.Eur; p.nu = 0.15;
    return find_model(k.G0 > 0 ? "HSsmall" : "HardeningSoil")->build(p);
}

// K0-consolidated state (axisymmetric: [r, z, rz] + hoop, tension positive), seeded by the same
// function the initial-stress procedure seeds with.
GaussState seed(const MaterialModel& m, double sv, double sh, const Clay& k) {
    GaussState gs;
    gs.stress = Eigen::Vector3d(-sh, -sv, 0.0);
    gs.stress_zz = -sh;
    const auto pe = hs_small_strain_params(m.hs, 0.0);
    const HsSeed s = k.pop > 0 ? hs_seed_from_history(pe, sh, sv, 0.0, sh, 2, 1.0, k.pop)
                               : hs_seed_from_history(pe, sh, sv, 0.0, sh, 1, k.ocr, 0.0);
    gs.pp = s.pp;
    gs.gamma_p = s.gamma_p;
    return gs;
}

struct Path { double q3, p3, max_asym, max_p_rise, p0; };

// Undrained (constant volume) axisymmetric triaxial to 3% axial strain: eps_r = eps_theta =
// -eps_z/2. dir = +1 compression (axial shortening), -1 extension.
Path undrained(const MaterialModel& m, const Clay& k, int dir) {
    GaussState gs = seed(m, k.sv, k.k0 * k.sv, k);
    const int n = 3000;
    const double d = 0.03 / n;
    const double p0 = k.sv * (1.0 + 2.0 * k.k0) / 3.0;
    Path out{0, 0, 0, 0, p0};
    double p_prev = p0;
    for (int i = 1; i <= n; ++i) {
        GaussState tr = gs;
        Eigen::Matrix4d T;
        const Eigen::Vector4d de(dir * d / 2, -dir * d, 0.0, dir * d / 2);
        integrate_point_axisym(m, gs, de, tr, T, TangentMode::kContinuum);
        gs = tr;
        const double sr = -gs.stress(0), sz = -gs.stress(1), sth = -gs.stress_zz;
        const double p = (sr + sz + sth) / 3.0;
        out.max_asym = std::max(out.max_asym, std::fabs(sr - sth) / p);
        out.max_p_rise = std::max(out.max_p_rise, (p - p_prev) / p0);
        p_prev = p;
        out.q3 = std::fabs(sz - sr);
        out.p3 = p;
    }
    return out;
}

void test_reference_clays() {
    // Clays 4-7 of a published deep-excavation study in soft clay (HSsmall, m = 1,
    // nu_ur = 0.15), each K0-consolidated at its own in-situ K0 to the sigma_v' of its depth,
    // then sheared undrained in compression and in extension. The start states and the (p', q)
    // reached at 3% axial strain are read off the published element-test stress paths of the same
    // parameter sets (the start points reproduce the in-situ K0 to 0.01, which is how the
    // consolidation K0 was identified).
    const Clay clays[] = {
        {"Clay4", 15400, 8300, 30800, 1, 32.0, 4e-4, 39230, 0.53, 1.00, 32, 0.61, 126.5, 80.4, 106.4, 81.0, 75.0},
        {"Clay5", 13100, 6600, 26200, 0, 30.5, 6e-4, 40960, 0.53, 1.00, 32, 0.58, 186.3, 117.7, 144.7, 118.6, 104.7},
        {"Clay6", 13100, 6600, 26200, 0, 30.5, 8e-4, 36060, 0.53, 1.14, 0, 0.57, 248.8, 154.7, 190.6, 156.6, 138.7},
        {"Clay7", 10500, 5285, 21000, 0, 30.5, 9e-4, 37460, 0.53, 1.14, 0, 0.56, 408.0, 252.4, 309.8, 255.3, 223.0},
    };
    for (const Clay& k : clays) {
        const MaterialModel m = build(k);
        const Path c = undrained(m, k, +1), e = undrained(m, k, -1);
        std::printf("  %s  p0' %.1f | compression (p', q) = (%.1f, %.1f) ref (%.1f, %.1f) |"
                    " extension (%.1f, %.1f) ref (%.1f, %.1f) | asym %.1e, p' rise %.1e\n",
                    k.name, c.p0, c.p3, c.q3, k.pc_ref, k.qc_ref, e.p3, e.q3, k.pe_ref, k.qe_ref,
                    e.max_asym, e.max_p_rise);
        char msg[160];
        std::snprintf(msg, sizeof msg, "%s: undrained compression q within 10%% of the reference", k.name);
        check(std::fabs(c.q3 / k.qc_ref - 1) < 0.10, msg);
        std::snprintf(msg, sizeof msg, "%s: undrained compression p' within 10%% of the reference", k.name);
        check(std::fabs(c.p3 / k.pc_ref - 1) < 0.10, msg);
        std::snprintf(msg, sizeof msg, "%s: undrained extension q within 10%% of the reference", k.name);
        check(std::fabs(e.q3 / k.qe_ref - 1) < 0.10, msg);
        std::snprintf(msg, sizeof msg, "%s: undrained extension p' within 10%% of the reference", k.name);
        check(std::fabs(e.p3 / k.pe_ref - 1) < 0.10, msg);
        std::snprintf(msg, sizeof msg, "%s: extension keeps sigma_r = sigma_theta", k.name);
        check(e.max_asym < 1e-6, msg);
        std::snprintf(msg, sizeof msg, "%s: psi = 0 undrained extension never raises p'", k.name);
        check(e.max_p_rise < 1e-6, msg);
        // The effective stress path cannot cross the Mohr-Coulomb extension line:
        // q <= 6 sin(phi) (p' + c cot phi) / (3 + sin phi).
        const double sphi = std::sin(k.phi * kPi / 180), ccot = k.c * std::cos(k.phi * kPi / 180) / sphi;
        std::snprintf(msg, sizeof msg, "%s: extension stays inside the Mohr-Coulomb extension line", k.name);
        check(e.q3 <= 6 * sphi * (e.p3 + ccot) / (3 + sphi) * (1 + 1e-3), msg);
    }
}

// The flow rule, measured on the plastic strain itself: d eps^p = d eps - C_e d sigma, with the
// elastic compliance at the start-of-increment E_ur (sigma3 held constant, so E_ur is constant).
// Drained principal-space states with the cap far away (OCR high) so only shear yields.
void test_flow_rule() {
    MaterialParams p;
    p.E50_ref = 30000; p.Eoed_ref = 30000; p.Eur_ref = 90000; p.m = 0.5; p.p_ref = 100;
    p.Rf = 0.9; p.nu_ur = 0.2; p.c = 0; p.phi_rad = 35 * kPi / 180; p.psi_rad = 10 * kPi / 180;
    p.k0nc_auto = true;
    const MaterialModel m = find_model("HardeningSoil")->build(p);
    const HardeningSoilParams& h = m.hs;
    const double sphi = std::sin(p.phi_rad), spsi = std::sin(p.psi_rad);
    const double scv = (sphi - spsi) / (1 - sphi * spsi);
    auto rowe = [&](double s1, double s3) {
        const double sm = (s1 - s3) / (s1 + s3);
        const double r = (sm - scv) / (1 - sm * scv);
        return sm < 0.75 * sphi ? 0.0 : std::max(r, 0.0);
    };
    struct St { const char* name; Eigen::Vector3d s; Eigen::Vector3d de; };
    const double s3 = 100, qf = h.q_failure(s3), q = 0.8 * qf, d = 1e-7;
    const St states[] = {
        {"compression corner (s2 = s3)", Eigen::Vector3d(s3 + q, s3, s3), Eigen::Vector3d(d, -0.3 * d, -0.3 * d)},
        {"extension corner (s1 = s2)", Eigen::Vector3d(s3 + q, s3 + q, s3), Eigen::Vector3d(0.3 * d, 0.3 * d, -d)},
        {"general state (s1 > s2 > s3)", Eigen::Vector3d(s3 + q, s3 + 0.4 * q, s3), Eigen::Vector3d(d, 0.0, -0.3 * d)},
    };
    for (const St& st : states) {
        const double gp0 = [&] {  // put the state ON its shear surface
            const double Ei = h.Ei(s3), qa = h.q_asymptote(s3), Eur = h.Eur(s3);
            return (2 / Ei) * q / (1 - q / qa) - 2 * q / Eur;
        }();
        const double pp = 1e6;  // cap far away
        const HsIntegrated r = hs_integrate(h, st.s, gp0, pp, st.de);
        const Eigen::Vector3d ds = r.stress - st.s;
        const double Eur = h.Eur(s3), nu = h.nu_ur;
        Eigen::Matrix3d Ce;
        Ce << 1, -nu, -nu, -nu, 1, -nu, -nu, -nu, 1;
        Ce /= Eur;
        const Eigen::Vector3d dep = st.de - Ce * ds;
        const double dev = dep.sum();
        const double dgp = r.gamma_p - gp0;
        const double s_expect = rowe(st.s(0), st.s(2));
        std::printf("  %-30s  d eps^p = (%.3e, %.3e, %.3e)  d gamma_p %.3e  eps_v/gamma %.5f"
                    " (-sin psi_m = %.5f)\n", st.name, dep(0), dep(1), dep(2), dgp, dev / dgp, -s_expect);
        char msg[160];
        std::snprintf(msg, sizeof msg, "%s: d eps_v^p = sin(psi_m) d gamma_p", st.name);
        check(dgp > 0 && std::fabs(dev / dgp + s_expect) < 2e-3, msg);
        if (st.s(1) < st.s(0) - 1 && st.s(1) > st.s(2) + 1) {
            std::snprintf(msg, sizeof msg, "%s: no plastic strain along the intermediate stress", st.name);
            check(std::fabs(dep(1)) < 1e-3 * dep.cwiseAbs().maxCoeff(), msg);
        }
        if (std::fabs(st.s(1) - st.s(2)) < 1e-9) {
            std::snprintf(msg, sizeof msg, "%s: the corner stays a corner (s2 = s3)", st.name);
            check(std::fabs(r.stress(1) - r.stress(2)) < 1e-9 * st.s(0), msg);
        }
        if (std::fabs(st.s(0) - st.s(1)) < 1e-9) {
            std::snprintf(msg, sizeof msg, "%s: the corner stays a corner (s1 = s2)", st.name);
            check(std::fabs(r.stress(0) - r.stress(1)) < 1e-9 * st.s(0), msg);
        }
    }
}
}  // namespace

int main() {
    std::printf("HS flow rule off and on the corners\n");
    test_flow_rule();
    std::printf("HSsmall undrained triaxial, compression and extension (reference clays)\n");
    test_reference_clays();
    if (g_failures) std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
