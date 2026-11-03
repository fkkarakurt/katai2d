// THE STRESS A MODEL RETURNS DOES NOT DEPEND ON WHETHER ITS TANGENT WAS ASKED FOR.
//
// Two of the solvers lean on this. The equilibrium iteration of a Plastic phase tries the full
// Newton step with the tangent assembled when the iteration before took its full step, and
// without it otherwise; a consolidation step tries the lengths of its line search without the
// element matrices and evaluates the length it accepts again, with them. Both are only the same
// iteration as before if the internal force and the state a point is left in are the same, to the
// bit, whether the point was integrated with TangentMode::kNone and no report, or with the
// continuum or the consistent tangent and a report. A model whose stress came to depend on the
// mode -- a cheaper walk when the tangent is not wanted, a different subdivision -- would make the
// iterates differ without anything failing, and the answer would move at the tolerance.
//
// NGI-ADP is the model that does not hold it: asked for its tangent, it subdivides its return one
// level less (adp::integrate_tangent calls adp::integrate at depth 1), and its too-small-to-return
// shortcut then applies at another level. Both solvers keep evaluating every trial with the tangent
// when a material of that model is present, so it is checked here and reported, not required.
//
// The registry is walked, so a model registered tomorrow is covered without anyone remembering
// to add it. Each model is integrated from a compressed state through three increments -- small,
// medium and large, the large one well into plastic flow for the soils -- in plane strain and in
// axisymmetry, and every field of the state it returns is compared byte for byte.
//
// verify: none -- an invariant of the integrators that the solvers' trial evaluations rely on;
// the correctness of each model's answer is its own test's job.
#include <katai/materials/registry.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace kc = katai::core;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

// The parameter block of test_material_axisym_coverage: every build function reads the fields it
// knows and ignores the rest.
kc::MaterialParams every_field() {
    kc::MaterialParams p;
    p.E = 3.0e4; p.nu = 0.25;
    p.c = 10.0; p.phi_rad = 0.5236; p.psi_rad = 0.0;    // 30 deg
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.p_ref = 100.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.G0_ref = 1.2e5; p.gamma07 = 1.5e-4;
    p.lam_star = 0.1; p.kap_star = 0.02; p.mu_star = 0.005;
    p.sig_ci = 5.0e4; p.mi = 10.0; p.gsi = 50.0; p.hb_D = 0.0; p.sig_psi = 0.0;
    p.k0nc_auto = true;
    p.adp_suP = 0.4; p.adp_suDSS = 0.7; p.adp_gfC = 0.02; p.adp_gfE = 0.06; p.adp_gfDSS = 0.04;
    return p;
}

kc::MaterialParams every_field_for(std::string_view name) {
    kc::MaterialParams p = every_field();
    if (name == "NgiAdp") { p.drainage = kc::DrainageClass::UndrainedC; p.nu = 0.495; }
    if (name == "JointedRock") {
        p.drainage = kc::DrainageClass::Drained;
        p.jr_aniso = true; p.jr_E2 = 0.5 * p.E; p.jr_nu2 = 0.2; p.jr_G2 = 0.3 * p.E;
        p.jr_sets = 2; p.jr_alpha[0] = 0.4; p.jr_alpha[1] = -0.9;
        p.jr_c[1] = 5.0; p.jr_phi[1] = 0.5; p.jr_psi[1] = 0.1;
    }
    return p;
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
template <int N>
bool same_bits(const Eigen::Matrix<double, N, 1>& a, const Eigen::Matrix<double, N, 1>& b) {
    for (int i = 0; i < N; ++i)
        if (!same_bits(a(i), b(i))) return false;
    return true;
}

// Every field of the state, byte for byte; the name of the first that differs, or "".
std::string first_difference(const kc::GaussState& a, const kc::GaussState& b) {
    if (!same_bits<3>(a.stress, b.stress)) return "stress";
    if (!same_bits(a.stress_zz, b.stress_zz)) return "stress_zz";
    if (!same_bits(a.eps_vol, b.eps_vol)) return "eps_vol";
    if (!same_bits(a.gamma_p, b.gamma_p)) return "gamma_p";
    if (!same_bits(a.pp, b.pp)) return "pp";
    if (!same_bits(a.gamma_hist, b.gamma_hist)) return "gamma_hist";
    if (!same_bits<4>(a.hs_hist, b.hs_hist)) return "hs_hist";
    if (!same_bits(a.hs_g_mono, b.hs_g_mono)) return "hs_g_mono";
    if (!same_bits(a.eps_vol_und, b.eps_vol_und)) return "eps_vol_und";
    if (!same_bits(a.pw_carried, b.pw_carried)) return "pw_carried";
    if (a.plasticity_excluded != b.plasticity_excluded) return "plasticity_excluded";
    if (!same_bits<4>(a.adp_s0, b.adp_s0)) return "adp_s0";
    if (a.adp_set != b.adp_set) return "adp_set";
    return "";
}

constexpr kc::TangentMode kModes[] = {kc::TangentMode::kContinuum, kc::TangentMode::kConsistent};
const char* mode_name(kc::TangentMode t) {
    return t == kc::TangentMode::kContinuum ? "continuum" : "consistent";
}

}  // namespace

int main() {
    std::printf("The stress a model returns does not depend on the tangent asked for\n\n");
    const std::vector<std::string_view> names = kc::model_names();
    check(!names.empty(), "the registry has models to walk");

    // A compressed state with a preconsolidation the soft-soil family can be loaded past, and
    // increments that keep the same direction -- compression with shear -- at three sizes.
    kc::GaussState committed;
    committed.stress = Eigen::Vector3d(-100.0, -120.0, -5.0);
    committed.stress_zz = -110.0;
    committed.pp = 200.0;
    const double scale[] = {1.0e-5, 1.0e-3, 2.0e-2};
    const char* size_name[] = {"small", "medium", "large"};
    const double dt_day = 1.0;   // the creep model reads it; the others do not

    int held_out = 0;
    for (std::string_view name : names) {
        const kc::ModelEntry* entry = kc::find_model(name);
        if (entry == nullptr) {
            check(false, std::string(name) + ": registered but not resolvable");
            continue;   // silent-drop-ok: the failed check above IS the report
        }
        const std::string why = entry->validate(every_field_for(name));
        if (!why.empty()) {
            check(false, std::string(name) + ": refuses the universal parameter block -- " + why);
            continue;   // silent-drop-ok: reported by the failed check above
        }
        const kc::MaterialModel m = entry->build(every_field_for(name));
        const bool required = name != "NgiAdp";
        std::string worst;   // the first difference met, for the report
        for (int s = 0; s < 3; ++s) {
            const Eigen::Vector3d de3 = scale[s] * Eigen::Vector3d(-0.4, -1.0, 0.3);
            Eigen::Vector4d de4;
            de4 << de3(0), de3(1), de3(2), -0.3 * scale[s];
            // Without the tangent and without a report: what a trial evaluation asks for.
            kc::GaussState lean3, lean4;
            Eigen::Matrix3d t3;
            Eigen::Matrix4d t4;
            kc::integrate_point(m, committed, de3, lean3, t3, kc::TangentMode::kNone, dt_day);
            kc::integrate_point_axisym(m, committed, de4, lean4, t4, kc::TangentMode::kNone, dt_day);
            for (const kc::TangentMode mode : kModes) {
                kc::GaussState full3, full4;
                kc::PointReport rep3;
                kc::PointReport4 rep4;
                kc::integrate_point(m, committed, de3, full3, t3, mode, dt_day, &rep3);
                kc::integrate_point_axisym(m, committed, de4, full4, t4, mode, dt_day, &rep4);
                const std::string d3 = first_difference(lean3, full3);
                const std::string d4 = first_difference(lean4, full4);
                const std::string where = std::string(name) + ", " + size_name[s] + " increment, " +
                                          mode_name(mode) + " tangent";
                if (required) {
                    check(d3.empty(), where + ", plane strain: the same state" +
                                          (d3.empty() ? "" : " (differs in " + d3 + ")"));
                    check(d4.empty(), where + ", axisymmetric: the same state" +
                                          (d4.empty() ? "" : " (differs in " + d4 + ")"));
                } else if (worst.empty() && (!d3.empty() || !d4.empty())) {
                    worst = where + (d3.empty() ? ", axisymmetric, " + d4 : ", plane strain, " + d3);
                }
            }
        }
        if (!required) {
            ++held_out;
            std::printf("note: %s is held out of the leaner trials by both solvers; here %s\n",
                        std::string(name).c_str(),
                        worst.empty() ? "it returned the same state on these increments"
                                      : ("it differed: " + worst).c_str());
        }
    }

    if (g_failures == 0) {
        std::printf("\nOK: %zu model(s), %d held out, three increments, both kinematics\n",
                    names.size(), held_out);
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
