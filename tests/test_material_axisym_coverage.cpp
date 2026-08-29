// EVERY registered constitutive model must actually be integrated, in BOTH modes.
//
// This gate exists because one was not. `integrate_point_axisym` dispatches on a closed enum
// with no `default:` branch, and when the Hoek-Brown model was added its axisymmetric case was
// simply not written. The switch then fell straight through: `trial` and `tangent` were left
// EXACTLY as the caller passed them in -- the previous iterate, handed back as if it had been
// integrated. Nothing said so. The equilibrium iteration merely failed to converge, and the run
// reported a stalled Newton search, which is what a difficult problem also looks like. The model
// had been verified twice at the material point and once through the whole FE path, and all
// three of those tests are plane strain.
//
// So the check is not "does the answer look right" -- other tests do that, model by model. It is
// the weaker and more general question a coverage gate should ask: was the output WRITTEN AT ALL.
// The outputs are poisoned with NaN first, so a branch that does not exist cannot leave behind
// something that merely looks plausible, and the loop walks the registry rather than a list
// maintained here, so a model registered tomorrow is covered without anyone remembering to.
//
// verify: none -- coverage mechanics (every registered model reaches an integrator in both
// modes and writes finite outputs). Correctness of each model's answer is its own test's job;
// this one cannot be satisfied by a plausible-looking wrong number.
#include <katai/materials/registry.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace kc = katai::core;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

// One parameter block that every model in the catalogue can be built from: each build function
// reads the fields it knows and ignores the rest, so filling all of them exercises the whole
// table without a per-model table here that could drift out of step with it.
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
    return p;
}

bool finite3(const Eigen::Vector3d& v) {
    return std::isfinite(v(0)) && std::isfinite(v(1)) && std::isfinite(v(2));
}
template <int N>
bool finite_mat(const Eigen::Matrix<double, N, N>& m) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (!std::isfinite(m(i, j))) return false;
    return true;
}

const double kNaN = std::nan("");

}  // namespace

int main() {
    std::printf("Every registered model is integrated in both modes\n\n");
    const std::vector<std::string_view> names = kc::model_names();
    check(!names.empty(), "the registry has models to walk");

    // A committed state under a modest all-round compression, and a strain increment with a
    // component in every direction -- including the HOOP, which is the one plane strain does
    // not have and the one an axisymmetric branch exists to carry.
    kc::GaussState committed;
    committed.stress = Eigen::Vector3d(-100.0, -120.0, -5.0);
    committed.stress_zz = -110.0;
    committed.pp = 200.0;                 // the soft-soil family needs a preconsolidation
    const Eigen::Vector3d de3(-2.0e-4, -3.0e-4, 1.0e-4);
    Eigen::Vector4d de4;
    de4 << -2.0e-4, -3.0e-4, 1.0e-4, -1.5e-4;

    for (std::string_view name : names) {
        const kc::ModelEntry* entry = kc::find_model(name);
        if (entry == nullptr) {
            check(false, std::string(name) + ": registered but not resolvable");
            continue;   // silent-drop-ok: the failed check above IS the report
        }
        const std::string why = entry->validate(every_field());
        if (!why.empty()) {
            check(false, std::string(name) + ": refuses the universal parameter block -- " + why);
            continue;   // silent-drop-ok: reported by the failed check above
        }
        const kc::MaterialModel m = entry->build(every_field());

        // POISON FIRST. A branch that never runs cannot then be mistaken for one that ran and
        // happened to return the previous stress -- which is exactly how the missing
        // Hoek-Brown case hid.
        kc::GaussState trial;
        trial.stress = Eigen::Vector3d(kNaN, kNaN, kNaN);
        trial.stress_zz = kNaN;
        Eigen::Matrix3d tangent3 = Eigen::Matrix3d::Constant(kNaN);
        kc::integrate_point(m, committed, de3, trial, tangent3);
        check(finite3(trial.stress) && std::isfinite(trial.stress_zz) &&
                  finite_mat<3>(tangent3),
              std::string(name) + ": plane strain writes a finite stress and tangent");

        kc::GaussState trial4;
        trial4.stress = Eigen::Vector3d(kNaN, kNaN, kNaN);
        trial4.stress_zz = kNaN;
        Eigen::Matrix4d tangent4 = Eigen::Matrix4d::Constant(kNaN);
        kc::integrate_point_axisym(m, committed, de4, trial4, tangent4);
        check(finite3(trial4.stress) && std::isfinite(trial4.stress_zz) &&
                  finite_mat<4>(tangent4),
              std::string(name) + ": axisymmetric writes a finite stress and tangent");

        // And the hoop must be CARRIED, not copied: a branch that quietly forwarded the
        // committed hoop stress would pass the finiteness check above. The increment has a
        // real hoop component, so sigma_theta cannot come back unchanged.
        check(std::fabs(trial4.stress_zz - committed.stress_zz) > 1e-9,
              std::string(name) + ": axisymmetric moves the HOOP stress the increment asked for");
    }

    if (g_failures == 0) {
        std::printf("\nOK: %zu model(s), both modes, every output written\n", names.size());
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
