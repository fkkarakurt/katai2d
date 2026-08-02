// Constitutive catalogue (katai/materials/registry.hpp) -- the construction seam.
//
// The registry replaced the analysis driver's inline SoilModel switch, and the
// contract of that replacement is bit-identity: build() must produce exactly
// the MaterialModel the driver used to assemble by hand. So the oracle here is
// that former inline construction, copied verbatim, and every comparison is
// == on doubles -- not a tolerance. The remaining cases pin the seam's other
// promises: honest drainage refusals with their exact user-facing messages,
// resolution by name with refusal of the unknown, and external registration
// that can never displace a shipped model (R2: composition must not be able
// to change an answer).
// verify: KV-NUM-001
//   oracle:   independent_path
//   source:   the analysis driver's former inline SoilModel construction (build_problem.hpp before commit b35bb9a), duplicated verbatim in this test as the oracle
//   locator:  oracle_common() and the per-model expected blocks in this file, field for field
//   quantity: every field of the constructed MaterialModel, per shipped model and drainage class
//   expected: bitwise equality with the former inline construction (== on doubles, no tolerance)
//   band:     exact -- construction is plain arithmetic; any deviation means the registry changed a validated number
#include <katai/materials/hardening_soil.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

using katai::core::DrainageClass;
using katai::core::MaterialModel;
using katai::core::MaterialParams;
using katai::core::MaterialType;
using katai::core::ModelEntry;
using katai::core::find_model;
using katai::core::model_names;
using katai::core::register_model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// A representative parameter set touching every field, angles pre-converted
// exactly as the driver converts them (acos(-1) * deg / 180).
MaterialParams sample_params() {
    const double kPi = std::acos(-1.0);
    MaterialParams p;
    p.E = 1.3e4; p.nu = 0.3; p.c = 5.0;
    p.phi_rad = 30.0 * kPi / 180.0;
    p.psi_rad = 5.0 * kPi / 180.0;
    p.tension_cutoff = true;
    p.tensile_strength = 2.5;
    p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.55; p.p_ref = 100.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.G0_ref = 1.2e5; p.gamma07 = 1.5e-4;
    p.lam_star = 0.10; p.kap_star = 0.02; p.mu_star = 0.005;
    p.k0nc_auto = true; p.k0nc = 0.5;
    p.drainage = DrainageClass::Drained;
    return p;
}

// The driver's former inline construction, kept as the oracle. Field for
// field, operation for operation -- see build_problem.hpp before the registry.
MaterialModel oracle_common(const MaterialParams& p) {
    MaterialModel mm;
    mm.youngs_modulus = p.E; mm.poisson_ratio = p.nu;
    mm.cohesion = p.c; mm.friction_angle = p.phi_rad; mm.dilatancy_angle = p.psi_rad;
    mm.tension_cutoff = p.tension_cutoff;
    mm.tensile_strength = std::max(0.0, p.tensile_strength);
    if (p.drainage == DrainageClass::UndrainedA || p.drainage == DrainageClass::UndrainedB) {
        mm.undrained = true; mm.undrained_poisson = 0.495;
    }
    return mm;
}

bool same_common(const MaterialModel& a, const MaterialModel& b) {
    return a.type == b.type && a.youngs_modulus == b.youngs_modulus &&
           a.poisson_ratio == b.poisson_ratio && a.cohesion == b.cohesion &&
           a.friction_angle == b.friction_angle && a.dilatancy_angle == b.dilatancy_angle &&
           a.tension_cutoff == b.tension_cutoff && a.tensile_strength == b.tensile_strength &&
           a.undrained == b.undrained && a.undrained_poisson == b.undrained_poisson;
}

void test_catalogue() {
    std::printf("catalogue: the six shipped models resolve by name, nothing else does\n");
    // profile_E / profile_c pin the former inline gating of build_profiles: E'_inc only
    // where E' itself is an input (LE, MC), c'_inc only where the MC strength is.
    const struct { const char* name; MaterialType type; bool nl, hard, soft, pe, pc; } expect[] = {
        {"LinearElastic", MaterialType::LinearElastic, false, false, false, true, false},
        {"MohrCoulomb", MaterialType::MohrCoulomb, true, false, false, true, true},
        {"HardeningSoil", MaterialType::HardeningSoil, true, true, false, false, false},
        {"HSsmall", MaterialType::HardeningSoil, true, true, false, false, false},
        {"SoftSoil", MaterialType::SoftSoil, true, false, true, false, false},
        {"SoftSoilCreep", MaterialType::SoftSoilCreep, true, false, true, false, false},
    };
    for (const auto& e : expect) {
        const ModelEntry* entry = find_model(e.name);
        check(entry != nullptr, "shipped model resolves by name");
        if (!entry) continue;
        check(entry->type == e.type, "entry carries the integrate_point tag");
        check(entry->nonlinear == e.nl, "nonlinear flag matches the driver's former switch");
        check(entry->hardening_family == e.hard, "hardening family flag");
        check(entry->softsoil_family == e.soft, "soft-soil family flag");
        check(entry->profile_E == e.pe, "E'_inc gradient flag matches the former build_profiles gate");
        check(entry->profile_c == e.pc, "c'_inc gradient flag matches the former build_profiles gate");
    }
    check(find_model("CamClay") == nullptr, "an unknown name resolves to nothing");
    check(find_model("") == nullptr, "the empty name resolves to nothing");
    check(model_names().size() >= 6, "the six built-ins are in the name list");
}

void test_le_mc_builds() {
    std::printf("build: LinearElastic and MohrCoulomb == the former inline construction\n");
    MaterialParams p = sample_params();

    MaterialModel exp = oracle_common(p);
    exp.type = MaterialType::LinearElastic;
    check(same_common(find_model("LinearElastic")->build(p), exp), "LE drained, field-exact");

    // Negative sigma_t clamps to zero exactly as the driver clamped it.
    MaterialParams pneg = p;
    pneg.tensile_strength = -5.0;
    check(find_model("LinearElastic")->build(pneg).tensile_strength == 0.0,
          "negative tensile strength clamps to 0");

    // Undrained (A): effective envelope untouched, pore-fluid flag + nu_u wired.
    MaterialParams pa = p;
    pa.drainage = DrainageClass::UndrainedA;
    MaterialModel mca = find_model("MohrCoulomb")->build(pa);
    MaterialModel expa = oracle_common(pa);
    expa.type = MaterialType::MohrCoulomb;
    check(same_common(mca, expa), "MC Undrained (A), field-exact");
    check(mca.undrained && mca.undrained_poisson == 0.495, "Undrained (A) wires Kw/n machinery");
    check(mca.friction_angle == p.phi_rad, "Undrained (A) keeps the effective envelope");

    // Undrained (B): c = su entered directly, Tresca phi = psi = 0 -- MC only.
    MaterialParams pb = p;
    pb.drainage = DrainageClass::UndrainedB;
    MaterialModel mcb = find_model("MohrCoulomb")->build(pb);
    check(mcb.undrained && mcb.friction_angle == 0.0 && mcb.dilatancy_angle == 0.0,
          "Undrained (B) is a Tresca su envelope");
    check(mcb.cohesion == p.c, "Undrained (B) keeps c = su");
}

void test_hs_builds() {
    std::printf("build: HardeningSoil / HSsmall == the former construction incl. cap calibration\n");
    MaterialParams p = sample_params();

    // Oracle: the driver's former HS branch, including the closed-form cap
    // calibration against Jaky's K0^NC.
    katai::core::HardeningSoilParams h;
    h.E50_ref = p.E50_ref; h.Eur_ref = p.Eur_ref; h.Eoed_ref = p.Eoed_ref;
    h.m = p.m; h.p_ref = p.p_ref; h.cohesion = p.c; h.friction = p.phi_rad;
    h.dilatancy = p.psi_rad; h.Rf = p.Rf; h.nu_ur = p.nu_ur;
    katai::core::HardeningSoilParams hss = h;
    hss.G0_ref = p.G0_ref; hss.gamma07 = p.gamma07;
    const double k0nc_auto = 1.0 - std::sin(p.phi_rad);
    katai::core::hs_calibrate_cap(h, k0nc_auto);
    katai::core::hs_calibrate_cap(hss, k0nc_auto);

    const MaterialModel hs = find_model("HardeningSoil")->build(p);
    check(hs.type == MaterialType::HardeningSoil, "HS tag");
    check(hs.hs.E50_ref == h.E50_ref && hs.hs.Eur_ref == h.Eur_ref &&
          hs.hs.Eoed_ref == h.Eoed_ref && hs.hs.m == h.m && hs.hs.p_ref == h.p_ref &&
          hs.hs.cohesion == h.cohesion && hs.hs.friction == h.friction &&
          hs.hs.dilatancy == h.dilatancy && hs.hs.Rf == h.Rf && hs.hs.nu_ur == h.nu_ur,
          "HS stiffness/strength wiring, field-exact");
    check(hs.hs.cap_alpha == h.cap_alpha && hs.hs.cap_beta == h.cap_beta,
          "HS cap calibration bit-identical to a direct hs_calibrate_cap call");
    check(hs.hs.G0_ref == 0.0, "plain HS carries no small-strain overlay");

    const MaterialModel sm = find_model("HSsmall")->build(p);
    check(sm.type == MaterialType::HardeningSoil, "HSsmall shares the HS tag");
    check(sm.hs.G0_ref == p.G0_ref && sm.hs.gamma07 == p.gamma07,
          "HSsmall small-strain overlay wired");
    check(sm.hs.cap_alpha == hss.cap_alpha && sm.hs.cap_beta == hss.cap_beta,
          "HSsmall cap calibration bit-identical");

    // A pinned K0^NC replaces Jaky exactly as the driver's ternary did.
    MaterialParams ppin = p;
    ppin.k0nc_auto = false; ppin.k0nc = 0.6;
    katai::core::HardeningSoilParams hpin = h;
    katai::core::hs_calibrate_cap(hpin, 0.6);
    const MaterialModel hs2 = find_model("HardeningSoil")->build(ppin);
    check(hs2.hs.cap_alpha == hpin.cap_alpha && hs2.hs.cap_beta == hpin.cap_beta,
          "pinned K0^NC feeds the calibration unchanged");
}

void test_ss_builds() {
    std::printf("build: SoftSoil / SoftSoilCreep == the former construction\n");
    MaterialParams p = sample_params();
    const double k0nc = 1.0 - std::sin(p.phi_rad);

    const MaterialModel ss = find_model("SoftSoil")->build(p);
    check(ss.type == MaterialType::SoftSoil, "SS tag");
    check(ss.ssoil.lam_star == p.lam_star && ss.ssoil.kap_star == p.kap_star &&
          ss.ssoil.nu_ur == p.nu_ur && ss.ssoil.c == p.c && ss.ssoil.phi == p.phi_rad &&
          ss.ssoil.psi == p.psi_rad && ss.ssoil.K0nc == k0nc,
          "SS parameter wiring, field-exact");

    const MaterialModel sc = find_model("SoftSoilCreep")->build(p);
    check(sc.type == MaterialType::SoftSoilCreep, "SSC tag");
    check(sc.ssc.lam_star == p.lam_star && sc.ssc.kap_star == p.kap_star &&
          sc.ssc.mu_star == p.mu_star && sc.ssc.nu_ur == p.nu_ur && sc.ssc.c == p.c &&
          sc.ssc.phi == p.phi_rad && sc.ssc.psi == p.psi_rad && sc.ssc.K0nc == k0nc,
          "SSC parameter wiring incl. mu*, field-exact");
}

void test_refusals() {
    std::printf("validate: honest refusals, exact user-facing messages\n");
    MaterialParams p = sample_params();

    check(find_model("LinearElastic")->validate(p).empty(), "LE drained accepted");
    check(find_model("MohrCoulomb")->validate(p).empty(), "MC drained accepted");
    p.drainage = DrainageClass::UndrainedB;
    check(find_model("MohrCoulomb")->validate(p).empty(), "MC Undrained (B) accepted");

    // HS refuses Undrained (B) -- and only (B) -- with the audited message.
    const std::string hs_msg =
        "Hardening Soil with Undrained (B) is not supported yet (the su/Tresca "
        "override is only implemented for Mohr-Coulomb; HS would silently follow "
        "its effective envelope like Undrained (A)). Use Undrained (A) with "
        "effective strength, or a Mohr-Coulomb set with c = su for Undrained (B).";
    check(find_model("HardeningSoil")->validate(p) == hs_msg, "HS + Undrained (B) refusal text");
    check(find_model("HSsmall")->validate(p) == hs_msg, "HSsmall shares the refusal");
    p.drainage = DrainageClass::UndrainedA;
    check(find_model("HardeningSoil")->validate(p).empty(), "HS + Undrained (A) accepted");

    // SS/SSC refuse both undrained classes with the audited message.
    const std::string ss_msg =
        "Soft Soil (Creep) with Undrained (A/B) drainage is not supported yet: the "
        "pore-fluid stiffness plumbing derives Kw/n from a constant E, and Soft Soil's "
        "stiffness is stress-dependent (K = p'/kappa*) -- it would silently behave "
        "drained. Use Drained (or model the short term with a Consolidation phase).";
    check(find_model("SoftSoil")->validate(p) == ss_msg, "SS + Undrained (A) refusal text");
    p.drainage = DrainageClass::UndrainedB;
    check(find_model("SoftSoilCreep")->validate(p) == ss_msg, "SSC + Undrained (B) refusal");
    p.drainage = DrainageClass::Drained;
    check(find_model("SoftSoil")->validate(p).empty(), "SS drained accepted");
    check(find_model("SoftSoilCreep")->validate(p).empty(), "SSC drained accepted");
}

std::string accept_all(const MaterialParams&) { return {}; }
MaterialModel build_stub(const MaterialParams& p) {
    MaterialModel mm;
    mm.type = MaterialType::LinearElastic;
    mm.youngs_modulus = 2.0 * p.E;   // recognisably not a shipped build
    return mm;
}

void test_registration() {
    std::printf("register: the external seam admits new names, never displaces shipped ones\n");
    const ModelEntry* mc_before = find_model("MohrCoulomb");

    ModelEntry user;
    user.name = "UserVerified";
    user.type = MaterialType::LinearElastic;
    user.nonlinear = false;
    user.hardening_family = false;
    user.softsoil_family = false;
    user.validate = accept_all;
    user.build = build_stub;
    check(register_model(user), "a new name registers");
    const ModelEntry* got = find_model("UserVerified");
    check(got && got->build(sample_params()).youngs_modulus == 2.0 * sample_params().E,
          "the registered model resolves and builds");

    ModelEntry imposter = user;
    imposter.name = "MohrCoulomb";
    check(!register_model(imposter), "a shipped name cannot be taken over");
    check(find_model("MohrCoulomb") == mc_before,
          "registration never moves an existing entry (pointer-stable)");
    check(find_model("MohrCoulomb")->build(sample_params()).type == MaterialType::MohrCoulomb,
          "the shipped MohrCoulomb still builds a MohrCoulomb");

    ModelEntry nameless = user;
    nameless.name = "";
    check(!register_model(nameless), "an empty name is refused");
    ModelEntry buildless = user;
    buildless.name = "Buildless";
    buildless.build = nullptr;
    check(!register_model(buildless), "an entry without a build function is refused");
}

} // namespace

int main() {
    std::printf("test_material_registry -- the constitutive construction seam\n");
    test_catalogue();
    test_le_mc_builds();
    test_hs_builds();
    test_ss_builds();
    test_refusals();
    test_registration();
    if (g_failures == 0) std::printf("all material-registry checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
