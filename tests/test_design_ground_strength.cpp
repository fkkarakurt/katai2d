// A material-factored design approach has to reach every strength the ground has -- and a joint has to
// take the strength its soil is actually solved with.
//
// EC7 DA1-C2 and DA3 divide the characteristic shear strength by the M2 partial factors and re-solve.
// Until 2026-09 four strengths were not divided, each silently and each on the unsafe side, while the
// report said the design approach had been applied:
//
//   - every INTERFACE: the joints are built from the characteristic strength before a phase applies
//     its approach, and nothing factored them afterwards. The sliding block carried 59.720232651139085
//     kN/m under EC7 DA3 -- its characteristic capacity, bit for bit.
//   - the depth GRADIENT of cohesion, c_inc: only c_ref was divided, so the deep ground kept its
//     characteristic strength. A footing on c_inc = 20 kPa/m read 125.46669113812283 kN/m under DA3,
//     bit for bit the run with c_ref alone factored by hand.
//   - SOFT SOIL and SOFT SOIL CREEP: both keep their failure line in their own parameter blocks,
//     which the factoring never touched. A footing pushed into each read 61.458466027030674 and
//     68.35659731218506 kN/m under DA3, bit for bit the characteristic runs.
//
// And one more, in the same place, without any design approach at all: a joint beside an Undrained
// (B) or (C) clay took its friction angle from the box the clay itself ignores (the soil is a Tresca
// material there, c = su, phi = 0, and the input contract says the entered phi "is ignored for
// strength"). With the box left at 26.6 deg the sliding block's Tresca joint carried 59.72 kN/m
// instead of B su = 10.
//
// verify: KV-STR-011
//   oracle:   closed_form
//   source:   EN 1997-1:2004 Eurocode 7, Geotechnical design -- Part 1: General rules, 2.4.6.2 (the design value of a geotechnical parameter, X_d = X_k / gamma_M, Eq. 2.2) and Annex A, Table A.4 (set M2: gamma_phi' = 1.25 applied to tan(phi'), gamma_c' = 1.25, gamma_cu = 1.4); Coulomb friction with adhesion on a planar joint, the same statics as KV-STR-002, and for an undrained joint the Tresca limit tau = R_inter s_u; KATAI 2D input contract (docs/k2d-format.md: drainage, Rinter, iface_material, phases[].design)
//   locator:  the checked-in KV-STR-002 block (B = 4 m, W = 100 kN/m, joint data set with R_inter = 1) pushed sideways until it slides, stated in full: (a) under EC7 DA3 and DA1-C2, F = B c_w / gamma_c' + W tan(phi_w) / gamma_phi' with c_w = 2.5 kN/m2, phi_w = 26.6 deg; (b) the joint data set Undrained (B), and then Undrained (C), with c = s_u = 2.5 kN/m2 and its friction box left at 26.6 deg, F = B s_u; (c) the same Undrained (B) joint under DA3, F = B s_u / gamma_cu
//   quantity: the horizontal failure force, as the sum of the reactions on the pushed edge, run from the checked-in tests/corpus/kv-str-002-sliding-block.k2d with the joint data set and the phase's design approach changed [kN/m]
//   expected: (a) 48.061016; (b) 10.0; (c) 7.142857
//   band:     2% vs the closed form, the band KV-STR-002 holds this mesh to, as asserted below -- measured -0.6627% for (a) (47.742501180826864 under DA3 and under DA1-C2), -1.0417% for (b) under both drainage types and for (c). (a) does not carry the characteristic run's bias (-0.5926%) because the adhesion and the friction terms have different discretisation errors on this mesh, which KV-STR-002 shows term by term; (c) has the adhesion term alone, and reproduces the characteristic undrained run's bias to 1e-9 relative. Two identities besides, each to 1e-9 relative: every design run equals the characteristic run on hand-factored joint strength (measured bit for bit), and an undrained joint gives the same force whatever its ignored friction box holds (bit for bit). Before 2026-09 (a) read 59.720232651139085, the characteristic capacity, and (b) read 59.72 against B s_u = 10
//
// verify: KV-FND-015
//   oracle:   independent_path
//   source:   EN 1997-1:2004, 2.4.6.2 (X_d = X_k / gamma_M, Eq. 2.2) and Annex A, Table A.4 (M2); the depth-varying cohesion of the input contract, c(y) = c_ref + c_inc (y_ref - y) (docs/k2d-format.md, c_inc / y_ref); the Soft Soil and Soft Soil Creep strength parameters c' and phi' (docs/k2d-format.md)
//   locator:  a strip footing pushed 0.05 m (Mohr-Coulomb) or 0.15 m (Soft Soil, Soft Soil Creep) into a 6 m x 4 m block, the footing reaction under EC7 DA3 against the SAME model run with no design approach on strength factored by hand: c_ref / 1.25, c_inc / 1.25, phi = atan(tan(phi) / 1.25); for an Undrained (B) Tresca material s_u / 1.4 and c_inc / 1.4
//   quantity: the vertical reaction under the footing [kN/m]
//   expected: the hand-factored run
//   band:     1e-9 relative, as asserted below -- the two runs solve the same numbers by two routes, measured 0, bit for bit, in all four runs here; the same comparison with the hand-factored angle carried through a differently rounded degree-radian conversion read 2.2e-11 for Soft Soil and 1.8e-13 for Soft Soil Creep, which is the scale the band leaves room for. The soft-soil pair is run with the design approach in BOTH phases: applied only in the loading phase the Soft Soil pair differs by 1.06e-5, because the preconsolidation that the first phase seeds is computed from the strength that phase was given (characteristic in one run, factored in the other), which is a different initial state and not a factoring error. Each case also shows that its design run is NOT its characteristic run, by more than the band: a factoring that reached nothing would pass the identity only if the characteristic and hand-factored runs agreed

#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kGammaPhi = 1.25, kGammaC = 1.25, kGammaCu = 1.4;   // EN 1997-1 Table A.4, M2

double factored_phi_deg(double phi_deg) {
    return std::atan(std::tan(phi_deg * kPi / 180.0) / kGammaPhi) * 180.0 / kPi;
}

double rel(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

m::Project corpus(const char* file) {
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/" + file;
    if (!m::load_project(path, pr, &err, nullptr)) check(false, "load " + path + ": " + err);
    return pr;
}

// ------------------------------------------------------------------ the joint (KV-STR-011) --
// The checked-in KV-STR-002 block; materials[1] is the joint's own data set.
m::Project block(m::DesignApproach design, m::Drainage drainage, double c, double phi) {
    m::Project pr = corpus("kv-str-002-sliding-block.k2d");
    if (pr.phases.empty() || pr.materials.size() < 2) return pr;
    pr.phases[0].design_approach = design;
    m::Material& joint = pr.materials[1];
    joint.drainage = drainage;
    joint.c = c;
    joint.phi = phi;
    if (drainage == m::Drainage::UndrainedC) { joint.nu = 0.495; joint.nu_u = 0.499; }
    return pr;
}

double sliding_force(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return -1.0;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return -1.0;
    }
    const auto& R = res[1];
    double rx = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.x[n]) < 1e-9) rx += R.reaction[2 * n];
    return std::fabs(rx);
}

void joint_takes_the_design_strength() {
    std::printf("\nKV-STR-011: the joint under a design approach, and beside an undrained clay\n");
    constexpr double B = 4.0, W = 100.0, cw = 2.5, phiw = 26.6;
    const double tanw = std::tan(phiw * kPi / 180.0);
    using DA = m::DesignApproach;
    using Dr = m::Drainage;

    // The characteristic run fixes the case's own discretisation bias for the drained joint.
    const double F_k = sliding_force(block(DA::None, Dr::Drained, cw, phiw));
    const double bias_drained = (F_k - (B * cw + W * tanw)) / (B * cw + W * tanw);

    // (a) EC7 DA3 and DA1-C2 factor the joint as they factor the soil.
    const double exact_a = B * cw / kGammaC + W * tanw / kGammaPhi;
    const double F_da3 = sliding_force(block(DA::EC7_DA3, Dr::Drained, cw, phiw));
    const double F_c2 = sliding_force(block(DA::EC7_DA1_C2, Dr::Drained, cw, phiw));
    const double F_hand = sliding_force(block(DA::None, Dr::Drained, cw / kGammaC, factored_phi_deg(phiw)));
    std::printf("      characteristic %.12f (%+.4f%%) | DA3 %.12f | DA1-C2 %.12f | hand-factored joint "
                "%.12f | design closed form %.6f\n",
                F_k, 100.0 * bias_drained, F_da3, F_c2, F_hand, exact_a);
    // The bias of the drained joint moves with its two terms (KV-STR-002 shows them apart), so the
    // design run is held to the closed form through the identity with hand-factored input, and to
    // the closed form itself within the case's own 2%.
    check(F_da3 > 0.0 && std::fabs(F_da3 - exact_a) < 0.02 * exact_a,
          "DA3: the failure force is the design closed form, within KV-STR-002's 2%");
    check(F_da3 > 0.0 && rel(F_da3, F_hand) <= 1e-9,
          "DA3 equals the characteristic run on hand-factored joint strength");
    check(F_c2 > 0.0 && rel(F_c2, F_hand) <= 1e-9, "and so does DA1-C2 (the same M2 set)");
    check(F_da3 < 0.9 * F_k, "and the design run is not the characteristic one");

    // (b) An undrained joint is a Tresca joint, whatever the friction box of its data set holds.
    const double exact_b = B * cw;
    for (const auto& [drain, name] : {std::pair{Dr::UndrainedB, "Undrained (B)"},
                                      std::pair{Dr::UndrainedC, "Undrained (C)"}}) {
        const double F_box = sliding_force(block(DA::None, drain, cw, phiw));
        const double F_zero = sliding_force(block(DA::None, drain, cw, 0.0));
        std::printf("      %s joint, s_u = %.1f: friction box 26.6 deg %.12f | box 0 %.12f | B s_u "
                    "%.1f (%+.4f%%)\n",
                    name, cw, F_box, F_zero, exact_b, 100.0 * (F_box - exact_b) / exact_b);
        check(F_box > 0.0 && std::fabs(F_box - exact_b) < 0.02 * exact_b,
              std::string(name) + ": the joint carries B s_u, within KV-STR-002's 2%");
        check(F_box > 0.0 && rel(F_box, F_zero) <= 1e-9,
              std::string(name) + ": the friction box the clay ignores, the joint ignores too");
    }

    // (c) ...and under DA3 it takes gamma_cu, not gamma_c'.
    const double exact_c = B * cw / kGammaCu;
    const double F_ud = sliding_force(block(DA::EC7_DA3, Dr::UndrainedB, cw, phiw));
    const double F_ud_hand = sliding_force(block(DA::None, Dr::UndrainedB, cw / kGammaCu, 0.0));
    const double F_ud_k = sliding_force(block(DA::None, Dr::UndrainedB, cw, 0.0));
    const double bias_ud = (F_ud_k - exact_b) / exact_b;
    std::printf("      Undrained (B) joint under DA3: %.12f | hand-factored s_u / 1.4 %.12f | B s_u / "
                "gamma_cu %.6f (%+.4f%%, characteristic bias %+.4f%%)\n",
                F_ud, F_ud_hand, exact_c, 100.0 * (F_ud - exact_c) / exact_c, 100.0 * bias_ud);
    check(F_ud > 0.0 && std::fabs((F_ud - exact_c) / exact_c - bias_ud) <= 1e-9,
          "DA3 on the undrained joint: B s_u / gamma_cu, with the characteristic run's exact bias");
    check(F_ud > 0.0 && rel(F_ud, F_ud_hand) <= 1e-9,
          "and it equals the run on s_u / 1.4 entered by hand");
}

// ------------------------------------------------------ gradient and soft soils (KV-FND-015) --
// The checked-in Prandtl footing geometry (6 m x 4 m, footing 2.4..3.6 m), pushed by a prescribed
// settlement instead of loaded, so a design approach changes the reaction and not the load.
m::Project footing(m::SoilModel model, m::Drainage drainage, double c, double c_inc, double phi,
                   double push, m::DesignApproach design, bool design_initial) {
    m::Project pr = corpus("kv-fnd-010-prandtl-strip-footing.k2d");
    if (pr.phases.empty() || pr.materials.empty()) return pr;
    pr.mesh.order = 6;
    pr.mesh.elem_size = 0.4;
    m::Material& s = pr.materials[0];
    s.model = model;
    s.drainage = drainage;
    s.c = c; s.c_inc = c_inc; s.y_ref = 4.0; s.phi = phi; s.psi = 0.0;
    if (model == m::SoilModel::SoftSoil || model == m::SoilModel::SoftSoilCreep) {
        // The soft soils need a stress level to be stiff at: weight, a pinned K0 and K0nc (so that
        // factoring phi cannot move either), and the geostatic seed.
        s.gamma_unsat = s.gamma_sat = 16.0;
        s.lam_star = 0.05; s.kap_star = 0.01; s.nu_ur = 0.15; s.mu_star = 0.002;
        s.k0_auto = false; s.k0 = 0.6;
        s.k0nc_auto = false; s.k0nc = 0.6;
        pr.initial_procedure = m::InitialProcedure::K0Procedure;
    }
    pr.loads.clear();
    m::PrescribedDisp d;
    d.name = "Footing";
    d.x1 = 2.4; d.y1 = 4.0; d.x2 = 3.6; d.y2 = 4.0;
    d.set_ux = false; d.set_uy = true; d.uy = -push;
    pr.disps = {d};
    pr.initial.load_active.clear();
    pr.initial.disp_active = {0};
    pr.initial.design_approach = design_initial ? design : m::DesignApproach::None;
    pr.phases.resize(1);
    pr.phases[0].load_active.clear();
    pr.phases[0].disp_active = {1};
    pr.phases[0].design_approach = design;
    return pr;
}

double footing_reaction(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return 0.0;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return 0.0;
    }
    const auto& R = res[1];
    double ry = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.y[n] - 4.0) < 1e-9 && R.mesh.x[n] > 2.4 - 1e-9 && R.mesh.x[n] < 3.6 + 1e-9)
            ry += R.reaction[2 * n + 1];
    return std::fabs(ry);
}

void identity(const char* what, const m::Project& design, const m::Project& hand,
              const m::Project& characteristic) {
    const double Fd = footing_reaction(design), Fh = footing_reaction(hand);
    const double Fk = footing_reaction(characteristic);
    std::printf("      %-38s design %.14f | hand-factored %.14f (%.1e) | characteristic %.6f\n", what,
                Fd, Fh, Fh > 0.0 ? rel(Fd, Fh) : 0.0, Fk);
    check(Fd > 0.0 && Fh > 0.0 && rel(Fd, Fh) <= 1e-9,
          std::string(what) + ": EC7 DA3 equals the run on strength factored by hand");
    check(Fk > 0.0 && rel(Fd, Fk) > 1e-4,
          std::string(what) + ": and the design run is not the characteristic one");
}

void gradient_and_soft_soils_take_the_design_strength() {
    std::printf("\nKV-FND-015: the cohesion gradient and the soft soils under EC7 DA3\n");
    using SM = m::SoilModel;
    using Dr = m::Drainage;
    constexpr auto DA3 = m::DesignApproach::EC7_DA3, NONE = m::DesignApproach::None;
    const double phi = 10.0, phid = factored_phi_deg(phi);

    identity("Mohr-Coulomb, c_inc = 20 kPa/m",
             footing(SM::MohrCoulomb, Dr::Drained, 10.0, 20.0, phi, 0.05, DA3, false),
             footing(SM::MohrCoulomb, Dr::Drained, 10.0 / kGammaC, 20.0 / kGammaC, phid, 0.05, NONE, false),
             footing(SM::MohrCoulomb, Dr::Drained, 10.0, 20.0, phi, 0.05, NONE, false));
    identity("Undrained (B) clay, c_inc = 20 kPa/m",
             footing(SM::MohrCoulomb, Dr::UndrainedB, 10.0, 20.0, 0.0, 0.05, DA3, false),
             footing(SM::MohrCoulomb, Dr::UndrainedB, 10.0 / kGammaCu, 20.0 / kGammaCu, 0.0, 0.05, NONE, false),
             footing(SM::MohrCoulomb, Dr::UndrainedB, 10.0, 20.0, 0.0, 0.05, NONE, false));
    const double phi_ss = 25.0, phid_ss = factored_phi_deg(phi_ss);
    identity("Soft Soil",
             footing(SM::SoftSoil, Dr::Drained, 5.0, 0.0, phi_ss, 0.15, DA3, true),
             footing(SM::SoftSoil, Dr::Drained, 5.0 / kGammaC, 0.0, phid_ss, 0.15, NONE, false),
             footing(SM::SoftSoil, Dr::Drained, 5.0, 0.0, phi_ss, 0.15, NONE, false));
    identity("Soft Soil Creep",
             footing(SM::SoftSoilCreep, Dr::Drained, 5.0, 0.0, phi_ss, 0.15, DA3, true),
             footing(SM::SoftSoilCreep, Dr::Drained, 5.0 / kGammaC, 0.0, phid_ss, 0.15, NONE, false),
             footing(SM::SoftSoilCreep, Dr::Drained, 5.0, 0.0, phi_ss, 0.15, NONE, false));
}

}  // namespace

int main() {
    std::printf("A design approach reaches every ground strength; a joint takes its soil's\n");
    joint_takes_the_design_strength();
    gradient_and_soft_soils_take_the_design_strength();
    if (g_failures == 0) {
        std::printf("\nOK: joints, gradients and soft soils carry the design strength\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
