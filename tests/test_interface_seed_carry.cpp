// THE NORMAL STRESS A WALL'S JOINT STARTS FROM, AND WHAT HAPPENS TO IT AFTERWARDS.
//
// A joint is a constitutive element: its normal stress is the stress it was installed with plus
// what its own opening and closing add, sigma_n = sigma_n0 + k_n du_n. Nothing else may change it
// -- no more than the soil beside it, whose stress changes only through its strain. Two rules
// follow, and each is checked here against a closed form or an identity, not a tolerance band:
//
// (A) THE SEED OF A WALL ACTIVE FROM THE START is the geostatic stress of its own side,
//     sigma_n0 = K0 sigma'_v with sigma'_v = gamma (surface of THAT side - y). A wall standing in a
//     step of the ground surface retains 4 m more on its high side than on its low one. Read on
//     the wall's line itself the overburden was the higher surface for both faces (the surface
//     at x is the highest top of the regions that reach x), and with the step drawn the other way
//     the half-open rule could pick either: the low face pressed as hard as the high one.
//     A step is not at rest under K0, so the K0 phase moves the ground a little and every station
//     carries k_n du_n besides its seed. The check is therefore the joint's own law: at the toe,
//     where the wall is pinned to the ground (du_n = 0), sigma_n IS the seed; everywhere else
//     (sigma_n - seed of its own side) / du_n is ONE stiffness k_n for every station of both faces
//     (the virtual thickness is a tenth of the model's element size, the same everywhere, and the
//     soil is one material). With the other side's seed the ratio scatters by more than 100%.
//     Both, with the step drawn on either side of the wall.
//
// (B) A JOINT'S STRESS IS CARRIED FROM PHASE TO PHASE. Until 2026-10-23 the phase chain carried a
//     joint's slip but rebuilt its sigma_n0 in every phase from the K0 formula on that phase's
//     geometry. The difference to the stress the parent phase was in equilibrium with went into
//     the phase's baseline and acted as a self-equilibrated load between wall and ground: an
//     excavation lowered the pit face's sigma_n0 by K0 gamma (depth dug) at once, without the
//     ground beside it having moved. The identity that exposes it is the definition of a nil
//     phase: a wall installed after a surcharge (seeded from the stress the ground carries then,
//     which includes the surcharge's lateral share) and followed by a phase that changes nothing
//     must not move, and its joint stresses must not change.
//
// verify: KV-STR-021
//   oracle:   closed_form
//   source:   Jaky (1944) / at-rest earth pressure sigma'_h = K0 sigma'_v with sigma'_v the overburden of the soil column above the point; KATAI 2D input contract (docs/k2d-format.md, structs[i].iface_pos / iface_neg): a wall's interfaces start from the at-rest stress of the ground beside them; the definition of a nil phase (a phase whose configuration equals its parent's applies nothing and changes nothing)
//   locator:  (A) a 1..6 m wall in a step of the ground surface (10 m on one side, 6 m on the other), linear elastic soil gamma = 18, K0 = 0.5, no water, active in the initial K0 phase, the step drawn on either side; (B) a wall installed in a Plastic phase after a one-sided 50 kPa surcharge, followed by a phase identical to it
//   quantity: (A) the joint normal stress of each face [kPa] at the toe, and the stiffness (sigma_n - seed) / du_n along it [kPa/m]; (B) the nil phase's displacement [m] and its change of the wall moment [kNm/m] and of the joint stresses [kPa]
//   expected: (A) at the toe 0.5 * 18 * (10 - 1) = 81 on the high face and 0.5 * 18 * (6 - 1) = 45 on the low face; one k_n for every station of both faces; (B) zero displacement and unchanged forces and joint stresses
//   band:     (A) 1e-9 relative at the toe (the formula itself); k_n within 1% across the stations (measured 7e-14: the identity holds to round-off), against a scatter above 100% with the other side's seed; (B) 1e-6 of the parent's values, as the existing nil-phase identity (round-off only)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

m::Material soil() {
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
    s.c = 5.0; s.phi = 30.0; s.rinter_rigid = false; s.Rinter = 1.0;
    s.k0_auto = false; s.k0 = 0.5;
    return s;
}

m::PlateMaterial wall_section() {
    m::PlateMaterial pm;
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = 0.0;
    return pm;
}

// (A) One face at a time (iface_pos or iface_neg), so each run's stations are one face's.
struct Station { double y, sn, gap; };
std::vector<Station> step_face(bool high_on_left, bool positive_face) {
    m::Project pr;
    pr.materials.push_back(soil());
    pr.has_water = false;
    const double hl = high_on_left ? 10.0 : 6.0, hr = high_on_left ? 6.0 : 10.0;
    m::SoilPolygon L; L.material = 0; L.name = "Left";
    L.x = {0, 10, 10, 0}; L.y = {0, 0, hl, hl};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free, (int)m::BCType::Free,
                 (int)m::BCType::HorizontallyFixed};
    m::SoilPolygon R; R.material = 0; R.name = "Right";
    R.x = {10, 20, 20, 10}; R.y = {0, 0, hr, hr};
    R.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::Free};
    pr.polygons = {L, R};
    pr.plates.push_back(wall_section());
    m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
    st.x1 = 10; st.y1 = 1; st.x2 = 10; st.y2 = 6; st.material = 0;
    st.iface_pos = positive_face; st.iface_neg = !positive_face;
    pr.structs.push_back(st);

    std::vector<Station> out;
    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    if (!M.ok) return out;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.empty() || !res[0].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res[0].message.c_str());
        return out;
    }
    for (const auto& ir : res[0].interface_forces)
        for (const auto& k : ir.stations) out.push_back({k.y, k.sigma_n, k.gap});
    return out;
}

void case_step() {
    std::printf("-- (A) a wall in a step of the ground: each face its own overburden --\n");
    const auto seed = [](double surface, double y) { return -0.5 * 18.0 * (surface - y); };
    for (bool high_left : {true, false}) {
        const auto a = step_face(high_left, true), b = step_face(high_left, false);
        check(!a.empty() && !b.empty(), "both faces solved with their stations");
        if (a.empty() || b.empty()) continue;
        // Which face is the high one: the one whose toe presses harder.
        const auto toe = [](const std::vector<Station>& f) {
            const Station* t = &f[0];
            for (const auto& k : f) if (k.y < t->y) t = &k;
            return *t;
        };
        const bool a_high = std::fabs(toe(a).sn) > std::fabs(toe(b).sn);
        const auto& hi = a_high ? a : b;
        const auto& lo = a_high ? b : a;
        const double toe_hi = std::fabs(toe(hi).sn - seed(10.0, 1.0)) / 81.0;
        const double toe_lo = std::fabs(toe(lo).sn - seed(6.0, 1.0)) / 45.0;
        // One k_n for every station that opened or closed, each read against a seed.
        const auto spread = [&](double s_hi, double s_lo) {
            double kmin = 1e300, kmax = -1e300, sum = 0.0;
            int n = 0;
            for (int f = 0; f < 2; ++f)
                for (const auto& k : f ? lo : hi) {
                    if (std::fabs(k.gap) < 1e-6) continue;
                    const double kn = (k.sn - seed(f ? s_lo : s_hi, k.y)) / k.gap;
                    kmin = std::fmin(kmin, kn); kmax = std::fmax(kmax, kn); sum += kn; ++n;
                }
            return n > 0 ? (kmax - kmin) / std::fabs(sum / n) : 1e300;
        };
        const double own = spread(10.0, 6.0), other = spread(6.0, 10.0);
        std::printf("   high side %s: toe %.6g / %.6g kPa (closed form 81 / 45); k_n spread %.2e "
                    "with each face's own seed, %.2e with the other's\n", high_left ? "left " : "right",
                    std::fabs(toe(hi).sn), std::fabs(toe(lo).sn), own, other);
        check(toe_hi < 1e-9 && toe_lo < 1e-9, "at the toe each face carries K0 gamma (its own surface - y)");
        check(own < 0.01, "along the wall: sigma_n - own seed = ONE k_n times the opening, both faces");
        check(other > 1.0, "the other side's seed would scatter it (the check has teeth)");
    }
}

// (B) Installed after the surcharge, then a nil phase.
void case_installed_then_nil() {
    std::printf("-- (B) a joint keeps its stress: a nil phase after a wall installed on loaded "
                "ground --\n");
    m::Project pr;
    pr.materials.push_back(soil());
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(P);
    pr.plates.push_back(wall_section());
    m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
    st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
    st.iface_pos = true; st.iface_neg = true;
    pr.structs.push_back(st);
    m::Load Ld; Ld.kind = m::LoadKind::Distributed;
    Ld.x1 = 0; Ld.y1 = 10; Ld.x2 = 8; Ld.y2 = 10;
    Ld.qx1 = 0; Ld.qy1 = -50; Ld.qx2 = 0; Ld.qy2 = -50;
    pr.loads.push_back(Ld);
    pr.initial.struct_active = {0};
    pr.initial.load_active = {0};
    m::Phase p1; p1.type = m::PhaseType::Plastic; p1.name = "Surcharge";
    p1.struct_active = {0}; p1.load_active = {1};
    m::Phase p2; p2.type = m::PhaseType::Plastic; p2.name = "Install";
    p2.struct_active = {1}; p2.load_active = {1};
    m::Phase p3 = p2; p3.name = "Nil";
    pr.phases = {p1, p2, p3};

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "model meshed");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 4 || !res[1].ok || !res[2].ok || !res[3].ok) {
        check(false, "K0 + Surcharge + Install + Nil all solved");
        for (const auto& r : res) std::printf("  (%s)\n", r.message.c_str());
        return;
    }
    const auto& par = res[2];
    const auto& nil = res[3];
    double difM = 0.0, teethM = 0.0;
    for (size_t i = 0; i < nil.struct_forces.size() && i < par.struct_forces.size(); ++i)
        for (size_t k = 0; k < nil.struct_forces[i].stations.size() &&
                           k < par.struct_forces[i].stations.size(); ++k) {
            difM = std::fmax(difM, std::fabs(nil.struct_forces[i].stations[k].M -
                                             par.struct_forces[i].stations[k].M));
            teethM = std::fmax(teethM, std::fabs(par.struct_forces[i].stations[k].M));
        }
    double difS = 0.0, teethS = 0.0;
    for (size_t i = 0; i < nil.interface_forces.size() && i < par.interface_forces.size(); ++i)
        for (size_t k = 0; k < nil.interface_forces[i].stations.size() &&
                           k < par.interface_forces[i].stations.size(); ++k) {
            const auto& a = nil.interface_forces[i].stations[k];
            const auto& b = par.interface_forces[i].stations[k];
            difS = std::fmax(difS, std::fmax(std::fabs(a.sigma_n - b.sigma_n), std::fabs(a.tau - b.tau)));
            teethS = std::fmax(teethS, std::fabs(b.sigma_n));
        }
    std::printf("   surcharge |u| %.4g m; nil |u| %.3g m; |dM| %.4g of %.4g kNm/m; joint |d| %.4g "
                "of %.4g kPa\n", res[1].max_disp, nil.max_disp, difM, teethM, difS, teethS);
    check(res[1].max_disp > 1e-4 && teethS > 10.0, "the surcharge deforms the ground and the joint carries stress (teeth)");
    check(nil.max_disp < 1e-6 * res[1].max_disp, "NIL phase: nothing moves");
    check(difM <= 1e-6 * std::fmax(teethM, 1.0), "NIL phase: the wall's moment is the parent's");
    check(difS < 1e-6 * teethS, "NIL phase: the joint's stresses are the parent's");
}
}  // namespace

int main() {
    std::printf("A wall's joint: the stress it starts from and the stress it keeps\n\n");
    case_step();
    case_installed_then_nil();
    if (g_failures == 0) {
        std::printf("\nOK: each face starts from its own side's at-rest stress and keeps it\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", g_failures);
    return 1;
}
