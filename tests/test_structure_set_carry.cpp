// A phase that activates or removes a structure must carry every other structure's state.
//
// The static chain carries the parent's structural state -- the total displacement datum and the
// committed plastic state -- into the next phase (test_staged_struct_carry pins the nil identity).
// Until 2026-09 it did so only when the structure set was IDENTICAL: activating one anchor dropped
// the whole carry, every other structure's force restarted from what the new phase alone put into
// it, and the only statement was a sentence in the phase message. Measured on the wall of
// test_staged_struct_carry, in a phase whose only change was an anchor of EA = 1 kN far from the
// wall: the wall's peak moment went from 1.914704 to 0.921088 kNm/m (-52%) with interfaces and from
// 17.889413 to 10.867339 (-39%) without. That is the standard anchored-excavation sequence.
//
// The carry is now matched by drawn structure, and a structure activated in a phase is installed on
// the ground as the parent left it: it reads the displacement since its installation. On a LINEAR
// model both rules have an exact independent path, superposition, because every phase is then a
// linear solve of its own change:
//
//   (a) activating a structure in a phase that changes nothing else changes nothing: no displacement,
//       the wall's moment unchanged, the new anchor's force zero;
//   (b) loading after that installation: the wall's moment is the first stage's (wall alone, first
//       load) plus the second stage's (wall and anchor, second load alone), and the anchor's force is
//       the second stage's alone -- it did not exist during the first;
//   (c) removing the anchor again: the change is the soil-and-wall model loaded by the force the anchor
//       carried, applied at its node;
//   (d) a PRESTRESSED anchor activated in the chain applies its lock-off force in that phase: the same
//       chain with a slack anchor and an equal point force along its axis gives the same answer.
//
// verify: KV-STR-012
//   oracle:   independent_path
//   source:   the principle of superposition for a linear elastic system (each phase of a chain of linear solves is the solve of that phase's change, and the solves add); the KATAI 2D input contract for staged activation (docs/k2d-format.md, phases[i].struct) and for the anchor lock-off force (anchors[i].prestress: applied when the anchor is activated, the anchor behaving as an elastic spring from that state); the prestress equivalence of KV-STR-001 (a constant internal force N0 along the anchor axis is, in the residual, an external force of N0 along that axis at the attached node)
//   locator:  weightless linear-elastic block 20 m x 10 m (E = 30000 kPa, nu = 0.3), base fixed, sides on rollers; a plate wall x = 10 m, y = 4..10 m (EA = 1.2e7 kN/m, EI = 1.6e5 kNm2/m) active from the first phase; a fixed-end anchor from (10, 9) to (-5, 12), EA = 2e5 kN, pulling the wall back against the surcharge; surcharge q1 = 50 kPa on 0 <= x <= 8, then a further 30 kPa on the same strip. Chain: q1 (wall) -> anchor activated -> +30 kPa -> anchor removed. Independent runs: the wall alone under q1; the wall and anchor from the start under 30 kPa alone; the wall alone under the anchor's released force as a point load; and the chain with a slack anchor plus a point force of N0 = 100 kN along its axis
//   quantity: the wall's bending moment at every station [kNm/m], the anchor's axial force [kN] and the largest displacement increment of each phase [m]
//   expected: (a) increment 0, moment unchanged, anchor force 0; (b) M_chain = M_first + M_second station by station, N_chain = N_second; (c) M_after - M_before = M_released, and the displacement increment equal to the released run's; (d) the two chains equal
//   band:     1e-9 relative on moments and forces and 1e-9 m on displacements, as asserted below -- the runs are linear solves converged to 1e-10 and differ only by solver round-off: measured (a) an increment of exactly 0 and the moment unchanged bit for bit; (b) 1.6e-13 on the moments and 8.7e-13 on the anchor force (3.487162 kN); (c) 1.5e-11 on the moment change and 5.3e-16 m on the displacements; (d) 8.4e-13. The equivalent point forces act along the line from the node the anchor ATTACHES to, (10, 8.875), and not from its drawn head: in this case's first draft, with the fixed point at (25, 12), the equivalences written from the drawn head missed by 1.6e-3 and 1.2e-3, and nothing but that direction separated them. The defect this case was built for is several orders outside the band: before the carry was matched by structure, (a) moved the wall's moment by 39% (measured on the same wall without interfaces) and (b) would read the second stage alone. A second, nonlinear check keeps the case honest where superposition does not hold: on the weighted interface wall, the phase activating an anchor of EA = 1 kN leaves the wall's moment and the joint's stresses unchanged to 1e-6 relative, the nil-phase tolerance of test_staged_struct_carry -- measured unchanged bit for bit, with and without interfaces

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;
using katai::app::SolveResult;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

constexpr double kAx = 10.0, kAy = 9.0, kFx = -5.0, kFy = 12.0;   // anchor head / fixed point
constexpr double kQ1 = 50.0, kDq = 30.0, kN0 = 100.0;

// Struct order: 0 wall, 1 anchor. Load order: 0 q1, 1 extra, 2 point load (off unless asked).
m::Project model(bool weightless_linear) {
    m::Project pr;
    pr.name = "structure set carry";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.mesh.elem_size = 1.0; pr.mesh.order = 6; pr.mesh.auto_refine = false;
    m::Material s;
    s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.e_init = 0.6; s.c = 5.0; s.phi = 30.0;
    s.gamma_unsat = weightless_linear ? 0.0 : 18.0;
    s.gamma_sat = weightless_linear ? 0.0 : 20.0;
    s.rinter_rigid = false; s.Rinter = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    m::PlateMaterial pm; pm.name = "Wall";
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = weightless_linear ? 0.0 : 24.0 * 0.4;
    pr.plates.push_back(pm);
    m::StructElement wall; wall.kind = m::StructKind::Plate; wall.name = "Wall"; wall.material = 0;
    wall.x1 = 10; wall.y1 = 4; wall.x2 = 10; wall.y2 = 10;
    pr.structs.push_back(wall);
    m::AnchorMaterial am; am.name = "Tie"; am.EA = 2.0e5; am.Lspacing = 1.0;
    pr.anchors.push_back(am);
    m::StructElement tie; tie.kind = m::StructKind::Anchor; tie.name = "Tie"; tie.material = 0;
    tie.x1 = kAx; tie.y1 = kAy; tie.x2 = kFx; tie.y2 = kFy;
    pr.structs.push_back(tie);
    m::Load q1; q1.kind = m::LoadKind::Distributed; q1.name = "q1";
    q1.x1 = 0; q1.y1 = 10; q1.x2 = 8; q1.y2 = 10; q1.qx1 = q1.qx2 = 0; q1.qy1 = q1.qy2 = -kQ1;
    pr.loads.push_back(q1);
    m::Load extra = q1; extra.name = "extra"; extra.qy1 = extra.qy2 = -kDq;
    pr.loads.push_back(extra);
    m::Load pt; pt.kind = m::LoadKind::Point; pt.name = "point";
    pt.x1 = pt.x2 = kAx; pt.y1 = pt.y2 = kAy; pt.qx1 = pt.qx2 = 0; pt.qy1 = pt.qy2 = 0;
    pr.loads.push_back(pt);
    pr.initial.struct_active = {1, 0};
    pr.initial.load_active = {0, 0, 0};
    return pr;
}

m::Phase phase(const char* name, std::vector<char> st, std::vector<char> ld) {
    m::Phase ph; ph.name = name; ph.type = m::PhaseType::Plastic;
    ph.struct_active = std::move(st); ph.load_active = std::move(ld);
    return ph;
}

std::vector<SolveResult> solve(const m::Project& pr, InitialPhase ip) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { check(false, "meshed: " + M.message); return {}; }
    auto res = katai::app::solve_phases(pr, M.mesh, ip);
    for (const auto& r : res)
        if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    return res;
}

const katai::core::StructForce* force_of(const SolveResult& r, const char* name) {
    for (const auto& f : r.struct_forces)
        if (f.name == name) return &f;
    return nullptr;
}

std::vector<double> moments(const SolveResult& r) {
    std::vector<double> M;
    if (const auto* f = force_of(r, "Wall"))
        for (const auto& st : f->stations) M.push_back(st.M);
    return M;
}

double anchor_N(const SolveResult& r) {
    const auto* f = force_of(r, "Tie");
    return f && !f->stations.empty() ? f->stations[0].N : 0.0;
}

double peak(const std::vector<double>& v) {
    double p = 0.0;
    for (double x : v) p = std::max(p, std::fabs(x));
    return p;
}

// Largest |a - b| over two equally long station lists, relative to the larger peak.
double rel_gap(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) return 1e300;
    double gap = 0.0;
    for (size_t k = 0; k < a.size(); ++k) gap = std::max(gap, std::fabs(a[k] - b[k]));
    return gap / std::max(peak(a), peak(b));
}

double max_abs_diff(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    if (a.size() != b.size()) return 1e300;
    return (a - b).cwiseAbs().maxCoeff();
}

void superposition_on_a_linear_model() {
    std::printf("\n(a)-(d) the chain against superposition, on a linear model\n");
    // The chain: q1 on the wall alone, then the anchor, then more load, then the anchor removed.
    m::Project chain = model(true);
    chain.phases = {phase("q1", {1, 0}, {1, 0, 0}), phase("anchor", {1, 1}, {1, 0, 0}),
                    phase("more", {1, 1}, {1, 1, 0}), phase("remove", {1, 0}, {1, 1, 0})};
    const auto C = solve(chain, InitialPhase::K0Procedure);
    check(C.size() == 5 && std::all_of(C.begin(), C.end(), [](const SolveResult& r) { return r.ok; }),
          "the chain solves: initial, q1, anchor, more, remove");
    if (C.size() != 5) return;

    // (a) Installing the anchor changes nothing.
    const auto M1 = moments(C[1]), M2 = moments(C[2]);
    std::printf("      q1: max|M| %.6f | anchor phase: max|M| %.6f, increment %.3e m, N %.3e kN\n",
                peak(M1), peak(M2), C[2].max_disp, anchor_N(C[2]));
    check(peak(M1) > 1.0, "the surcharge bends the wall (teeth)");
    check(C[2].max_disp <= 1e-9, "(a) activating the anchor moves nothing");
    check(rel_gap(M2, M1) <= 1e-9, "(a) the wall's moment is carried through the activation");
    check(std::fabs(anchor_N(C[2])) <= 1e-9 * kQ1, "(a) the new anchor starts with no force");
    check(C[2].message.find("installed in this phase") != std::string::npos,
          "(a) and the phase says the anchor was installed in it");

    // (b) More load: the first stage plus the second stage.
    m::Project second = model(true);
    second.initial.struct_active = {1, 1};
    second.phases = {phase("extra", {1, 1}, {0, 1, 0})};
    const auto S = solve(second, InitialPhase::K0Procedure);
    check(S.size() == 2 && S[1].ok, "the second stage alone (wall and anchor under 30 kPa) solves");
    if (S.size() != 2 || !S[1].ok) return;
    const auto M3 = moments(C[3]), Ms = moments(S[1]);
    std::vector<double> sum(M1.size());
    for (size_t k = 0; k < M1.size() && k < Ms.size(); ++k) sum[k] = M1[k] + Ms[k];
    std::printf("      more: max|M| %.9f, N %.9f | superposed max|M| %.9f, N %.9f | gaps %.2e / %.2e\n",
                peak(M3), anchor_N(C[3]), peak(sum), anchor_N(S[1]), rel_gap(M3, sum),
                std::fabs(anchor_N(C[3]) - anchor_N(S[1])) / std::fabs(anchor_N(S[1])));
    check(std::fabs(anchor_N(S[1])) > 1.0, "the anchor carries the second stage (teeth)");
    check(rel_gap(M3, sum) <= 1e-9, "(b) the wall's moment is the sum of the two stages");
    check(std::fabs(anchor_N(C[3]) - anchor_N(S[1])) <= 1e-9 * std::fabs(anchor_N(S[1])),
          "(b) the anchor carries the second stage alone");
    check(max_abs_diff(C[3].disp, S[1].disp) <= 1e-9,
          "(b) and the displacement increment is the second stage's");

    // (c) Removing the anchor releases its force onto the wall and the soil.
    const double N3 = anchor_N(C[3]);
    // The anchor attaches to the mesh node nearest its drawn head, and pulls along the line from
    // THAT node to its fixed point; so does the equivalent point force, which lands on the same node.
    const auto& mesh = C[3].mesh;
    int head = 0;
    for (int n = 1; n < mesh.node_count; ++n)
        if (std::hypot(mesh.x[n] - kAx, mesh.y[n] - kAy) <
            std::hypot(mesh.x[head] - kAx, mesh.y[head] - kAy))
            head = n;
    const Eigen::Vector2d dir =
        Eigen::Vector2d(kFx - mesh.x[head], kFy - mesh.y[head]).normalized();
    std::printf("      the anchor attaches at (%.4f, %.4f)\n", mesh.x[head], mesh.y[head]);
    m::Project released = model(true);
    released.initial.struct_active = {1, 0};
    // The anchor pulled its head towards the fixed point with N3; removed, that pull is gone, which
    // is the soil-and-wall model loaded by N3 along the axis AWAY from the fixed point.
    released.loads[2].qx1 = released.loads[2].qx2 = -N3 * dir(0);
    released.loads[2].qy1 = released.loads[2].qy2 = -N3 * dir(1);
    released.phases = {phase("release", {1, 0}, {0, 0, 1})};
    const auto Rl = solve(released, InitialPhase::K0Procedure);
    check(Rl.size() == 2 && Rl[1].ok, "the released force alone solves");
    if (Rl.size() != 2 || !Rl[1].ok) return;
    const auto M4 = moments(C[4]), Mr = moments(Rl[1]);
    std::vector<double> change(M4.size());
    for (size_t k = 0; k < M4.size() && k < M3.size(); ++k) change[k] = M4[k] - M3[k];
    std::printf("      remove: change in max|M| %.9f | released run %.9f | gap %.2e, increment gap %.2e m\n",
                peak(change), peak(Mr), rel_gap(change, Mr), max_abs_diff(C[4].disp, Rl[1].disp));
    check(force_of(C[4], "Tie") == nullptr, "(c) the removed anchor is no longer reported");
    check(rel_gap(change, Mr) <= 1e-9, "(c) removing the anchor releases exactly its force");
    check(max_abs_diff(C[4].disp, Rl[1].disp) <= 1e-9,
          "(c) and the displacement increment is the released run's");

    // (d) A prestressed anchor activated in the chain applies its lock-off force in that phase.
    m::Project pre = model(true);
    pre.anchors[0].prestress = kN0;
    pre.phases = {phase("q1", {1, 0}, {1, 0, 0}), phase("anchor", {1, 1}, {1, 0, 0})};
    m::Project slack = model(true);
    slack.loads[2].qx1 = slack.loads[2].qx2 = kN0 * dir(0);
    slack.loads[2].qy1 = slack.loads[2].qy2 = kN0 * dir(1);
    slack.phases = {phase("q1", {1, 0}, {1, 0, 0}), phase("anchor", {1, 1}, {1, 0, 1})};
    const auto Pp = solve(pre, InitialPhase::K0Procedure);
    const auto Ps = solve(slack, InitialPhase::K0Procedure);
    check(Pp.size() == 3 && Pp[2].ok && Ps.size() == 3 && Ps[2].ok,
          "the prestressed chain and its slack-plus-force twin solve");
    if (Pp.size() != 3 || Ps.size() != 3 || !Pp[2].ok || !Ps[2].ok) return;
    const auto Mp = moments(Pp[2]), Mq = moments(Ps[2]);
    std::printf("      lock-off %.0f kN: increment %.6e vs %.6e m, max|M| %.9f vs %.9f (gap %.2e)\n", kN0,
                Pp[2].max_disp, Ps[2].max_disp, peak(Mp), peak(Mq), rel_gap(Mp, Mq));
    check(Pp[2].max_disp > 1e-5, "(d) the lock-off force moves the wall (teeth)");
    check(max_abs_diff(Pp[2].disp, Ps[2].disp) <= 1e-9 && rel_gap(Mp, Mq) <= 1e-9,
          "(d) the prestressed anchor equals a slack anchor plus its lock-off as a point force");
}

void nil_activation_on_the_interface_wall() {
    std::printf("\nthe nonlinear check: an inert anchor activated beside an interface wall\n");
    for (bool interfaces : {true, false}) {
        m::Project pr = model(false);
        pr.structs[0].iface_pos = pr.structs[0].iface_neg = interfaces;
        pr.anchors[0].EA = 1.0;
        pr.structs[1].x1 = 16; pr.structs[1].y1 = 9; pr.structs[1].x2 = 19; pr.structs[1].y2 = 9;
        pr.phases = {phase("q1", {1, 0}, {1, 0, 0}), phase("anchor", {1, 1}, {1, 0, 0})};
        const auto R = solve(pr, InitialPhase::K0Procedure);
        const std::string tag = interfaces ? "wall with interfaces" : "wall without interfaces";
        check(R.size() == 3 && R[1].ok && R[2].ok, tag + ": the chain solves");
        if (R.size() != 3 || !R[2].ok) continue;
        const auto Ma = moments(R[1]), Mb = moments(R[2]);
        double joint_gap = 0.0, joint_teeth = 0.0;
        for (size_t i = 0; i < R[2].interface_forces.size() && i < R[1].interface_forces.size(); ++i)
            for (size_t k = 0; k < R[2].interface_forces[i].stations.size() &&
                               k < R[1].interface_forces[i].stations.size(); ++k) {
                const auto& a = R[1].interface_forces[i].stations[k];
                const auto& b = R[2].interface_forces[i].stations[k];
                joint_gap = std::max({joint_gap, std::fabs(a.tau - b.tau), std::fabs(a.sigma_n - b.sigma_n)});
                joint_teeth = std::max({joint_teeth, std::fabs(a.tau), std::fabs(a.sigma_n)});
            }
        std::printf("      %s: max|M| %.6f -> %.6f (gap %.2e), increment %.3e / %.3e m, joint gap %.2e kPa\n",
                    tag.c_str(), peak(Ma), peak(Mb), rel_gap(Ma, Mb), R[1].max_disp, R[2].max_disp,
                    joint_gap);
        check(rel_gap(Ma, Mb) <= 1e-6, tag + ": the wall's moment is carried through the activation");
        check(R[2].max_disp <= 1e-6 * R[1].max_disp, tag + ": nothing moves");
        if (interfaces)
            check(joint_teeth > 10.0 && joint_gap <= 1e-6 * joint_teeth,
                  tag + ": the joint's stresses are carried through the activation");
    }
}

}  // namespace

int main() {
    std::printf("KV-STR-012: activating or removing a structure carries every other one\n");
    superposition_on_a_linear_model();
    nil_activation_on_the_interface_wall();
    if (g_failures == 0) {
        std::printf("\nOK: structure-set changes carry the structural state\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
