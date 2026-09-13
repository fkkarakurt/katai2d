// An interface inside a CONSOLIDATION phase -- the case a real retaining wall is, and the one the
// previous increment deliberately refused. The refusal said the split seam "would silently make the
// joint impermeable, which is a modelling claim". That was right about the danger and wrong about
// the default: the model already carries the interface's cross permeability as an input
// (`flow_barrier`), and its default is FULLY PERMEABLE -- the flow net runs through the line. So
// the fix was not to keep refusing but to read the input that was already there.
//
// WHAT AN INTERFACE ADDS TO A COUPLED PHASE. It splits the mesh, so the joint carries two pore
// pressures where the ground carried one. Which of them the water sees is an input, not a
// consequence of the splitting:
//   fully permeable  -- one pressure. The two nodes SHARE a pore equation, so continuity and the
//                       flux balance across the joint hold by construction rather than by a
//                       constraint that could be assembled slightly wrong;
//   impermeable      -- two pressures, which is what the bare split already gives;
//   semi-permeable   -- a conductance q = dh/R, which is neither, and is refused.
//
// THE TEST IS THE DIFFERENCE THE FLAG MAKES. One model, one surcharge on one side of a wall, and
// the observable is the pore pressure across the seam: tied to machine zero when the joint is
// permeable, and a real head difference when it is not. An implementation that ignored the flag
// would produce one of those in both cases, and nothing else in the result would look wrong.
//
// The mechanical half comes with the same limit the rest of this phase's structural branch has --
// the joint is ELASTIC and cannot slip -- and the same treatment: it is measured and reported
// (K2D-A016) rather than declared, because a joint that cannot slip is stiffer than the real one,
// so the wall deflects less and attracts more load, which is the unsafe direction.
//
// verify: KV-STR-007
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, structs[i].flow_barrier / hyd_res) for the three cross-permeability cases and what each means (fully permeable, the default: no effect on flow; impermeable: separate pore-pressure degrees of freedom on the two sides; semi-permeable: a hydraulic resistance d/k); and the drained limit of Biot consolidation, as in KV-STR-006 -- once the excess pore pressure is gone the coupled problem IS the drained problem, which this program solves by an independent path (solve_nonlinear, with a real Coulomb return on the joint)
//   locator:  a permeable joint is ONE pressure: p_left - p_right = 0 identically, because the two nodes share the pore equation. An impermeable joint is two, and a surcharge on one side cannot cross it. Coulomb capacity at a station: tau_max = c_i - sigma_n tan(phi_i), with c_i = R_inter c and tan(phi_i) = R_inter tan(phi)
//   quantity: the maximum |p_left - p_right| across the seam EARLY in the dissipation (Tv = 0.05) for a permeable and an impermeable joint [kPa]; the settlement and the joint's peak shear at the end of a consolidation phase run to Tv = 4, against the same model as a drained Plastic phase [m, kPa]; the Coulomb demand/capacity the elastic joint reaches; and whether a joint past its capacity is reported
//   expected: permeable -> 0 (one equation, so the difference is a number minus itself); impermeable -> of the order of the surcharge (measured 12.62 kPa of 50 kPa; 12.83 when the case was first recorded, and the assertion is the order, above 10% of the surcharge). With a joint below its capacity the coupled phase must reach the drained answer in BOTH the settlement and the joint's shear; with one above it (1.19x here) it must not, and the run must say so; semi-permeable must be refused with what it actually is
//   band:     1e-12 kPa on the permeable seam, measured 0.000e+00 exactly. 5e-4 relative on the drained-limit pair, measured -0.0042% on the settlement and 0.0000% on the joint's shear (53.2745 kPa against 53.2745). Both are asserted on a LINEAR-ELASTIC soil on purpose: with Mohr-Coulomb the same pair comes out at +0.15% and the cause is not the joint but the soil's stress path -- the coupled phase loads it undrained and the Plastic phase drained, and a yield surface remembers the difference. That figure is measured and printed beside the assertion rather than used to widen it, because a band wide enough to hold it would also be wide enough to hide a joint that was not in the system at all
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kW = 20.0, kH = 10.0;
constexpr double kXw = 10.0;        // the wall line
constexpr double kYtoe = 3.0;       // it stops short of the base, so water can go under it
constexpr double kE = 5000.0, kNu = 0.2;
constexpr double kPerm = 1.0e-3;
constexpr double kQ = 50.0;         // surcharge, on the LEFT half only

double eoed() { return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }
double cv() { return kPerm * eoed() / katai::app::kGammaWater; }

// The block with a vertical interface at x = kXw. `barrier`: the model's cross-permeability
// vocabulary (0 permeable, 1 impermeable, 2 semi-permeable). `rinter`: 1 = as strong as the soil
// (it will not slip), small = a weak joint that will.
m::Project block(int barrier, double rinter, bool plastic_phase, double tv = 4.0,
                 bool elastic_soil = true) {
    m::Project pr;
    m::Material s; s.model = elastic_soil ? m::SoilModel::LinearElastic : m::SoilModel::MohrCoulomb;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 17.0; s.gamma_sat = 20.0; s.e_init = 0.7;
    s.c = 25.0; s.phi = 25.0; s.psi = 0.0;
    s.kx = kPerm; s.ky = kPerm;
    s.rinter_rigid = rinter >= 1.0;
    s.Rinter = rinter;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};
    P.edge_head = {0.0, 0.0, kH, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::StructElement e; e.kind = m::StructKind::Interface; e.name = "Joint";
    e.x1 = kXw; e.y1 = kYtoe; e.x2 = kXw; e.y2 = kH; e.material = -1;
    e.flow_barrier = barrier;
    pr.structs.push_back(e);

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 2.0; L.y1 = kH; L.x2 = kXw; L.y2 = kH;       // the LEFT side only
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};

    m::Phase ph; ph.load_active = {1};
    if (plastic_phase) {
        ph.name = "Drained (long term)"; ph.type = m::PhaseType::Plastic;
    } else {
        ph.name = "Consolidation"; ph.type = m::PhaseType::Consolidation;
        ph.duration = tv * kH * kH / cv();
        ph.time_steps = 60;
    }
    pr.phases.push_back(ph);
    return pr;
}

struct Answer {
    bool ok = false;
    double settle = 0.0, tau = 0.0, util = 0.0, seam_gap = 0.0;
    bool warned_slip = false, any_slip = false;
    bool warned_barrier = false;   // K2D-A010: a declared barrier the phase cannot read
    std::vector<double> excess_pore;
    std::string msg;
};

// The largest pore-pressure difference between two nodes that sit at the same place on the wall
// line -- i.e. across the seam. The split leaves coincident duplicates, so "the same place" is a
// coordinate match, which is what a reader of the result would do too.
double seam_pore_gap(const katai::core::SolveResult& R) {
    double worst = 0.0;
    for (int a = 0; a < R.mesh.node_count; ++a) {
        if (std::fabs(R.mesh.x[a] - kXw) > 1e-9 || R.mesh.y[a] < kYtoe - 1e-9) continue;
        for (int b = a + 1; b < R.mesh.node_count; ++b) {
            if (std::fabs(R.mesh.x[b] - kXw) > 1e-9) continue;
            if (std::fabs(R.mesh.y[b] - R.mesh.y[a]) > 1e-9) continue;
            if (a < (int)R.excess_pore.size() && b < (int)R.excess_pore.size())
                worst = std::fmax(worst, std::fabs(R.excess_pore[a] - R.excess_pore[b]));
        }
    }
    return worst;
}

Answer run(const m::Project& pr) {
    Answer a;
    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    if (!M.ok) { a.msg = M.message; return a; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2) { a.msg = "phases did not run"; return a; }
    const auto& R = res[1];
    a.ok = R.ok; a.msg = R.message;
    if (!R.ok) return a;
    for (int n = 0; n < R.mesh.node_count; ++n)
        a.settle = std::fmax(a.settle, std::fabs(R.disp(2 * n + 1)));
    if (!R.interface_forces.empty()) {
        a.tau = R.interface_forces.front().max_abs_tau;
        a.util = R.interface_forces.front().max_utilisation;
        a.any_slip = R.interface_forces.front().any_slip;
    }
    a.seam_gap = seam_pore_gap(R);
    a.excess_pore = R.excess_pore;
    for (const auto& d : R.diagnostics) {
        if (d.code == std::string("K2D-A016")) a.warned_slip = true;
        if (d.code == std::string("K2D-A010")) a.warned_barrier = true;
    }
    return a;
}

// ---------------------------------------------------------------------------------------------
void test_cross_permeability_is_read() {
    std::printf("-- (1) the joint's cross permeability decides what the water does at it --\n");
    // EARLY in the dissipation (Tv = 0.05), because that is when there is a difference to see:
    // by the end of a phase run to Tv = 4 the excess pore pressure has gone from both sides and
    // the two joints agree about nothing in particular. Measuring at the end would have made this
    // test pass for a build that ignored the flag entirely.
    const Answer perm = run(block(0, 1.0, false, 0.05));
    const Answer imp  = run(block(1, 1.0, false, 0.05));
    check(perm.ok, "a consolidation phase with an interface runs (it used to be refused)");
    if (!perm.ok) { std::printf("   (%s)\n", perm.msg.c_str()); return; }
    check(imp.ok, "and so does the same model with the joint impermeable");
    if (!imp.ok) { std::printf("   (%s)\n", imp.msg.c_str()); return; }

    std::printf("   max |p_left - p_right| across the seam:  permeable %.3e kPa,  impermeable %.4f kPa\n",
                perm.seam_gap, imp.seam_gap);
    check(perm.seam_gap < 1e-12,
          "permeable: one pressure at the joint -- the two sides share the equation");
    check(imp.seam_gap > 0.10 * kQ,
          "impermeable: two pressures, and the surcharge on one side cannot cross");
    check(imp.seam_gap > 1e6 * std::fmax(perm.seam_gap, 1e-300),
          "so the flag is READ, rather than the split deciding on its own");
    // K2D-A010 exists for a barrier the phase CANNOT read. It used to fire here too -- for every
    // plate or interface in a coupled phase -- on the very seam this test shows being read.
    check(!perm.warned_barrier && !imp.warned_barrier,
          "and K2D-A010 stays silent on both: a split seam's cross permeability is read");

    // The same line as a WALL WITH INTERFACES: the plate's own flow_barrier rides on the wall's
    // seam, which is a different list in the driver from a bare interface's.
    const auto wall = [](m::Project pr) {
        m::PlateMaterial pm; pm.EA = 1.0e6; pm.EI = 1.0e4;
        pr.plates.push_back(pm);
        pr.structs[0].kind = m::StructKind::Plate;
        pr.structs[0].name = "Wall";
        pr.structs[0].material = 0;
        pr.structs[0].iface_pos = pr.structs[0].iface_neg = true;
        return pr;
    };
    const Answer wimp = run(wall(block(1, 1.0, false, 0.05)));
    check(wimp.ok, "a wall with interfaces, impermeable, in a consolidation phase runs");
    if (wimp.ok) {
        std::printf("   wall with interfaces, impermeable: %.4f kPa\n", wimp.seam_gap);
        check(wimp.seam_gap > 1e6 * std::fmax(perm.seam_gap, 1e-300) && !wimp.warned_barrier,
              "a wall's seam is read the same way, and raises no K2D-A010");
    }

    // The FULLY-COUPLED phase is handed the same seam ties, and nothing else checked that it reads
    // them. Asked at Tv = 0.5, where a head difference is still there to see.
    const auto coupled = [](m::Project pr) {
        pr.phases[0].name = "Fully coupled";
        pr.phases[0].type = m::PhaseType::FullyCoupled;
        return pr;
    };
    const Answer fperm = run(coupled(block(0, 1.0, false, 0.5)));
    const Answer fimp  = run(coupled(block(1, 1.0, false, 0.5)));
    const Answer fwall = run(coupled(wall(block(1, 1.0, false, 0.5))));
    check(fperm.ok && fimp.ok && fwall.ok,
          "in a fully-coupled phase the permeable joint, the impermeable joint and the wall run");
    if (!fperm.ok || !fimp.ok || !fwall.ok) {
        std::printf("   (%s | %s | %s)\n", fperm.msg.c_str(), fimp.msg.c_str(), fwall.msg.c_str());
        return;
    }
    std::printf("   fully coupled, Tv = 0.5:  permeable %.3e kPa,  impermeable %.4f kPa,  "
                "wall with interfaces %.4f kPa\n", fperm.seam_gap, fimp.seam_gap, fwall.seam_gap);
    check(fperm.seam_gap < 1e-12, "fully coupled, permeable: one pressure at the joint");
    check(fimp.seam_gap > 1e6 * std::fmax(fperm.seam_gap, 1e-300) &&
              fwall.seam_gap > 1e6 * std::fmax(fperm.seam_gap, 1e-300),
          "fully coupled, impermeable joint and wall: two -- the flag is read in this phase too");
    check(!fperm.warned_barrier && !fimp.warned_barrier && !fwall.warned_barrier,
          "and K2D-A010 stays silent there too");
}

void test_unsplit_barrier_is_said() {
    std::printf("\n-- (1b) a barrier on a line the mesh was NOT split along is not read, and says so --\n");
    // The same block, the same line, drawn as a PLATE WITHOUT INTERFACES. Nothing splits the mesh,
    // so there is no second pore equation for "impermeable" to act on.
    const auto plate_line = [](int barrier) {
        m::Project pr = block(barrier, 1.0, false, 0.05);
        m::PlateMaterial pm; pm.EA = 1.0e6; pm.EI = 1.0e4;
        pr.plates.push_back(pm);
        pr.structs[0].kind = m::StructKind::Plate;
        pr.structs[0].name = "Wall";
        pr.structs[0].material = 0;
        return pr;
    };
    const Answer perm = run(plate_line(0));
    const Answer imp  = run(plate_line(1));
    check(perm.ok && imp.ok, "both plate models ran");
    if (!perm.ok || !imp.ok) return;
    check(!imp.excess_pore.empty() && imp.excess_pore == perm.excess_pore,
          "an impermeable plate without interfaces gives the permeable plate's pore field bit for bit");
    check(imp.warned_barrier, "and the run says the declaration was not read (K2D-A010)");
    check(!perm.warned_barrier,
          "while the fully permeable plate, which behaves exactly as declared, raises nothing");

    // The same pair in a fully-coupled phase, at Tv = 0.5 like the joints above.
    const auto coupled_plate = [&plate_line](int barrier) {
        m::Project pr = plate_line(barrier);
        pr.phases[0].type = m::PhaseType::FullyCoupled;
        pr.phases[0].duration = 0.5 * kH * kH / cv();
        return pr;
    };
    const Answer fperm = run(coupled_plate(0));
    const Answer fimp  = run(coupled_plate(1));
    check(fperm.ok && fimp.ok, "both plate models ran in a fully-coupled phase");
    if (!fperm.ok || !fimp.ok) return;
    check(!fimp.excess_pore.empty() && fimp.excess_pore == fperm.excess_pore && fimp.warned_barrier &&
              !fperm.warned_barrier,
          "fully coupled: the same -- not read, bit for bit, and said only where it was declared");
}

void test_drained_limit_with_a_joint() {
    std::printf("\n-- (2) with a joint that does not slip, the drained limit still holds --\n");
    const Answer drained = run(block(0, 1.0, true));
    const Answer consol  = run(block(0, 1.0, false));
    check(drained.ok && consol.ok, "both paths solved");
    if (!drained.ok || !consol.ok) return;
    std::printf("   settlement %.6e vs %.6e m  (%+.4f%%);  interface |tau|max %.4f vs %.4f kPa\n",
                consol.settle, drained.settle, 100.0 * (consol.settle / drained.settle - 1.0),
                consol.tau, drained.tau);
    std::printf("   peak Coulomb demand/capacity: coupled %.3f, drained reference %.3f (it slipped: %d)\n",
                consol.util, drained.util, (int)drained.any_slip);
    check(consol.util < 1.0, "the joint does not reach its capacity here (so the two are comparable)");
    check(!consol.warned_slip, "and nothing is reported about slip, because there is none");
    check(std::fabs(consol.settle / drained.settle - 1.0) < 5e-4,
          "settlement at the end of consolidation = the drained answer");
    check(std::fabs(consol.tau / drained.tau - 1.0) < 5e-4,
          "and so does the shear the joint carries -- the interface is IN the coupled system");

    // WHAT THE FIRST VERSION OF THIS COMPARISON MEASURED, AND WHY IT IS NOT ASSERTED. With a
    // Mohr-Coulomb soil the same two runs differ by ~0.15% rather than ~0.004%, and it is not the
    // joint: the drained reference does not slip either way. It is the SOIL's stress path. The
    // coupled phase applies the load undrained -- the water takes it first, the effective stress
    // arrives gradually -- while the Plastic phase applies it drained from the start, and a
    // material with a yield surface remembers the difference. The sibling test of the coupled
    // plastic core says the same thing about the same limit. So the assertion above is made on a
    // linear-elastic soil, where the paths coincide, and the plastic figure is measured and
    // printed rather than folded into a band wide enough to hide the joint.
    const Answer d_mc = run(block(0, 1.0, true, 4.0, /*elastic_soil=*/false));
    const Answer c_mc = run(block(0, 1.0, false, 4.0, /*elastic_soil=*/false));
    if (d_mc.ok && c_mc.ok)
        std::printf("   the same pair on a Mohr-Coulomb soil: %+.4f%% -- the SOIL's path "
                    "dependence, measured, not asserted\n",
                    100.0 * (c_mc.settle / d_mc.settle - 1.0));
    check(d_mc.ok && c_mc.ok, "the Mohr-Coulomb pair also solves");
}

void test_a_joint_that_would_slip_says_so() {
    std::printf("\n-- (3) a joint that WOULD slip is reported, not quietly held --\n");
    const Answer weak = run(block(0, 0.02, false));
    check(weak.ok, "the phase with a weak joint still solves");
    if (!weak.ok) { std::printf("   (%s)\n", weak.msg.c_str()); return; }
    std::printf("   peak Coulomb demand/capacity = %.2fx\n", weak.util);
    check(weak.util > 1.0, "the elastic joint carries more shear than it could");
    check(weak.warned_slip, "K2D-A016: and the run says so, naming the joint and the exceedance");
    // The reference the reader is sent to: the same stage as a Plastic phase, where the Coulomb
    // return actually runs. It must give a DIFFERENT answer -- otherwise the warning would be
    // about nothing.
    const Answer drained_weak = run(block(0, 0.02, true));
    if (drained_weak.ok)
        std::printf("   the Plastic phase it points at: |tau|max %.4f vs %.4f kPa in the elastic one\n",
                    drained_weak.tau, weak.tau);
    check(drained_weak.ok && drained_weak.tau < weak.tau,
          "and the Coulomb path really does carry less shear there");
}

void test_semi_permeable_is_refused() {
    std::printf("\n-- (4) what is between the two ends of the scale is refused, with the reason --\n");
    const Answer semi = run(block(2, 1.0, false));
    check(!semi.ok, "a semi-permeable joint in a consolidation phase is refused");
    if (!semi.ok) {
        std::printf("   (%s)\n", semi.msg.c_str());
        check(semi.msg.find("conductance") != std::string::npos,
              "and the message says what it is, not merely that it is unsupported");
    }
}

}  // namespace

int main() {
    std::printf("Interfaces inside a consolidation phase\n\n");
    test_cross_permeability_is_read();
    test_unsplit_barrier_is_said();
    test_drained_limit_with_a_joint();
    test_a_joint_that_would_slip_says_so();
    test_semi_permeable_is_refused();
    if (g_failures == 0) {
        std::printf("\nOK: the joint carries shear, the flag carries the water, and neither is silent\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
