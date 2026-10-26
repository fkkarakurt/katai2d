// An embedded pile row inside a CONSOLIDATION phase, its shaft slipping.
//
// A consolidation phase refused an embedded beam: its structures were solved as one elastic
// matrix, and a pile's skin and base are defined by the resistance at which they slip. The
// structures of the phase are now evaluated by their own laws at every iterate (StructureLaws),
// so the pile can be loaded while the ground around and beneath it consolidates.
//
// THE ORACLE IS THE SHAFT AT ITS CAPACITY. Where the skin slips it carries exactly its capacity,
// t = T_skin / Ls per metre of pile and per metre of row, so the axial force changes along that
// part of the pile by t per metre -- whatever the soil, whatever the time. The pile here is pushed
// down at a head tied to the ground surface; its lower half slips throughout, and there the force
// must fall linearly from the toe by t (y - y_toe). The head itself moves with the ground it is
// tied to, so the skin at the top holds, and the load the pile takes at the head is the system's
// to decide; the closed form is asserted where it applies.
//
// THE DRAINED LIMIT, as in KV-STR-006: at the end of its dissipation a consolidation phase must
// land where the drained Plastic phase does. It does so exactly where the end state does not depend
// on the path -- the CONTROL, the same pile with an unlimited skin, where nothing slips (to
// 0.0014%). With the shaft slipping the two are 0.24% apart, and that is the skin that holds at
// the head: a point that sticks and then slips remembers the order its loads came in, and the
// consolidation phase loads the ground undrained first and drained after, the Plastic phase drained
// throughout. The control's toe force is 26% from the slipping pile's, so the slip is what the
// first comparison measures.
//
// verify: KV-STR-023
//   oracle:   closed_form
//   source:   the embedded beam (pile row) as KATAI 2D defines it, KATAI 2D input contract (docs/k2d-format.md, embedded: Tskin_max / Fmax_base / Lspacing): an axial skin resistance per metre of pile, divided by the out-of-plane spacing for the row, at which the shaft slips; and the drained limit of Biot consolidation, as in KV-STR-006 -- once the excess pore pressure is gone the coupled problem IS the drained one, which this program solves by an independent path (solve_nonlinear)
//   locator:  a linear-elastic clay block 20 x 10 m, E = 5000 kPa, nu = 0.2, k = 1e-3 m/day, drained at the top; a pile at x = 10 m from y = 2 to 10 m (D = 0.4 m, Ls = 2.5 m, hinged at the head), loaded at the head with 100 kN/m; T_skin = 10 kN/m per pile, t = 4 kN/m per metre of row. Consolidation to Tv = 4 (cv = k Eoed / gamma_w)
//   quantity: the axial force along the lower half of the pile (y < 6 m) at the end of the consolidation phase [kN/m]; the head settlement [m] and the toe force [kN/m] against the drained Plastic phase, with the skin slipping and with an unlimited one
//   expected: N(y) = N(toe) - t (y - y_toe), t = 4, where every skin point has slipped; the drained phase's settlement and toe force
//   band:     1e-6 of the head load on the slipping shaft (measured 9e-13: the capacity is enforced exactly); 5e-4 on the control's settlement and toe force against the drained phase (measured -0.0014%, -0.0010%); 0.5% on the slipping pile's (measured +0.24%, +0.28%: the path the stick-slip point at the head takes); the control's toe force more than 10% from the slipping pile's (26%)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kW = 20.0, kH = 10.0;
constexpr double kE = 5000.0, kNu = 0.2, kPerm = 1.0e-3;
constexpr double kXp = 10.0, kYtoe = 2.0;
constexpr double kP = 100.0;          // head load [kN/m of row]
constexpr double kLs = 2.5;           // pile spacing [m]
constexpr double kTskin = 10.0;       // skin resistance per pile [kN/m]
constexpr double kLower = 6.0;        // the lower half of the shaft: from the toe at 2 m to 6 m
constexpr int kLowerPoints = 9;       // ... three beam elements of the six, three skin points each

double eoed() { return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }
double cv() { return kPerm * eoed() / katai::app::kGammaWater; }

m::Project block(bool plastic_phase, double tskin) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 17.0; s.gamma_sat = 20.0; s.e_init = 0.7;
    s.kx = kPerm; s.ky = kPerm;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0}; P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};
    P.edge_head = {0.0, 0.0, kH, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::EmbeddedBeamMaterial em; em.name = "Pile";
    em.E = 3.0e7; em.diameter = 0.4; em.Lspacing = kLs; em.gamma = 0.0;
    em.Tskin_max = tskin; em.Fmax_base = 0.0;
    pr.embedded.push_back(em);
    m::StructElement e; e.kind = m::StructKind::EmbeddedBeam; e.name = "Pile";
    e.x1 = kXp; e.y1 = kYtoe; e.x2 = kXp; e.y2 = kH; e.material = 0;
    pr.structs.push_back(e);

    m::Load L; L.kind = m::LoadKind::Point; L.name = "Head load";
    L.x1 = kXp; L.y1 = kH; L.qx1 = 0.0; L.qy1 = -kP;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};

    m::Phase ph; ph.load_active = {1};
    if (plastic_phase) {
        ph.name = "Drained (long term)"; ph.type = m::PhaseType::Plastic;
    } else {
        ph.name = "Consolidation"; ph.type = m::PhaseType::Consolidation;
        ph.duration = 4.0 * kH * kH / cv();
        ph.time_steps = 60;
        ph.consol_load_at_start = true;
    }
    pr.phases.push_back(ph);
    return pr;
}

struct Answer {
    bool ok = false;
    std::string msg;
    double head_settle = 0.0, n_toe = 0.0, n_head = 0.0, pore = 0.0;
    double worst_slope = 0.0;      // the slipping lower shaft against N(y) = N(toe) - t (y - toe)
    int skin_points = 0, slipped = 0, lower_points = 0, lower_slipped = 0;
};

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
        if (std::fabs(R.mesh.x[n] - kXp) < 1e-9 && std::fabs(R.mesh.y[n] - kH) < 1e-9)
            a.head_settle = std::fabs(R.disp(2 * n + 1));
    const double t = kTskin / kLs;
    for (const auto& f : R.struct_forces) {
        if (f.name != "Pile") continue;
        double ytop = -1e300, ybot = 1e300;
        for (const auto& st : f.stations) {
            if (st.y > ytop) { ytop = st.y; a.n_head = st.N; }
            if (st.y < ybot) { ybot = st.y; a.n_toe = st.N; }
        }
        for (const auto& st : f.stations)   // the lower half of the shaft
            if (st.y < kLower - 1e-9)
                a.worst_slope = std::max(a.worst_slope,
                                         std::fabs(st.N - (a.n_toe - t * (st.y - kYtoe))) / kP);
    }
    // The skin points run from the toe up, three to a beam element; the lower half of the shaft is
    // its first kLowerPoints.
    for (double s : R.struct_state.embedded_skin_slip) {
        const bool lower = a.skin_points < kLowerPoints;
        ++a.skin_points;
        if (s < -1e-12) { ++a.slipped; if (lower) ++a.lower_slipped; }   // the pile moved down past the soil
        if (lower) ++a.lower_points;
    }
    if (!R.consol_excess_pore.empty()) a.pore = R.consol_excess_pore.back();
    return a;
}
}  // namespace

int main() {
    std::printf("== an embedded pile row slipping in a consolidation phase (KV-STR-023) ==\n");
    std::printf("   cv = %.4f m2/day, Tv = 4 at t = %.0f day; skin capacity t = %.1f kN/m per metre of row\n",
                cv(), 4.0 * kH * kH / cv(), kTskin / kLs);
    const Answer c = run(block(false, kTskin));
    check(c.ok, "a consolidation phase with an embedded beam solves (it used to be refused)");
    if (!c.ok) { std::printf("   (%s)\n", c.msg.c_str()); return 1; }
    std::printf("   consolidation: N head %.4f, toe %.4f kN/m; skin points slipped %d of %d (%d of %d in the lower half);"
                " the lower half against N(toe) - t (y - toe): %.2e of P; excess pore left %.2e kPa\n",
                c.n_head, c.n_toe, c.slipped, c.skin_points, c.lower_slipped, c.lower_points, c.worst_slope, c.pore);
    check(c.lower_points > 0 && c.lower_slipped == c.lower_points, "the lower half of the shaft has slipped throughout");
    check(c.worst_slope < 1e-6, "and there the axial force falls by exactly the skin capacity per metre (1e-6 of P)");

    const Answer d = run(block(true, kTskin));
    check(d.ok, "the drained Plastic phase solves");
    if (d.ok) {
        std::printf("   slipping: head settlement consolidation %.8e, drained %.8e m (%+.4f%%); N toe %.5f vs %.5f (%+.4f%%)\n",
                    c.head_settle, d.head_settle, 100.0 * (c.head_settle / d.head_settle - 1.0), c.n_toe, d.n_toe,
                    100.0 * (c.n_toe / d.n_toe - 1.0));
        check(d.worst_slope < 1e-6 && d.lower_slipped == d.lower_points,
              "the drained phase slips the same lower half at the same capacity (1e-6)");
        check(std::fabs(c.head_settle / d.head_settle - 1.0) < 5e-3 && std::fabs(c.n_toe / d.n_toe - 1.0) < 5e-3,
              "and ends within 0.5% of the consolidation phase: the history of the skin that holds at the head");
    }

    // The control: an unlimited skin. Nothing slips, so the end state cannot depend on the path,
    // and the consolidation phase must land on the drained one as the elastic structures of
    // KV-STR-006 do.
    const Answer uc = run(block(false, 0.0)), ud = run(block(true, 0.0));
    check(uc.ok && ud.ok, "the control, an unlimited skin, solves on both paths");
    if (uc.ok && ud.ok) {
        std::printf("   control: head settlement %+.5f%%, N toe %.5f vs %.5f (%+.5f%%); against the slipping toe %.1f%%\n",
                    100.0 * (uc.head_settle / ud.head_settle - 1.0), uc.n_toe, ud.n_toe,
                    100.0 * (uc.n_toe / ud.n_toe - 1.0), 100.0 * std::fabs(uc.n_toe / c.n_toe - 1.0));
        check(uc.slipped == 0, "nothing slips in the control");
        check(std::fabs(uc.head_settle / ud.head_settle - 1.0) < 5e-4 && std::fabs(uc.n_toe / ud.n_toe - 1.0) < 5e-4,
              "the control's consolidation lands on its drained answer (5e-4)");
        check(std::fabs(uc.n_toe / c.n_toe - 1.0) > 0.10, "and its toe force is far from the slipping shaft's -- the slip is what was measured above");
    }
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
