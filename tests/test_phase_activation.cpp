// A wall with interfaces, or an interface, activated and deactivated per phase.
//
// Both split the mesh along their line, and until 2026-09 the split -- and with it the element --
// had to be present in every phase: a phase that left one out was refused ("per phase is not
// supported yet"). So the everyday sequence of an excavation, a K0 initial phase and the wall
// installed in the first construction stage, could not be modelled; the wall had to exist in the
// geostatic phase itself.
//
// The mesh is now split in every phase whether the element is active or not, so that every phase has
// the same nodes. In a phase where the element is inactive nothing is built on the seam and its two
// sides are tied into one unknown each: the ground is continuous. An element activated in a later
// phase is installed on the ground as the phase before left it (the installation cohorts of
// KV-STR-012), and a joint takes over the stress the ground carried across its line: the recovered
// committed stress, normal n.sigma.n as its sigma_n0 and shear t.sigma.n as an initial slip offset.
//
//   (a) an element inactive in every phase is the same model as no split at all: the displacement
//       field equals the one of the same line drawn as an inactive geogrid, which constrains the mesh
//       identically and splits nothing;
//   (b) the usual sequence -- K0, then the wall activated -- equals the wall present from the K0 phase:
//       after a level K0 phase the recovered stress IS K0 sigma'_v, the joint seed the wall has always
//       been given;
//   (c) activated after a stage that is not geostatic, the joint's seed is the recovered stress, which
//       is not the discrete traction exactly: the activation phase is not a perfect nil, and the case
//       states how far from one it is -- and that the shear seed is what brings it there;
//   (d) activation in a consolidation phase is refused (the seeding runs in Plastic phases), and so is
//       deactivating a joint whose two sides have moved apart while a structure standing on its line
//       is carried on.
//
// verify: KV-STR-013
//   oracle:   independent_path
//   source:   the KATAI 2D input contract for staged activation (docs/k2d-format.md, phases[i].struct) -- an inactive structure is absent from the phase's model; for (b) the geostatic stress field of the K0 procedure (horizontal effective stress K0 sigma'_v, vertical sigma'_v) and the wished-in-place normal stress it gives a wall's interfaces, sigma_n0 = (K0 n_x^2 + n_y^2) sigma'_v; for (c) the Cauchy traction on a line of unit normal n and tangent t, sigma_n = n.sigma.n and tau = t.sigma.n
//   locator:  a Mohr-Coulomb or linear-elastic block 20 m x 10 m (E = 30000 kPa, nu = 0.3, c = 5 kPa, phi = 30 deg, gamma = 18 kN/m3, K0 = 0.5), base fixed, sides on rollers; a wall x = 10 m, y = 4..10 m (EA = 1.2e7 kN/m, EI = 1.6e5 kNm2/m) with interfaces both sides (R_inter = 1); a 50 kPa surcharge on 0 <= x <= 8 m. (a) K0 -> surcharge with the wall (or a standalone interface on the same line) inactive in both phases, against the same line as an inactive geogrid. (b) K0 (wall inactive) -> wall activated -> surcharge, against K0 (wall active) -> nil -> surcharge, with the wall weightless and at w = 9.6 kN/m/m. (c) gravity loading -> surcharge (wall inactive) -> wall activated, nothing else
//   quantity: (a) the nodal displacement field at the original nodes [m]; (b) the surcharge phase's displacement field [m] and the wall's moment at every station [kNm/m]; (c) the largest displacement increment of the activation phase against the surcharge phase's, and the wall's largest moment in it [m; kNm/m]
//   expected: (a) equal; (b) equal; (c) a displacement increment below 1e-3 of the surcharge phase's and a wall moment below 0.05 kNm/m
//   band:     (a) 1e-12 m, measured 0 -- the tie gives the twins their original's equation, and the equation numbering of the unsplit mesh, so the two runs assemble the same system; (b) 1e-12 m and 1e-9 relative, measured 0 and 0 for the weightless wall and 2.3e-15 m and 7.4e-13 for w = 9.6 (whose own weight settles it in the activation phase exactly as it settles the wall active from the K0 phase, 1.349e-3 m in both); (c) measured 2.15e-4 of the surcharge increment (2.7509e-6 m against 1.2821e-2) and 0.00996 kNm/m, where a seed without its shear part gives 3.6e-3 and 0.47 kNm/m and one with the shear part reversed 7.2e-3 and 0.95 -- the band separates them. On a linear-elastic block the same activation moves 2.6e-4 m whatever the shear seed, because that ground carried tension across the seam and a joint with no tensile strength opens: that is the joint, not the seed, and it is printed rather than asserted

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

// Struct order: 0 the line (wall / interface / geogrid), 1 an anchor (inactive unless asked).
m::Project block(m::SoilModel model, m::StructKind kind, bool interfaces, double w = 0.0) {
    m::Project pr;
    pr.name = "phase activation";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.mesh.elem_size = 1.0; pr.mesh.order = 6; pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "soil"; s.model = model;
    s.E = 30000; s.nu = 0.3; s.e_init = 0.6; s.c = 5.0; s.phi = 30.0; s.psi = 0.0;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.rinter_rigid = false; s.Rinter = 1.0;
    s.k0_auto = false; s.k0 = 0.5;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    m::PlateMaterial pm; pm.name = "Wall";
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = w;
    pr.plates.push_back(pm);
    m::GeogridMaterial gm; gm.name = "Grid"; gm.EA = 1.0;
    pr.geogrids.push_back(gm);
    m::AnchorMaterial am; am.name = "Tie"; am.EA = 2.0e5;
    pr.anchors.push_back(am);
    m::StructElement line; line.kind = kind; line.name = "Line"; line.material = kind == m::StructKind::Interface ? -1 : 0;
    line.x1 = 10; line.y1 = 4; line.x2 = 10; line.y2 = 10;
    line.iface_pos = line.iface_neg = interfaces;
    pr.structs.push_back(line);
    m::StructElement tie; tie.kind = m::StructKind::Anchor; tie.name = "Tie"; tie.material = 0;
    tie.x1 = 10; tie.y1 = 10; tie.x2 = -5; tie.y2 = 12;
    pr.structs.push_back(tie);
    m::Load q; q.kind = m::LoadKind::Distributed; q.name = "Surcharge";
    q.x1 = 0; q.y1 = 10; q.x2 = 8; q.y2 = 10; q.qx1 = q.qx2 = 0; q.qy1 = q.qy2 = -50.0;
    pr.loads.push_back(q);
    return pr;
}

m::Phase phase(const char* name, m::PhaseType type, std::vector<char> st, std::vector<char> ld) {
    m::Phase ph; ph.name = name; ph.type = type;
    ph.struct_active = std::move(st); ph.load_active = std::move(ld);
    return ph;
}

std::vector<SolveResult> solve(m::Project pr, InitialPhase ip, std::vector<char> initial_st,
                               std::vector<m::Phase> phases) {
    pr.initial.struct_active = std::move(initial_st);
    pr.initial.load_active = {0};
    pr.phases = std::move(phases);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { check(false, "meshed: " + M.message); return {}; }
    return katai::app::solve_phases(pr, M.mesh, ip);
}

bool all_ok(const std::vector<SolveResult>& res, size_t n) {
    if (res.size() != n) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return false;
    }
    for (const auto& r : res)
        if (!r.ok) { std::printf("      (%s)\n", r.message.c_str()); return false; }
    return true;
}

std::vector<double> wall_moments(const SolveResult& r) {
    std::vector<double> M;
    for (const auto& f : r.struct_forces)
        if (f.name == "Line")
            for (const auto& st : f.stations) M.push_back(st.M);
    return M;
}

double peak(const std::vector<double>& v) {
    double p = 0.0;
    for (double x : v) p = std::max(p, std::fabs(x));
    return p;
}

// Largest |a - b| over the first `nodes` nodes of two displacement fields.
double field_gap(const SolveResult& a, const SolveResult& b, int nodes) {
    double g = 0.0;
    for (int i = 0; i < 2 * nodes && i < a.disp.size() && i < b.disp.size(); ++i)
        g = std::max(g, std::fabs(a.disp[i] - b.disp[i]));
    return g;
}

void inactive_everywhere_is_no_split() {
    std::printf("\n(a) an element inactive in every phase is the model without it\n");
    using K = m::StructKind;
    const auto grid = solve(block(m::SoilModel::MohrCoulomb, K::Geogrid, false), InitialPhase::K0Procedure,
                            {0, 0}, {phase("Surcharge", m::PhaseType::Plastic, {0, 0}, {1})});
    check(all_ok(grid, 2), "the inactive geogrid line solves");
    if (!all_ok(grid, 2)) return;
    const int n = grid[1].mesh.node_count;
    for (const auto& [kind, ifc, what] :
         {std::tuple{K::Plate, true, "a wall with interfaces"}, std::tuple{K::Interface, false, "an interface"}}) {
        const auto res = solve(block(m::SoilModel::MohrCoulomb, kind, ifc), InitialPhase::K0Procedure, {0, 0},
                               {phase("Surcharge", m::PhaseType::Plastic, {0, 0}, {1})});
        check(all_ok(res, 2), std::string(what) + " inactive in every phase: the chain solves");
        if (!all_ok(res, 2)) continue;
        const double gap = field_gap(res[1], grid[1], n);
        std::printf("      %s: %d nodes (split) against %d, largest gap at the original nodes %.3e m, "
                    "max|u| %.12e / %.12e\n",
                    what, res[1].mesh.node_count, n, gap, res[1].max_disp, grid[1].max_disp);
        check(res[1].mesh.node_count > n, std::string(what) + ": the mesh is split although it is inactive");
        check(gap <= 1e-12, std::string(what) + ": the displacement field is the unsplit model's");
    }
}

void activated_after_k0_equals_present_from_k0() {
    std::printf("\n(b) K0, then the wall activated == the wall present from the K0 phase\n");
    for (double w : {0.0, 9.6}) {
        const auto late = solve(block(m::SoilModel::MohrCoulomb, m::StructKind::Plate, true, w), InitialPhase::K0Procedure,
                                {0, 0}, {phase("Wall", m::PhaseType::Plastic, {1, 0}, {0}),
                                         phase("Surcharge", m::PhaseType::Plastic, {1, 0}, {1})});
        const auto early = solve(block(m::SoilModel::MohrCoulomb, m::StructKind::Plate, true, w), InitialPhase::K0Procedure,
                                 {1, 0}, {phase("Nil", m::PhaseType::Plastic, {1, 0}, {0}),
                                          phase("Surcharge", m::PhaseType::Plastic, {1, 0}, {1})});
        const std::string tag = "w = " + std::to_string(w).substr(0, 3);
        check(all_ok(late, 3) && all_ok(early, 3), tag + ": both chains solve");
        if (!all_ok(late, 3) || !all_ok(early, 3)) continue;
        const auto Ml = wall_moments(late[2]), Me = wall_moments(early[2]);
        double mgap = Ml.size() == Me.size() && !Ml.empty() ? 0.0 : 1e300;
        for (size_t k = 0; k < Ml.size() && k < Me.size(); ++k) mgap = std::max(mgap, std::fabs(Ml[k] - Me[k]));
        const double ugap = field_gap(late[2], early[2], late[2].mesh.node_count);
        const double settle_late = late[1].max_disp, settle_early = w > 0.0 ? early[0].max_disp : 0.0;
        std::printf("      %s: activation phase |u| %.6e (the wall present from K0 settles %.6e there); "
                    "surcharge phase gaps %.3e m, %.3e of max|M| %.6f\n",
                    tag.c_str(), settle_late, settle_early, ugap, mgap / std::max(peak(Me), 1e-30), peak(Me));
        check(peak(Me) > 1.0, tag + ": the surcharge bends the wall (teeth)");
        check(ugap <= 1e-12 && mgap <= 1e-9 * peak(Me),
              tag + ": the surcharge phase is the same whether the wall came in with K0 or after it");
        check(std::fabs(settle_late - settle_early) <= 1e-12,
              tag + ": activating the wall does what its presence in the K0 phase did (its own weight)");
    }
}

void activated_after_a_non_geostatic_stage() {
    std::printf("\n(c) the wall activated after a stage that is not geostatic\n");
    for (const auto model : {m::SoilModel::MohrCoulomb, m::SoilModel::LinearElastic}) {
        const auto res = solve(block(model, m::StructKind::Plate, true), InitialPhase::GravityLoading, {0, 0},
                               {phase("Surcharge", m::PhaseType::Plastic, {0, 0}, {1}),
                                phase("Wall", m::PhaseType::Plastic, {1, 0}, {1})});
        const bool mc = model == m::SoilModel::MohrCoulomb;
        const std::string tag = mc ? "Mohr-Coulomb" : "linear elastic";
        check(all_ok(res, 3), tag + ": gravity -> surcharge -> wall activated solves");
        if (!all_ok(res, 3)) continue;
        const double ratio = res[2].max_disp / res[1].max_disp;
        const double M = peak(wall_moments(res[2]));
        std::printf("      %s: activation |u| %.4e against the surcharge phase's %.4e (%.2e), wall max|M| %.5f\n",
                    tag.c_str(), res[2].max_disp, res[1].max_disp, ratio, M);
        if (mc) {
            check(ratio <= 1e-3, tag + ": the activation moves the ground by less than 1e-3 of the stage before");
            check(M <= 0.05, tag + ": and the wall takes almost no moment from its own installation");
        }
    }
}

void refusals() {
    std::printf("\n(d) what is refused\n");
    const auto cons = solve(block(m::SoilModel::MohrCoulomb, m::StructKind::Plate, true), InitialPhase::K0Procedure,
                            {0, 0}, {phase("Consolidate", m::PhaseType::Consolidation, {1, 0}, {0})});
    check(cons.size() == 2 && cons[0].ok && !cons[1].ok &&
              cons[1].message.find("Plastic (staged construction) phase only") != std::string::npos,
          "activating a wall in a consolidation phase is refused, with the remedy");
    // The wall carries the surcharge with its interfaces; the anchor at its head is carried on while
    // the wall is deactivated -- the joint's two sides have moved apart under the anchor.
    const auto off = solve(block(m::SoilModel::MohrCoulomb, m::StructKind::Plate, true), InitialPhase::K0Procedure,
                           {1, 1}, {phase("Surcharge", m::PhaseType::Plastic, {1, 1}, {1}),
                                    phase("Wall off", m::PhaseType::Plastic, {0, 1}, {1})});
    check(off.size() == 3 && off[1].ok && !off[2].ok &&
              off[2].message.find("cannot be continued from both") != std::string::npos,
          "deactivating the wall under a carried anchor is refused");
    const auto both = solve(block(m::SoilModel::MohrCoulomb, m::StructKind::Plate, true), InitialPhase::K0Procedure,
                            {1, 1}, {phase("Surcharge", m::PhaseType::Plastic, {1, 1}, {1}),
                                     phase("Wall off", m::PhaseType::Plastic, {0, 0}, {1})});
    check(all_ok(both, 3), "and deactivating the anchor with it runs");
}

}  // namespace

int main() {
    std::printf("KV-STR-013: walls and interfaces activated and deactivated per phase\n");
    inactive_everywhere_is_no_split();
    activated_after_k0_equals_present_from_k0();
    activated_after_a_non_geostatic_stage();
    refusals();
    if (g_failures == 0) {
        std::printf("\nOK: walls and interfaces follow the phase's activation\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
