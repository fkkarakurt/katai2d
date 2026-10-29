// The water at a wall's joint. A joint's law is a law of the soil skeleton: it holds and rubs with the
// EFFECTIVE normal stress, and the pore water goes past it -- back onto the ground face it stands
// against, and onto the wall.
//
// The pore load leaves a term p n on each face of a wall's split seam, as on any boundary of the
// ground. Until 2026-10-29 the joint stood against it with its own normal stress: installed carrying
// the effective stress, it came to the total stress by closing by p / k_n -- the ground moved at
// every wall installed below the water table -- and then rubbed with the total stress, a joint 5 m
// below the water table holding three times the friction its effective stress allows. With one side
// dug away and the water left standing, the exposed seam face took the pit's water on the ground,
// not on the wall.
//
// THE ORACLE IS THE SAME GROUND DRY. Effective stress is what the skeleton feels: the ground under a
// water table at its surface and the same ground dry with the buoyant unit weight gamma' = gamma_sat
// - gamma_w have the same effective stresses everywhere, the same stiffness and strength, and -- with
// water standing in a pit dug in it -- the same displacements: the water against the wall's two
// faces cancels, and the water on the pit floor balances the floor's own pore water.
//
// verify: KV-STR-033
//   oracle:   independent_path
//   source:   Terzaghi's effective stress (the skeleton's strength and stiffness read sigma' = sigma - u); Coulomb friction of a soil-structure interface on the effective normal stress, |tau| <= c_i - sigma'_n tan(phi_i) (Terzaghi 1943; Potts & Zdravkovic 2001, interface elements); KATAI 2D input contract (docs/k2d-format.md, materials[i].Rinter: c_i = R c, tan(phi_i) = R tan(phi); structs[i].iface_pos / iface_neg)
//   locator:  a Mohr-Coulomb block 40 x 20 m (E = 20000 kPa, nu = 0.3, c = 1 kPa, phi = 25 deg, psi = 0, R_inter = 0.5, K0 = 0.5), a 10 m wall (EA = 1e7 kN/m, EI = 1e5 kNm2/m, weightless) at x = 20 with interfaces on both sides, installed in a Plastic phase after the K0 phase; wet: gamma_sat = 18 kN/m3 under a water table at the surface; dry: gamma = 18 - 9.81 and no water. (3): the wet wall pulled up at its head with 1.5 times the friction of its joints; (4): a 4 m pit on the right of the wall, the water left standing in it
//   quantity: the displacement of the installation [m]; the joints' normal stress at every station [kPa]; the shear of every slipping joint point against c_i - sigma'_n tan(phi_i) [kPa]; the wall's largest moment [kNm/m] and displacement [m] after the pit
//   expected: no movement; the wet joints' normal stress = the dry ones'; every joint point slipping between the head and the toe, at its effective friction; the wet wall's moment and displacement = the dry ones'
//   band:     1e-9 m on the installation; 1e-6 relative on the normal stress (the two are the same equations with the water moved across; measured 0); 90% of the points slipping and 1e-3 relative on their shear (measured 128 of 128, 1e-13); 1e-3 relative on the moment and the displacement (measured 0.0000%). Before the fix, on the same model with Jaky's K0 (python/tests/joint_water.py): 0.43 mm of installation, 89 kPa on the normal stress, no point slipping, the moment 10.7% and the displacement 26.7% under the dry ones
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

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kGw = 9.81;
constexpr double kRinter = 0.5, kC = 1.0, kPhi = 25.0;

// pull: an upward point load on the wall's head [kN/m] in a phase after the installation; pit: a
// phase digging the pit with the water left standing in it.
std::vector<katai::app::SolveResult> run(bool wet, double pull, bool pit) {
    m::Project pr;
    pr.x_min = 0; pr.x_max = 40; pr.y_min = 0; pr.y_max = 20;
    m::Material s;
    s.name = "Soil"; s.model = m::SoilModel::MohrCoulomb; s.drainage = m::Drainage::Drained;
    s.E = 2.0e4; s.nu = 0.3; s.c = kC; s.phi = kPhi; s.psi = 0.0;
    s.gamma_sat = 18.0; s.gamma_unsat = wet ? 18.0 : 18.0 - kGw;
    if (!wet) s.gamma_sat = s.gamma_unsat;
    s.rinter_rigid = false; s.Rinter = kRinter;
    s.k0_auto = false; s.k0 = 0.5;
    pr.materials.push_back(s);
    const int F = (int)m::BCType::Free, H = (int)m::BCType::HorizontallyFixed, X = (int)m::BCType::FullyFixed;
    m::SoilPolygon behind; behind.name = "Behind"; behind.material = 0;
    behind.x = {0, 20, 20, 20, 0}; behind.y = {0, 0, 16, 20, 20}; behind.edge_bc = {X, F, F, F, H};
    m::SoilPolygon dig; dig.name = "Pit"; dig.material = 0;
    dig.x = {20, 40, 40, 20}; dig.y = {16, 16, 20, 20}; dig.edge_bc = {F, H, F, F};
    m::SoilPolygon below; below.name = "Below"; below.material = 0;
    below.x = {20, 40, 40, 20}; below.y = {0, 0, 16, 16}; below.edge_bc = {X, H, F, F};
    pr.polygons = {behind, dig, below};
    pr.has_water = wet;
    pr.wx = {0, 40}; pr.wy = {20, 20};
    m::PlateMaterial pm; pm.name = "Wall"; pm.EA = 1.0e7; pm.EI = 1.0e5; pm.w = 0.0;
    pr.plates.push_back(pm);
    m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall";
    w.x1 = 20; w.y1 = 20; w.x2 = 20; w.y2 = 10; w.material = 0;
    w.iface_pos = true; w.iface_neg = true;
    pr.structs.push_back(w);
    m::Load up; up.kind = m::LoadKind::Point; up.name = "Pull";
    up.x1 = 20; up.y1 = 20; up.qx1 = 0.0; up.qy1 = pull;
    pr.loads.push_back(up);
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;   // the identities below hold on any mesh; no refinement needed
    pr.initial.struct_active = {0};
    pr.initial.load_active = {0};
    m::Phase install; install.type = m::PhaseType::Plastic; install.name = "Install";
    install.struct_active = {1}; install.load_active = {0};
    pr.phases = {install};
    if (pull != 0.0) {
        m::Phase p = install; p.name = "Pull"; p.load_active = {1};
        pr.phases.push_back(p);
    }
    if (pit) {
        m::Phase p = install; p.name = "Dig"; p.poly_active = {1, 0, 1};
        pr.phases.push_back(p);
    }
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return {};
    return katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
}

bool all_ok(const std::vector<katai::app::SolveResult>& r, size_t n) {
    if (r.size() != n) return false;
    for (const auto& x : r)
        if (!x.ok) { std::printf("  (%s)\n", x.message.c_str()); return false; }
    return true;
}

}  // namespace

int main() {
    std::printf("== the water at a wall's joint (KV-STR-033) ==\n");
    const double tan_i = kRinter * std::tan(kPhi * 3.14159265358979323846 / 180.0), c_i = kRinter * kC;

    // (1) + (2) the installation
    const auto wet = run(true, 0.0, false), dry = run(false, 0.0, false);
    check(all_ok(wet, 2) && all_ok(dry, 2), "both installations solve");
    if (wet.size() == 2 && dry.size() == 2) {
        std::printf("      installation below the water table: max |u| = %.3e m\n", wet[1].max_disp);
        check(wet[1].max_disp < 1e-9, "installing a wall below the water table moves nothing");
        double dev = 0.0, big = 0.0;
        size_t n = 0;
        bool same = wet[1].interface_forces.size() == dry[1].interface_forces.size();
        for (size_t i = 0; same && i < wet[1].interface_forces.size(); ++i) {
            const auto& a = wet[1].interface_forces[i].stations;
            const auto& b = dry[1].interface_forces[i].stations;
            same = a.size() == b.size();
            for (size_t k = 0; same && k < a.size(); ++k, ++n) {
                dev = std::max(dev, std::fabs(a[k].sigma_n - b[k].sigma_n));
                big = std::max(big, std::fabs(b[k].sigma_n));
            }
        }
        std::printf("      %zu joint stations: max |sigma_n(wet) - sigma_n(dry, gamma')| = %.3e of %.2f kPa\n",
                    n, dev, big);
        check(same && n > 0 && dev < 1e-6 * big,
              "the joint holds the effective normal stress: the dry ground of buoyant weight's");

        // (3) pulled out with 1.5 times the friction of its joints
        double cap = 0.0;
        for (const auto& f : wet[1].interface_forces)
            for (size_t k = 0; k + 1 < f.stations.size(); ++k) {
                const auto& a = f.stations[k];
                const auto& b = f.stations[k + 1];
                const double ta = std::max(0.0, c_i - a.sigma_n * tan_i), tb = std::max(0.0, c_i - b.sigma_n * tan_i);
                cap += 0.5 * (ta + tb) * std::fabs(b.s - a.s);
            }
        const auto pulled = run(true, 1.5 * cap, false);
        int slipping = 0, points = 0;
        double worst = 0.0;
        if (all_ok(pulled, 3))
            for (const auto& f : pulled[2].interface_forces)
                for (const auto& q : f.stations) {
                    if (!(q.y > 10.5 && q.y < 19.5)) continue;
                    ++points;
                    if (!q.slipping) continue;
                    ++slipping;
                    const double lim = c_i - q.sigma_n * tan_i;
                    worst = std::max(worst, std::fabs(std::fabs(q.tau) - lim) / lim);
                }
        std::printf("      friction of the joints on their effective normal stress %.2f kN/m; pulled with 1.5 x: "
                    "%d of %d points slipping, worst |tau| against c_i - sigma'_n tan(phi_i) %.2e\n",
                    cap, slipping, points, worst);
        check(points > 0 && slipping >= 0.9 * points && worst < 1e-3,
              "pulled out, the joint slips with the friction of its effective normal stress");
    }

    // (4) a pit dug below the water table, the water left standing in it
    const auto pw = run(true, 0.0, true), pd = run(false, 0.0, true);
    check(all_ok(pw, 3) && all_ok(pd, 3), "both pits solve");
    if (pw.size() == 3 && pd.size() == 3 && !pw[2].struct_forces.empty() && !pd[2].struct_forces.empty()) {
        double mw = 0.0, md = 0.0, uw = 0.0, ud = 0.0;
        for (const auto& s : pw[2].struct_forces[0].stations) { mw = std::max(mw, std::fabs(s.M)); uw = std::max(uw, std::fabs(s.ux)); }
        for (const auto& s : pd[2].struct_forces[0].stations) { md = std::max(md, std::fabs(s.M)); ud = std::max(ud, std::fabs(s.ux)); }
        std::printf("      wall after the pit: max |M| %.4f (wet) / %.4f (dry, gamma') kNm/m, max |ux| %.4f / %.4f mm\n",
                    mw, md, 1000.0 * uw, 1000.0 * ud);
        check(md > 1.0 && std::fabs(mw / md - 1.0) < 1e-3 && std::fabs(uw / ud - 1.0) < 1e-3,
              "a pit dug below the water table bends the wall as the dry ground of buoyant weight does");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
