// A joint's strength is that of the ground where the joint stands: R times the strength at its depth, in NGI-ADP on its
// orientation and sense as well, and from the side that has a strength whichever way the line was drawn. And a fixity
// on a line inside the model holds the whole line.
//
// verify: KV-STR-036
//   oracle:   closed_form
//   source:   Grimstad, G., Andresen, L. & Jostad, H.P. (2012). NGI-ADP: anisotropic shear strength model for clay. Int. J. Numer. Anal. Methods Geomech. 36(4), 483-497 -- the traction criterion of the model on a plane, with no volume change in the plane; the Coulomb limit of a joint for the block sliding on it
//   locator:  (a) the criterion by itself, suA 60, suP 24, suDSS 42, R 0.5, on lines at 0, 90 and +/-45 degrees, for both senses; (b) a rigid block 4 m wide sliding on NGI-ADP clay over a joint of R 0.5, drawn both ways, 6-noded, 0.5 m; (c) the same block on Tresca clay of c 20 and 5 per metre below a reference 10 m above the joint, R 0.5; (d) an elastic column 1 x 2 m fixed vertically along its top, a line inside the model, under a block that is then removed
//   quantity: (a) the limits [kPa]; (b), (c) the horizontal force on the block per metre of its width at 150 mm of slip [kPa]; (d) the vertical displacement of the line [m]
//   expected: (a) R suDSS on the level and vertical lines either sense, R suP and R suA at 45 degrees by sense; (b) R suDSS = 21, both ways of drawing; (c) R (20 + 5 x 10) = 35; (d) zero at every node of the line
//   band:     (a) 1e-12 relative; (b), (c) 1% (the elastic slip of the last increment is in it); (d) exact
#include <katai/fem/elements/interface.hpp>
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
namespace iface = katai::core::iface;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

constexpr double kR = 0.5, kSuA = 60.0, kSuP = 24.0, kSuD = 42.0;

void criterion() {
    std::printf("-- (a) the traction criterion --\n");
    double cp = 0.0, cn = 0.0;
    iface::ngi_adp_joint_strength(kR, kSuA, kSuP, kSuD, 1.0, 0.0, cp, cn);
    check(rel(cp, kR * kSuD) < 1e-12 && rel(cn, kR * kSuD) < 1e-12, "a level line: R suDSS, either sense");
    iface::ngi_adp_joint_strength(kR, kSuA, kSuP, kSuD, 0.0, 1.0, cp, cn);
    check(rel(cp, kR * kSuD) < 1e-12 && rel(cn, kR * kSuD) < 1e-12, "a vertical line: R suDSS, either sense");
    const double r2 = std::sqrt(0.5);
    iface::ngi_adp_joint_strength(kR, kSuA, kSuP, kSuD, r2, r2, cp, cn);
    check(rel(cp, kR * kSuP) < 1e-12, "45 degrees, the upper body going up the slope: R suP");
    check(rel(cn, kR * kSuA) < 1e-12, "and coming down it: R suA");
    double cp2 = 0.0, cn2 = 0.0;
    iface::ngi_adp_joint_strength(kR, kSuA, kSuP, kSuD, -r2, -r2, cp2, cn2);
    check(rel(cp2, cp) < 1e-12 && rel(cn2, cn) < 1e-12, "whichever way the line was drawn");
    iface::ngi_adp_joint_strength(kR, kSuA, kSuP, kSuD, r2, -r2, cp, cn);
    check(rel(cp, kR * kSuA) < 1e-12 && rel(cn, kR * kSuP) < 1e-12, "and on the mirror slope the senses swap");
}

void stress_level() {
    std::printf("-- (a2) the stiffness follows the stress the joint was installed with --\n");
    iface::InterfaceProps p;
    p.kn = 1000.0; p.ks = 100.0; p.c_i = 1e6; p.k_power = 0.5; p.k_pref = 100.0;
    check(rel(iface::stiffness_factor(p, -100.0), 1.0) < 1e-12, "at the reference stress: the stiffness at that stress");
    check(rel(iface::stiffness_factor(p, -400.0), 2.0) < 1e-12, "at four times that, by the power 0.5: twice");
    check(rel(iface::stiffness_factor(p, -1.0), std::sqrt(0.1)) < 1e-12, "at next to nothing: the floor, a tenth of the reference stress");
    check(iface::stiffness_factor(p, 0.0) == 1.0, "a new contact, with no stress to read: the reference stiffness");
    const auto r = iface::coulomb_return(p, 1e-3, -1e-3, 0.0, -400.0);
    check(rel(r.tau, 2.0 * 100.0 * 1e-3) < 1e-12 && rel(r.sigma_n, -400.0 - 2.0 * 1000.0 * 1e-3) < 1e-12 && rel(r.Ds, 200.0) < 1e-12 && rel(r.Dn, 2000.0) < 1e-12,
          "and the law reads it, in its stresses and in its tangents");
    p.k_power = 0.0;
    check(iface::stiffness_factor(p, -400.0) == 1.0, "with no power: constant");
}

// The block is drawn on the ground, the joint between them; returns the horizontal force on the block per metre.
double block_on(const m::Material& ground, bool line_left_to_right) {
    constexpr double kB0 = 8.0, kB1 = 12.0, kY = 8.0;
    m::Project pr;
    pr.name = "block on a joint";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = kY + 0.5;
    pr.materials.push_back(ground);
    m::Material block; block.name = "Block"; block.model = m::SoilModel::LinearElastic; block.E = 3e7; block.nu = 0.2;
    block.gamma_unsat = block.gamma_sat = 0.0;
    pr.materials.push_back(block);
    const int F = (int)m::BCType::Free, FF = (int)m::BCType::FullyFixed, H = (int)m::BCType::HorizontallyFixed;
    m::SoilPolygon soil; soil.name = "Clay"; soil.material = 0;
    soil.x = {0, 20, 20, kB1, kB0, 0}; soil.y = {0, 0, kY, kY, kY, kY};
    soil.edge_bc = {FF, H, F, F, F, H};
    pr.polygons.push_back(soil);
    m::SoilPolygon top; top.name = "Block"; top.material = 1;
    top.x = {kB0, kB1, kB1, kB0}; top.y = {kY, kY, kY + 0.5, kY + 0.5};
    top.edge_bc = {F, F, F, F};
    pr.polygons.push_back(top);
    m::StructElement joint; joint.kind = m::StructKind::Interface; joint.name = "Base";
    joint.x1 = line_left_to_right ? kB0 : kB1; joint.y1 = kY;
    joint.x2 = line_left_to_right ? kB1 : kB0; joint.y2 = kY;
    pr.structs.push_back(joint);
    m::PrescribedDisp d; d.name = "Push"; d.x1 = kB0; d.y1 = kY + 0.5; d.x2 = kB1; d.y2 = kY + 0.5;
    d.set_ux = true; d.ux = 0.15; d.set_uy = false;
    pr.disps.push_back(d);
    // Held down, as a joint with no normal stress opens and carries nothing: the strength read is that of a joint in contact.
    m::Load v; v.kind = m::LoadKind::Distributed; v.name = "V"; v.x1 = kB0; v.y1 = kY + 0.5; v.x2 = kB1; v.y2 = kY + 0.5;
    v.qx1 = v.qx2 = 0.0; v.qy1 = v.qy2 = -50.0;
    pr.loads.push_back(v);
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;
    pr.initial.disp_active = {0};
    pr.initial.load_active = {0};
    m::Phase hold; hold.name = "Hold"; hold.type = m::PhaseType::Plastic; hold.disp_active = {0}; hold.load_active = {1};
    pr.phases.push_back(hold);
    m::Phase ph; ph.name = "Push"; ph.type = m::PhaseType::Plastic; ph.disp_active = {1}; ph.load_active = {1};
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return -1.0;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.empty() || !res.back().ok) return -2.0;
    const auto& r = res.back();
    double hsum = 0.0;
    for (int n = 0; n < r.mesh.node_count; ++n)
        if (std::fabs(r.mesh.y[n] - (kY + 0.5)) < 1e-6) hsum += r.reaction[2 * n];
    return std::fabs(hsum) / (kB1 - kB0);
}

void sliding_block() {
    std::printf("-- (b) a block sliding on NGI-ADP clay, the line drawn both ways --\n");
    m::Material clay; clay.name = "Clay"; clay.model = m::SoilModel::NgiAdp; clay.drainage = m::Drainage::UndrainedC;
    clay.c = kSuA; clay.phi = 0.0; clay.psi = 0.0; clay.nu = 0.495; clay.E = 2.0 * 1.495 * 200.0 * kSuA;
    clay.adp_suP = kSuP / kSuA; clay.adp_suDSS = kSuD / kSuA;
    clay.adp_gfC = 4.0; clay.adp_gfE = 15.0; clay.adp_gfDSS = 7.0;
    clay.k0_auto = false; clay.k0 = 1.0;
    clay.gamma_unsat = clay.gamma_sat = 0.0;
    clay.rinter_rigid = false; clay.Rinter = kR;
    const double a = block_on(clay, true), b = block_on(clay, false);
    std::printf("      H/B %.3f and %.3f against R suDSS %.3f\n", a, b, kR * kSuD);
    check(a > 0.0 && rel(a, kR * kSuD) < 0.01, "drawn left to right: R suDSS");
    check(b > 0.0 && rel(b, kR * kSuD) < 0.01, "drawn right to left: the same, not the elastic block's boxes");

    std::printf("-- (c) the same block on Tresca clay stronger with depth --\n");
    m::Material inc; inc.name = "Clay"; inc.model = m::SoilModel::MohrCoulomb; inc.drainage = m::Drainage::UndrainedB;
    inc.E = 2.0e4; inc.nu = 0.3; inc.c = 20.0; inc.phi = 0.0; inc.psi = 0.0; inc.c_inc = 5.0; inc.y_ref = 18.0;
    inc.gamma_unsat = inc.gamma_sat = 0.0;
    inc.rinter_rigid = false; inc.Rinter = kR;
    const double c = block_on(inc, true);
    std::printf("      H/B %.3f against R (20 + 5 x 10) = %.3f\n", c, kR * 70.0);
    check(c > 0.0 && rel(c, kR * 70.0) < 0.01, "the joint reads the strength at its own depth");
}

void fixed_inner_line() {
    std::printf("-- (d) a fixity on a line inside the model --\n");
    m::Project pr;
    pr.name = "column under a block";
    pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 3;
    m::Material soil; soil.name = "Soil"; soil.model = m::SoilModel::LinearElastic; soil.E = 1000.0; soil.nu = 0.3;
    soil.gamma_unsat = soil.gamma_sat = 0.0;
    pr.materials.push_back(soil);
    m::Material dead = soil; dead.name = "Block"; dead.E = 1e5; dead.gamma_unsat = dead.gamma_sat = 100.0;
    pr.materials.push_back(dead);
    const int F = (int)m::BCType::Free, FF = (int)m::BCType::FullyFixed, H = (int)m::BCType::HorizontallyFixed,
              V = (int)m::BCType::VerticallyFixed;
    m::SoilPolygon s; s.name = "S"; s.material = 0; s.x = {0, 1, 1, 0}; s.y = {0, 0, 2, 2}; s.edge_bc = {FF, H, V, H};
    pr.polygons.push_back(s);
    m::SoilPolygon b; b.name = "B"; b.material = 1; b.x = {0, 1, 1, 0}; b.y = {2, 2, 3, 3}; b.edge_bc = {F, H, F, H};
    pr.polygons.push_back(b);
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;
    m::Phase ph; ph.name = "Remove"; ph.type = m::PhaseType::Plastic; ph.poly_active = {1, 0};
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "the model meshes");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(!res.empty() && res.back().ok, "the phase solves");
    if (res.empty() || !res.back().ok) return;
    const auto& r = res.back();
    double worst = 0.0;
    int on_line = 0;
    for (int n = 0; n < r.mesh.node_count; ++n)
        if (std::fabs(r.mesh.y[n] - 2.0) < 1e-6) {
            worst = std::max(worst, std::fabs(r.disp[2 * n + 1]));
            ++on_line;
        }
    std::printf("      %d nodes on the line, largest movement %.3e m\n", on_line, worst);
    check(on_line > 2 && worst == 0.0, "every node of the line stays");
}
}  // namespace

int main() {
    std::printf("== a joint's strength where it stands, a fixity on an inner line (KV-STR-036) ==\n");
    criterion();
    stress_level();
    sliding_block();
    fixed_inner_line();
    if (g_failures) {
        std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
        return 1;
    }
    std::printf("OK: the joint's strength verified\n");
    return 0;
}
