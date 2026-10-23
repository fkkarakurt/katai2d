// test_axisym_interface -- an interface in an axisymmetric model: a pile shaft.
//
// A pile load test is modelled in axisymmetry as a soil region for the pile and an interface
// along its shaft (and base). Until 2026-10-22 every structural element was refused in an
// axisymmetric model, the joint included. A joint is the same element there, integrated per
// radian: its traction acts on a ring of radius r, so each Newton-Cotes point (a node of the
// joint) weighs r.
//
// Oracle (closed form): a weightless cylinder of radius a and length H, free at its base, inside
// weightless ground, joined to it along r = a by an interface with cohesion c_i and no friction,
// is pulled up far enough for the whole joint to slip. The force it then carries per radian is
//     F = c_i a H,
// the shear strength times the lateral area of the shaft divided by 2 pi. A joint integrated as in
// plane strain would carry c_i H, a factor of 1/a away; a = 0.5 m makes that a factor of two.
// The elastic part of the joint before it slips is checked against its own stiffness too: at a
// pull small enough for the joint to stay elastic, an elastic cylinder of axial stiffness E a^2/2
// per radian on a joint of shear stiffness k_s carries F = k_s u a H tanh(lambda H)/(lambda H),
// lambda^2 = 2 k_s/(E a) -- the load-transfer (t-z) solution of an elastic pile.
//
// verify: KV-STR-016
//   oracle:   closed_form
//   source:   the shaft capacity of a cylinder in a cohesive joint, F = 2 pi a H c_i (per radian c_i a H), and the elastic joint's shear force k_s u (2 pi a H); KATAI 2D input contract (docs/k2d-format.md, axisymmetric, structs[i].iface_ks)
//   locator:  the vertical reaction on the prescribed displacement of the cylinder's top, summed over its nodes (per radian in an axisymmetric model)
//   quantity: pull-out force per radian [kN/rad] at full slip and in the elastic range
//   expected: c_i a H within 0.5% at full slip, on 6- and 15-noded elements; the load-transfer closed form within 1% in the elastic range
//   band:     0.5% / 1% -- measured 0.000% at full slip; -0.3% in the elastic range, where the ground's own shear compliance is left out of the closed form

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

struct Pull { bool ok = false; double F = 0.0; std::string msg; };

static Pull pull(int order, double uy, double ci, double ks) {
    constexpr double a = 0.5, H = 4.0, R = 4.0, E = 1.0e7, nu = 0.2;
    m::Project pr; pr.axisymmetric = true;
    m::Material soil; soil.name = "Ground"; soil.model = m::SoilModel::LinearElastic;
    soil.E = E; soil.nu = nu; soil.gamma_unsat = soil.gamma_sat = 0.0; soil.e_init = 0.5;
    m::Material pile = soil; pile.name = "Pile";
    m::Material joint = soil; joint.name = "Joint"; joint.model = m::SoilModel::MohrCoulomb;
    joint.c = ci; joint.phi = 0.0; joint.psi = 0.0; joint.rinter_rigid = false; joint.Rinter = 1.0;
    joint.tension_cutoff = false;
    pr.materials = {soil, pile, joint};
    m::SoilPolygon P; P.name = "Pile"; P.material = 1;
    P.x = {0, a, a, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    m::SoilPolygon G; G.name = "Ground"; G.material = 0;
    G.x = {a, R, R, a}; G.y = {0, 0, H, H};
    // The ground is held at its outer radius; its base is on rollers, so the corner where the
    // shaft meets the base takes no vertical fixity that the cylinder's own free base would share.
    G.edge_bc = {(int)m::BCType::HorizontallyFixed, (int)m::BCType::FullyFixed, (int)m::BCType::Free, (int)m::BCType::Free};
    pr.polygons = {P, G};
    pr.has_water = false;
    m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Shaft";
    j.x1 = a; j.y1 = 0; j.x2 = a; j.y2 = H; j.iface_material = 2; j.iface_ks = ks; j.iface_kn = 10.0 * ks;
    pr.structs.push_back(j);
    // The pull acts on the cylinder's top and stops short of the joint, whose two sides share the
    // corner's coordinates: a line through the corner would pull the ground's node too.
    m::PrescribedDisp d; d.name = "Pull"; d.x1 = 0; d.y1 = H; d.x2 = a - 0.01; d.y2 = H;
    d.set_ux = false; d.set_uy = true; d.uy = uy;
    pr.disps.push_back(d);
    pr.initial.disp_active = {0};
    m::Phase ph; ph.name = "Pull"; ph.type = m::PhaseType::Plastic; ph.disp_active = {1}; ph.load_steps = 20;
    pr.phases.push_back(ph);
    Pull out;
    const auto M = katai::app::mesh_from_project(pr, 0.25, order);
    if (!M.ok) { out.msg = M.message; return out; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2 || !res[1].ok) { out.msg = res.empty() ? "no result" : res.back().message; return out; }
    const auto& Rz = res[1];
    if (Rz.reaction.size() != 2 * Rz.mesh.node_count) { out.msg = "no reactions"; return out; }
    for (int n = 0; n < Rz.mesh.node_count; ++n)
        if (std::fabs(Rz.mesh.y[n] - H) < 1e-9 && Rz.mesh.x[n] <= a + 1e-9) out.F += Rz.reaction[2 * n + 1];
    out.F = std::fabs(out.F);
    out.ok = true;
    return out;
}

int main() {
    std::printf("An interface in an axisymmetric model: a cylinder pulled through a cohesive joint\n\n");
    constexpr double a = 0.5, H = 4.0, ci = 20.0, ks = 1.0e4;
    for (int order : {6, 15}) {
        const std::string tag = std::to_string(order) + "-noded";
        const Pull full = pull(order, 0.05, ci, ks);
        if (!full.ok) { check(false, tag + ": the pull-out solves (" + full.msg + ")"); continue; }
        const double want = ci * a * H;
        std::printf("  %s: force at full slip %.4f kN/rad, closed form c_i a H = %.4f (plane strain would be %.4f)\n",
                    tag.c_str(), full.F, want, ci * H);
        check(std::fabs(full.F - want) / want < 0.005, tag + ": the shaft carries c_i a H per radian at full slip");
        const double u_el = 0.2 * ci / ks;   // a fifth of the slip displacement
        const Pull el = pull(order, u_el, ci, ks);
        if (!el.ok) { check(false, tag + ": the elastic pull solves (" + el.msg + ")"); continue; }
        // An elastic cylinder (axial stiffness E a^2/2 per radian) on an elastic joint (k_s a per
        // unit length per radian), pulled at its top: the load-transfer closed form
        // F = k_s u a H tanh(lambda H) / (lambda H), lambda^2 = 2 k_s / (E a).
        constexpr double Ep = 1.0e7;
        const double lam = std::sqrt(2.0 * ks / (Ep * a));
        const double want_el = ks * u_el * a * H * std::tanh(lam * H) / (lam * H);
        std::printf("  %s: elastic force %.4f kN/rad, load-transfer closed form %.4f (rigid cylinder %.4f)\n",
                    tag.c_str(), el.F, want_el, ks * u_el * a * H);
        check(std::fabs(el.F - want_el) / want_el < 0.01, tag + ": in the elastic range the joint carries the load-transfer closed form");
    }
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
