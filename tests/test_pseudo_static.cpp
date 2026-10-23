// PSEUDO-STATIC ACCELERATION: a phase's horizontal and vertical seismic coefficients as body forces.
//
// The equivalent static method of a seismic slope or wall check applies the inertia of the ground
// as a body force: kh W horizontally and kv W vertically, W the TOTAL weight (the pore water moves
// with the soil), kv positive upward -- it takes kv W off the weight, as (1 - kv) does in TBDY 2018
// Eq. 16.31. Each Plastic phase states its own; a Safety phase after it reduces the strength under
// its parent's.
//
// Oracles (closed form):
//   (1) the support reactions under kh, kv: the base pushes back with -kh W horizontally and
//       (1 - kv) W vertically, W the weight above it. The reported reaction leaves out the share of
//       the weight that lands on the supported nodes themselves -- in both components alike -- so
//       the RATIO -kh : (1 - kv) is exact;
//   (2) a block resting on a frictional joint (c_i, phi_i) on a horizontal plane, stiff enough to
//       fail all at once, pushed by kh and lightened by kv: the factor by which the joint's strength
//       has to be divided to make it slide,
//           FS = (c_i B + (1 - kv) W tan phi_i) / (kh W);
//   (3) the same block on a plane inclined at beta, kh acting down the slope:
//           FS = (c_i B + W [(1 - kv) cos beta - kh sin beta] tan phi_i)
//                / (W [(1 - kv) sin beta + kh cos beta]),
//       the planar case of TBDY 2018 Eq. 16.31 with the factor on both c and tan phi (16.13.11);
//   (4) the code's verdict on that factor, gamma_Rk >= 1.10 (16.13.11 c), on either side of it.
//
// verify: KV-SLP-005
//   oracle:   closed_form
//   source:   limit equilibrium of a rigid block on a planar Mohr-Coulomb joint under horizontal and vertical body forces kh W and kv W (the pseudo-static method; Kramer, S.L. (1996). Geotechnical Earthquake Engineering. Prentice Hall, sec. 10.6.1.1); TBDY 2018 (Resmi Gazete 30364 bis, 18.03.2018) 16.13.9 Eq. 16.28 for the inertia forces and 16.13.11 Eq. 16.31 for the factor on c' and tan phi'
//   locator:  a stiff 4 x 1 m block (gamma = 20) on a joint over a fixed base, started by gravity loading: (1) and (2) on a horizontal joint, phi_i = 30 deg, c_i = 0 and 5 kPa, kh = 0.2, kv = 0, +0.1 and -0.1; (3) on a plane at 10 deg, kh = 0.15; 6- and 15-noded elements
//   quantity: the ratio of the horizontal to the vertical base reaction [-]; factor of safety [-]
//   expected: (1) -kh / (1 - kv); (2) (c_i B + (1 - kv) W tan phi_i) / (kh W) = 2.8868, 2.5981, 4.4254; (3) the inclined form above = 1.7224; (4) met at kh = 0.50 (1.155), not met at kh = 0.55 (1.050)
//   band:     (1) 1e-9 relative (the body force is the weight integrated in another direction); (2), (3) 1% (the strength-reduction search's resolution on a block that fails at once; measured 0.00% on both element orders)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
constexpr double kPi = 3.14159265358979323846;

// A stiff block (B x 1 m, gamma = 20) on a joint over a fixed base, the whole model turned by beta
// (the plane of the joint inclined at beta, rising to the left so that +x along the plane is down
// the slope). Phases: K0 (or gravity for an inclined plane), a Plastic phase with kh, kv, then a
// Safety phase. Returns the factor of safety and, from the Plastic phase, the base reactions.
struct Out {
    bool ok = false;
    double fos = 0.0, Rx = 0.0, Ry = 0.0, W = 0.0;
    std::string msg;
    int verdict = 0;   // K2D-A039 in the Safety phase: 1 note (met), -1 warning (not met), 0 none
};
Out block(double beta_deg, double c_i, double phi_i, double kh, double kv, int order,
          bool tbdy = false) {
    constexpr double B = 4.0, hb = 1.0, hk = 1.0, gamma = 20.0;
    const double b = beta_deg * kPi / 180.0, cb = std::cos(b), sb = std::sin(b);
    // the four corners of a rectangle (x along the plane, y normal to it), turned by -beta about
    // the origin: the plane descends towards +x
    const auto X = [&](double x, double y) { return x * cb + y * sb; };
    const auto Y = [&](double x, double y) { return -x * sb + y * cb; };
    m::Project pr;
    m::Material base; base.name = "Base"; base.model = m::SoilModel::LinearElastic;
    base.E = 2.0e7; base.nu = 0.2; base.gamma_unsat = base.gamma_sat = gamma; base.e_init = 0.5;
    m::Material joint = base; joint.name = "Joint"; joint.model = m::SoilModel::MohrCoulomb;
    joint.c = c_i; joint.phi = phi_i; joint.psi = 0.0; joint.rinter_rigid = false; joint.Rinter = 1.0;
    m::Material blk = base; blk.name = "Block";
    pr.materials = {base, blk, joint};
    m::SoilPolygon Pb; Pb.name = "Base"; Pb.material = 0;
    Pb.x = {X(0, 0), X(B, 0), X(B, hb), X(0, hb)}; Pb.y = {Y(0, 0), Y(B, 0), Y(B, hb), Y(0, hb)};
    Pb.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free};
    m::SoilPolygon Pk; Pk.name = "Block"; Pk.material = 1;
    Pk.x = {X(0, hb), X(B, hb), X(B, hb + hk), X(0, hb + hk)};
    Pk.y = {Y(0, hb), Y(B, hb), Y(B, hb + hk), Y(0, hb + hk)};
    Pk.edge_bc.assign(4, (int)m::BCType::Free);
    pr.polygons = {Pb, Pk};
    pr.has_water = false;
    m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Joint";
    j.x1 = X(0, hb); j.y1 = Y(0, hb); j.x2 = X(B, hb); j.y2 = Y(B, hb); j.iface_material = 2;
    pr.structs.push_back(j);
    m::Phase p1; p1.type = m::PhaseType::Plastic; p1.name = "Shaken";
    p1.pseudo_kh = kh; p1.pseudo_kv = kv;
    m::Phase p2; p2.type = m::PhaseType::Safety; p2.name = "Safety";
    if (tbdy) p2.design_approach = m::DesignApproach::TBDY2018_Seismic;
    pr.phases = {p1, p2};
    Out o;
    o.W = gamma * B * hk;
    const auto M = katai::app::mesh_from_project(pr, 0.25, order);
    if (!M.ok) { o.msg = M.message; return o; }
    const auto res = katai::app::solve_phases(pr, M.mesh, beta_deg == 0.0 ? InitialPhase::K0Procedure
                                                                           : InitialPhase::GravityLoading);
    if (res.size() != 3 || !res[1].ok || !res[2].ok) {
        o.msg = res.empty() ? "no result" : res.back().message;
        return o;
    }
    o.ok = true;
    o.fos = res[2].fos;
    for (const auto& d : res[2].diagnostics)
        if (d.code == "K2D-A039")
            o.verdict = d.severity == katai::core::DiagnosticSeverity::Warning ? -1 : 1;
    // The base's reactions in the Plastic phase (the soil's share; the weight lumped on the
    // supported nodes themselves is in neither component, which is why the RATIO is the identity).
    const auto& Rc = res[1].reaction;
    for (int n = 0; 2 * n + 1 < Rc.size(); ++n) { o.Rx += Rc[2 * n]; o.Ry += Rc[2 * n + 1]; }
    return o;
}
}  // namespace

int main() {
    std::printf("Pseudo-static acceleration: kh W and kv W as body forces\n\n");
    const double tan30 = std::tan(30.0 * kPi / 180.0);
    for (int order : {6, 15}) {
        std::printf("-- %d-noded --\n", order);
        // (2) horizontal plane
        struct Case { double c, kh, kv; };
        for (const Case cs : {Case{0.0, 0.2, 0.0}, Case{0.0, 0.2, 0.1}, Case{5.0, 0.2, -0.1}}) {
            const Out o = block(0.0, cs.c, 30.0, cs.kh, cs.kv, order);
            const double want = (cs.c * 4.0 + (1.0 - cs.kv) * o.W * tan30) / (cs.kh * o.W);
            char label[160];
            std::snprintf(label, sizeof(label), "horizontal plane, c_i = %.0f, kh = %.2f, kv = %+.2f", cs.c, cs.kh, cs.kv);
            if (!o.ok) { check(false, std::string(label) + ": solved (" + o.msg + ")"); continue; }
            std::printf("   %s: FoS %.4f, closed form %.4f (%+.2f%%); base reaction ratio %.6f, -kh/(1-kv) %.6f\n",
                        label, o.fos, want, 100.0 * (o.fos - want) / want, o.Rx / o.Ry, -cs.kh / (1.0 - cs.kv));
            check(std::fabs(o.fos - want) / want < 0.01, std::string(label) + ": factor of safety = closed form within 1%");
            // The base pushes back: against the inertia toward +x, and up under the weight.
            check(std::fabs(o.Rx / o.Ry + cs.kh / (1.0 - cs.kv)) < 1e-9 * (cs.kh / (1.0 - cs.kv)),
                  std::string(label) + ": base reactions in the ratio -kh : (1 - kv)");
        }
        // (3) inclined plane, kh down the slope
        {
            const double beta = 10.0, kh = 0.15, kv = 0.0, c = 0.0;
            const double b = beta * kPi / 180.0;
            const Out o = block(beta, c, 30.0, kh, kv, order);
            const double W = o.W;
            const double want = (c * 4.0 + W * ((1.0 - kv) * std::cos(b) - kh * std::sin(b)) * tan30) /
                                (W * ((1.0 - kv) * std::sin(b) + kh * std::cos(b)));
            if (!o.ok) check(false, "inclined plane: solved (" + o.msg + ")");
            else {
                std::printf("   plane at %.0f deg, kh = %.2f: FoS %.4f, closed form %.4f (%+.2f%%)\n", beta, kh,
                            o.fos, want, 100.0 * (o.fos - want) / want);
                check(std::fabs(o.fos - want) / want < 0.01, "inclined plane: factor of safety = closed form within 1%");
            }
        }
    }
    // (4) the TBDY 2018 verdict (16.13.11 c) on either side of 1.10: kh = 0.50 gives
    //     tan 30 / 0.50 = 1.155, kh = 0.55 gives 1.050
    std::printf("-- the TBDY 2018 verdict --\n");
    for (const double kh : {0.50, 0.55}) {
        const Out o = block(0.0, 0.0, 30.0, kh, 0.0, 6, true);
        const double want = tan30 / kh;
        if (!o.ok) { check(false, "verdict case solved (" + o.msg + ")"); continue; }
        std::printf("   kh = %.2f: FoS %.4f (closed form %.4f), verdict %s\n", kh, o.fos, want,
                    o.verdict > 0 ? "met" : o.verdict < 0 ? "not met" : "none");
        check(want >= 1.10 ? o.verdict == 1 : o.verdict == -1,
              want >= 1.10 ? "gamma_Rk = 1.155 >= 1.10: met (a note)" : "gamma_Rk = 1.050 < 1.10: not met (a warning)");
    }
    if (g_failures == 0) {
        std::printf("\nOK: the pseudo-static body force gives the closed-form factors of safety\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", g_failures);
    return 1;
}
