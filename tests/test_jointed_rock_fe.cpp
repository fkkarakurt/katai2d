// Jointed rock through the whole program: a specimen, a specimen's factor of safety, a rock slope.
//
// The point (test_jointed_rock, KV-CST-025) answers what stress the model returns. This answers
// whether the program around it -- the file's fields, the catalogue, the anisotropic operator in
// the global stiffness, the strength reduction of every joint set -- delivers that to a model:
//   (a) a plane strain block of bedded rock (bedding horizontal, transversely isotropic) pushed
//       down by its top, elastic: its stiffness is the transversely isotropic body's in plane
//       strain, sigma_yy / eps_yy = 1 / (1/E2 - nu2^2 E1 / E2^2);
//   (b) the same block with its joints dipping beta, pushed past failure, uni-axially (smooth
//       base, one corner held, sides free), with non-associated slip (psi = 0): it carries
//       Jaeger's sigma_1 = 2 c / ((1 - tan phi cot beta) sin 2 beta), and the run says that its
//       joints are non-associated (K2D-M008);
//   (c) the block under (p3, p1), and its factor of safety: the factor F at which Jaeger's stress
//       with c / F and tan(phi) / F is p1;
//   (d) a rock slope, 10 m high at 60 degrees, on one set of cohesionless joints dipping out of
//       the face, associated: every plane of the set that daylights in the face cuts a block whose
//       factor is tan(phi) / tan(dip) (Hoek, E. & Bray, J.W. (1981). Rock Slope Engineering, 3rd
//       ed. IMM, London, ch. 7, without cohesion, water or a tension crack). Cohesionless,
//       because with cohesion the smallest factor is the plane through the toe's, and there the
//       plane runs into the ground in front: a limit a mesh approaches only from above.
//       Associated, because with psi < phi the same slope stops under its weight: at 45% with psi 0
//       and 65% with psi 10 (15-noded, 2 m, phi 45) -- which K2D-M008 says to the user.
//
// verify: KV-CST-026
//   oracle:   closed_form
//   source:   Jaeger, J.C. (1960). Shear failure of anisotropic rocks. Geological Magazine 97(1), 65-72 -- the plane of weakness; Hoek, E. & Bray, J.W. (1981). Rock Slope Engineering, 3rd ed. Institution of Mining and Metallurgy, London -- plane failure of a rock slope; Lekhnitskii, S.G. (1963). Theory of Elasticity of an Anisotropic Elastic Body. Holden-Day -- the transversely isotropic compliance
//   locator:  (a)-(c) a 1 x 1 m plane strain block, tri6 at 0.25 m, smooth base, one corner held horizontally, sides free; E1 = 2 GPa, nu1 = 0.25, E2 = 0.5 GPa, nu2 = 0.2, G2 = 0.3 GPa; joints c = 100 kPa, phi = 30 deg, psi = 0; (a) top displaced 1e-4 m with the bedding horizontal; (b) joints dipping 50 and 65 deg, top displaced 5 mm; (c) joints dipping 55 deg under p3 = 200 kPa and p1 = 700 kPa (Jaeger at full strength: 969 kPa); (d) a 10 m slope at 60 deg on a 70 x 30 m section, gamma = 25 kN/m3, gravity loading, cohesionless associated joints dipping 35 deg (phi = 45) and 30 deg (phi = 40), tri15 at 1 m
//   quantity: (a) the mean sigma_yy over the block [kPa]; (b) the mean sigma_yy at the end [kPa]; (c), (d) the factor of safety [-]
//   expected: (a) the plane strain stiffness of the transversely isotropic body; (b) Jaeger's sigma_1; (c) the F at which the reduced Jaeger stress equals p1; (d) tan(phi) / tan(dip) = 1.4281 and 1.4533
//   band:     (a), (b) 1e-6 relative (a uniform state: only the global equilibrium tolerance); (c) 1% (the search's resolution); (d) 2% (measured -0.00% on the first)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdlib>
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
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }
constexpr double PI = 3.14159265358979323846;
double rad(double d) { return d * PI / 180.0; }

constexpr double E1 = 2.0e6, NU1 = 0.25, E2 = 0.5e6, NU2 = 0.2, G2 = 0.3e6;

m::Material rock(double dip, double c, double phi, double gamma, bool aniso) {
    m::Material r;
    r.name = "Bedded rock"; r.model = m::SoilModel::JointedRock; r.drainage = m::Drainage::Drained;
    r.gamma_unsat = r.gamma_sat = gamma; r.e_init = 0.5;
    r.E = E1; r.nu = NU1;
    r.c = c; r.phi = phi; r.psi = 0.0;
    r.tension_cutoff = true; r.tensile_strength = 0.0;
    r.jr_sets = 1; r.jr_alpha[0] = dip;
    r.jr_aniso = aniso; r.jr_E2 = E2; r.jr_nu2 = NU2; r.jr_G2 = G2;
    r.k0_auto = false; r.k0 = 1.0;
    return r;
}

double jaeger(double s3, double c, double phi, double beta) {
    const double t = std::tan(phi);
    return s3 + 2.0 * (c + s3 * t) / ((1.0 - t / std::tan(beta)) * std::sin(2.0 * beta));
}

// A 1 x 1 m block on a smooth base, its top displaced by uy, held horizontally at one corner only.
// Its sides are free: slip on an inclined joint shears the block, and a side held horizontally
// would forbid the uniform shear that slip needs (the block would be propped, not tested); a top
// held horizontally would forbid the lateral strain. With one corner held -- a zero-length
// prescribed line -- the uniform state with sigma_xx = sigma_xy = 0 is admissible: the block
// shears as a parallelogram over its base.
bool g_said_m008 = false;   // the last push said K2D-M008 (non-associated joints)
std::pair<bool, double> push(const m::Material& mat, double uy) {
    g_said_m008 = false;
    m::Project pr;
    pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 1;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 0.25;
    pr.mesh.auto_refine = false;
    pr.materials.push_back(mat);
    m::SoilPolygon g;
    g.name = "Block"; g.material = 0; g.x = {0, 1, 1, 0}; g.y = {0, 0, 1, 1};
    g.edge_bc = {(int)m::BCType::VerticallyFixed, 0, 0, 0};
    pr.polygons.push_back(g);
    m::PrescribedDisp D;
    D.name = "Top"; D.x1 = 0; D.y1 = 1; D.x2 = 1; D.y2 = 1; D.set_uy = true; D.uy = uy;
    pr.disps.push_back(D);
    m::PrescribedDisp C;
    C.name = "Corner"; C.x1 = C.x2 = 0.0; C.y1 = C.y2 = 0.0; C.set_ux = true; C.ux = 0.0; C.set_uy = false;
    pr.disps.push_back(C);
    pr.initial.disp_active = {0, 1};
    m::Phase ph = pr.initial;
    ph.name = "Push"; ph.type = m::PhaseType::Plastic; ph.disp_active = {1, 1}; ph.load_steps = 20;
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return {false, 0.0};
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    for (const auto& r : res)
        for (const auto& d : r.diagnostics) g_said_m008 = g_said_m008 || d.code == "K2D-M008";
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.substr(0, 300).c_str());
        return {false, 0.0};
    }
    double s = 0.0;
    for (const auto& v : res[1].stress.stress) s += v(1);
    return {true, s / (double)res[1].stress.stress.size()};
}

}  // namespace

int main() {
    std::printf("== Jointed rock through the program (KV-CST-026) ==\n");
    {
        // (a) the anisotropic stiffness, bedding horizontal
        const auto r = push(rock(0.0, 100.0, 30.0, 0.0, true), -1e-4);
        const double want = -1e-4 / (1.0 / E2 - NU2 * NU2 * E1 / (E2 * E2));
        std::printf("      bedded block, elastic: sigma_yy %.6f kPa, transversely isotropic plane strain %.6f kPa\n", r.second, want);
        check(r.first && rel(r.second, want) < 1e-6, "(a) a bedded block is as stiff as the transversely isotropic body in plane strain");
    }
    for (double dip : {50.0, 65.0}) {
        // (b) Jaeger, uni-axially, with non-associated slip (psi = 0): a uniform state, where the
        // point's verification carries over whole -- and the run says the joints are non-associated
        const auto r = push(rock(dip, 100.0, 30.0, 0.0, true), -5e-3);
        const bool said = g_said_m008;
        const double want = -jaeger(0.0, 100.0, rad(30.0), rad(dip));
        std::printf("      joints dipping %.0f deg: sigma_yy %.6f kPa, Jaeger %.6f kPa\n", dip, r.second, want);
        char lbl[160];
        std::snprintf(lbl, sizeof lbl, "(b) joints dipping %.0f deg: the block carries Jaeger's uni-axial strength", dip);
        check(r.first && rel(r.second, want) < 1e-6, lbl);
        std::snprintf(lbl, sizeof lbl, "(b) joints dipping %.0f deg, psi < phi: the run says the joints are non-associated (K2D-M008)", dip);
        check(said, lbl);
    }
    {
        m::Material assoc = rock(50.0, 100.0, 30.0, 0.0, true);
        assoc.psi = assoc.phi;
        (void)push(assoc, -1e-4);
        check(!g_said_m008, "(b) ...and with psi = phi it says nothing of the kind");
    }
    {
        // (c) a block's factor of safety
        constexpr double p3 = 200.0, p1 = 700.0, dip = 55.0, c = 100.0, phi = 30.0;   // Jaeger at full strength: 969 kPa
        double lo = 0.5, hi = 5.0;
        for (int it = 0; it < 200; ++it) {
            const double F = 0.5 * (lo + hi);
            if (jaeger(p3, c / F, std::atan(std::tan(rad(phi)) / F), rad(dip)) > p1) lo = F; else hi = F;
        }
        const double want = 0.5 * (lo + hi);
        m::Project pr;
        pr.materials.push_back(rock(dip, c, phi, 0.0, true));
        m::SoilPolygon B; B.name = "Rock"; B.material = 0;
        B.x = {0, 1, 1, 0}; B.y = {0, 0, 1, 1};
        // the sides free and loaded alike, one corner held horizontally -- here by a very stiff
        // horizontal ground anchor, since neither the K0 procedure nor a Safety phase takes a
        // prescribed displacement; the uniform state leaves that corner where it is, so the
        // anchor carries nothing: the block may shear (above)
        B.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::Free, (int)m::BCType::Free,
                     (int)m::BCType::Free};
        pr.polygons.push_back(B);
        pr.has_water = false;
        m::Load top; top.kind = m::LoadKind::Distributed; top.name = "p1";
        top.x1 = 0; top.y1 = 1; top.x2 = 1; top.y2 = 1; top.qx1 = top.qx2 = 0.0; top.qy1 = top.qy2 = -p1;
        m::Load side; side.kind = m::LoadKind::Distributed; side.name = "p3";
        side.x1 = 1; side.y1 = 0; side.x2 = 1; side.y2 = 1; side.qx1 = side.qx2 = -p3; side.qy1 = side.qy2 = 0.0;
        m::Load left = side; left.name = "p3 left";
        left.x1 = 0; left.x2 = 0; left.qx1 = left.qx2 = p3;
        pr.loads = {top, side, left};
        m::AnchorMaterial am;
        am.name = "Hold"; am.EA = 1.0e10; am.Lspacing = 1.0;
        pr.anchors.push_back(am);
        m::StructElement hold;
        hold.kind = m::StructKind::Anchor; hold.name = "Corner";
        hold.x1 = 0.0; hold.y1 = 0.0; hold.x2 = -1.0; hold.y2 = 0.0; hold.material = 0;
        pr.structs.push_back(hold);
        pr.initial.load_active = {0, 0, 0};
        pr.initial.struct_active = {1};
        m::Phase load; load.name = "Load"; load.type = m::PhaseType::Plastic; load.load_active = {1, 1, 1};
        load.struct_active = {1}; load.load_steps = 10;
        m::Phase safety; safety.name = "Safety"; safety.type = m::PhaseType::Safety; safety.load_active = {1, 1, 1};
        safety.struct_active = {1};
        pr.phases = {load, safety};
        const auto M = katai::app::mesh_from_project(pr, 0.25, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = M.ok && res.size() == 3 && res[1].ok && res[2].ok && res[2].fos > 0.0;
        if (!ok) {
            std::printf("      (mesh %s; %zu phase(s))\n", M.ok ? "ok" : M.message.c_str(), res.size());
            for (const auto& r : res) std::printf("      (%s)\n", r.message.substr(0, 300).c_str());
        }
        if (ok)
            std::printf("      block under p3 = %.0f, p1 = %.0f kPa, joints at %.0f deg: factor of safety %.4f, Jaeger reduced %.4f (%+.2f%%)\n",
                        p3, p1, dip, res[2].fos, want, 100.0 * (res[2].fos - want) / want);
        check(ok && rel(res[2].fos, want) < 0.01, "(c) the block's factor of safety is the F at which the reduced joint strength carries the state");
    }
    {
        // (d) a rock slope sliding on a daylighting joint set. Cohesionless joints: every plane of
        // the set that daylights in the face cuts a block of the same factor, tan(phi) / tan(dip),
        // so the answer does not hang on the plane through the toe. (With cohesion the smallest
        // factor is the toe plane's, and there the plane runs into the ground in front of the toe:
        // a limit the mesh can only approach from above, so it is no sharp oracle.)
        constexpr double H = 10.0, face = 60.0, gamma = 25.0;
        const double xt = 20.0, yt = 20.0;                       // toe
        const double xc = xt + H / std::tan(rad(face));          // crest
        // Associated slip: K2D-M008 says why (non-associated joints stalled in this slope).
        struct Case { double dip, phi, psi; };
        for (const Case k : {Case{35.0, 45.0, 45.0}, Case{30.0, 40.0, 40.0}}) {
            const double want = std::tan(rad(k.phi)) / std::tan(rad(k.dip));
            m::Project pr;
            m::Material r = rock(k.dip, 0.0, k.phi, gamma, false);
            r.psi = k.psi;
            pr.materials.push_back(r);
            m::SoilPolygon P; P.name = "Rock"; P.material = 0;
            P.x = {0.0, 70.0, 70.0, xc, xt, 0.0};
            P.y = {0.0, 0.0, yt + H, yt + H, yt, yt};
            P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, 0, 0, 0,
                         (int)m::BCType::HorizontallyFixed};
            pr.polygons.push_back(P);
            pr.has_water = false;
            pr.initial_procedure = m::InitialProcedure::GravityLoading;
            m::Phase safety; safety.name = "Safety"; safety.type = m::PhaseType::Safety;
            pr.phases = {safety};
            const auto M = katai::app::mesh_from_project(pr, 1.0, 15);
            const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::GravityLoading);
            const bool ok = M.ok && res.size() == 2 && res[0].ok && res[1].ok && res[1].fos > 0.0;
            if (!ok && !res.empty()) {
                std::printf("      (%s)\n", res.back().message.substr(0, 300).c_str());
                for (const auto& d : res[0].diagnostics) std::printf("      [%s] %s\n", d.code.c_str(), d.message.substr(0, 160).c_str());
            }
            if (ok)
                std::printf("      slope on joints dipping %.0f deg, phi = %.0f, psi = %.0f deg: factor of safety %.4f, tan(phi)/tan(dip) %.4f (%+.2f%%)\n",
                            k.dip, k.phi, k.psi, res[1].fos, want, 100.0 * (res[1].fos - want) / want);
            char lbl[200];
            std::snprintf(lbl, sizeof lbl, "(d) joints at %.0f deg, phi = %.0f, psi = %.0f: the slope slides at tan(phi)/tan(dip), within 2%%", k.dip, k.phi, k.psi);
            check(ok && rel(res[1].fos, want) < 0.02, lbl);
        }
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
