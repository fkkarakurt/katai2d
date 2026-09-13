// Hoek-Brown through the calculation path -- the step that makes the model REACHABLE. The core
// and its closed forms are KV-CST-014; what is asked here is whether a `.k2d` file that names the
// model produces the rock the published criterion describes, which is a different question and
// the one a user of the program depends on.
//
// THE EXPERIMENT IS A TRIAXIAL ELEMENT TEST, run twice at two cell pressures, because that is
// where a rock model differs from a soil model and where a Mohr-Coulomb fit cannot follow. A block
// of rock is brought to an isotropic cell pressure and then squeezed axially past yield; the axial
// stress it settles on must be the criterion's own envelope
//     sigma_1 = sigma_3 - |sigma_ci| ( m_b (-sigma_3 / |sigma_ci|) + s )^a
// evaluated at the lateral stress the specimen is actually carrying. Doing it at TWO cell
// pressures is the point: the envelope is CURVED, so the ratio of the two failure stresses is a
// number no straight line can also produce, and a model that had quietly fallen back to
// Mohr-Coulomb would miss the second one even if it were tuned to hit the first.
//
// IT IS DRIVEN BY DISPLACEMENT, and that is a finding rather than a preference. The first working
// version pushed the specimen with a LOAD and read the highest level the phase could equilibrate.
// It reproduced the envelope to -0.99% at the high cell pressure and only -5.72% at the low one,
// and the deficit was not the load step: cutting the increment from 395 to 182 kPa moved the error
// from -5.65% to -5.72%, which is to say not at all. The cause is the experiment. This specimen is
// HOMOGENEOUS and the model is PERFECTLY PLASTIC, so every point reaches the envelope at the same
// instant and there is no reserve anywhere to redistribute into; a load-controlled Newton can only
// approach that limit from below and stops at the last increment it happened to equilibrate. Under
// displacement control the same state is an equilibrium the solver can stand on, and the stress is
// then READ rather than inferred -- which also removes the assumption the old fixture had to make
// about which principal stress the confinement was.
//
// verify: KV-CST-015
//   oracle:   closed_form
//   source:   Hoek, E., Carranza-Torres, C. & Corkum, B. (2002). Hoek-Brown failure criterion -- 2002 edition. Proc. NARMS-TAC Conference, Toronto, 267-273 -- the generalised criterion with its rock-mass constants m_b, s and a; the same criterion KV-CST-014 checks at the material point, asked here of the assembled FE path from a project file (KATAI 2D input contract, docs/k2d-format.md, sigci / mi / gsi / hbD)
//   locator:  sigma_1 = sigma_3 - |sigma_ci| (m_b (-sigma_3/|sigma_ci|) + s)^a, compression negative, with m_b = m_i exp((GSI-100)/(28-14D)), s = exp((GSI-100)/(9-3D)), a = 1/2 + (exp(-GSI/15) - exp(-20/3))/6
//   quantity: the axial effective stress a displacement-driven triaxial specimen carries on the plastic plateau [kPa], at two cell pressures, and the ratio between the two
//   expected: the Hoek-Brown envelope evaluated at the lateral stress the specimen is carrying, and a ratio that is the curve's, not a line's
//   band:     0.1% on each plateau and on their ratio -- MEASURED: -0.0000% at both cell pressures (sigma_1 = -10490.7 against the envelope's -10490.7, and -33891.4 against -33891.4), with the block 0.05% / 0.02% from homogeneous and the phase carrying the full imposed settlement (lambda = 1.000). The band is the margin over that measurement, not a tolerance the answer needs; the load-controlled version of this same fixture needed 3% and used all of it.
#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
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

constexpr double kW = 1.0, kH = 1.0;      // the specimen [m] -- see the note in main()
constexpr double kSigCi = 50000.0;        // intact rock, 50 MPa
constexpr double kMi = 10.0, kGsi = 50.0, kD = 0.0;

// The published envelope, written here rather than taken from the header under test.
double envelope_sigma1(double sigma3) {
    const double mb = kMi * std::exp((kGsi - 100.0) / (28.0 - 14.0 * kD));
    const double s = std::exp((kGsi - 100.0) / (9.0 - 3.0 * kD));
    const double a = 0.5 + (std::exp(-kGsi / 15.0) - std::exp(-20.0 / 3.0)) / 6.0;
    return sigma3 - kSigCi * std::pow(mb * (-sigma3 / kSigCi) + s, a);
}

// A rock specimen under a cell pressure `conf` (a compression, so a positive number here is the
// magnitude) on its two free faces, then squeezed from the top by an imposed settlement `squeeze`.
// A quarter of a smooth-platen specimen: the base is a roller and the left edge the symmetry line,
// so the block is free to expand sideways and the stress state stays homogeneous.
m::Project specimen(double conf, double squeeze) {
    m::Project pr;
    m::Material r; r.name = "Sandstone"; r.model = m::SoilModel::HoekBrown;
    r.E = 5.0e6; r.nu = 0.25;
    r.gamma_unsat = 25.0; r.gamma_sat = 25.0;   // a real rock weight; see the note in main()
    r.sig_ci = kSigCi; r.mi = kMi; r.gsi = kGsi; r.hb_D = kD;
    r.psi = 0.0; r.sig_psi = 0.0;
    pr.materials.push_back(r);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    // bottom roller, right free (the cell pressure acts there), top free (cell pressure, then the
    // imposed settlement), left on the symmetry line.
    P.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed};
    P.edge_head = {0, 0, 0, 0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::Load side; side.kind = m::LoadKind::Distributed; side.name = "Cell pressure (lateral)";
    side.x1 = kW; side.y1 = 0.0; side.x2 = kW; side.y2 = kH;
    side.qx1 = side.qx2 = -conf; side.qy1 = side.qy2 = 0.0;   // pushing inwards (-x)
    pr.loads.push_back(side);

    // THE CELL PRESSURE ACTS ON THE TOP TOO. An earlier version of this fixture squeezed the
    // specimen laterally with the top free -- which is not a confinement at all but a deviatoric
    // state of the confining magnitude, and at 8 MPa it is five times what this rock mass can
    // carry unconfined (its sigma_c is 3.0 MPa). The specimen failed in the INITIAL phase and the
    // run said so; the mistake was the experiment's, not the model's. A triaxial test starts
    // isotropic, and so does this one.
    m::Load cell; cell.kind = m::LoadKind::Distributed; cell.name = "Cell pressure (axial)";
    cell.x1 = 0.0; cell.y1 = kH; cell.x2 = kW; cell.y2 = kH;
    cell.qx1 = cell.qx2 = 0.0; cell.qy1 = cell.qy2 = -conf;
    pr.loads.push_back(cell);

    // The platen. In the loading phase the axial cell pressure is switched off and this takes its
    // place: the top face is driven down by a fixed amount, and whatever stress that needs is the
    // answer. The lateral cell pressure stays on throughout -- it is the sigma_3 of the test.
    m::PrescribedDisp plat; plat.name = "Platen";
    plat.x1 = 0.0; plat.y1 = kH; plat.x2 = kW; plat.y2 = kH;
    plat.set_ux = false; plat.set_uy = true; plat.uy = -squeeze;
    pr.disps.push_back(plat);

    pr.initial.load_active = {1, 1};      // isotropic cell pressure
    pr.initial.disp_active = {0};         // the platen is not down yet
    m::Phase ph; ph.name = "Axial squeeze"; ph.type = m::PhaseType::Plastic;
    ph.load_active = {1, 0};              // lateral pressure stays; the axial one hands over
    ph.disp_active = {1};
    // THE COST OF AN ELASTIC TANGENT, PAID HERE. This model hands the solver the elastic operator
    // rather than a consistent one (hoek_brown.hpp says why), so on the plastic plateau Newton
    // converges linearly and needs the iterations.
    ph.load_steps = 40;
    ph.max_iterations = 600;
    pr.phases.push_back(ph);
    return pr;
}

// Run the element test and read the stress in the middle of the specimen. Returns false only if
// the run itself did not happen; `spread` reports how far from homogeneous the block ended up,
// which is what says the reading is a MATERIAL answer and not a boundary artefact.
bool element_test(double conf, double squeeze, double& sxx, double& syy, double& spread,
                  double& lam) {
    // THROUGH THE FILE, not past it. The headline of this test is that a `.k2d` naming
    // Hoek-Brown gets the published criterion's rock, and a programmatic project would only have
    // shown that the ENGINE does -- the five rock keys are written only by a material that uses
    // the model, so a writer or reader that dropped one would leave a default-parameter rock
    // behind and every number below would still be a rock's. So the project is serialised,
    // parsed back, and it is the PARSED one that is meshed and solved.
    const auto built = specimen(conf, squeeze);
    const std::string text = m::project_to_json(built);
    m::Project pr;
    std::string err;
    std::vector<katai::io::Issue> notes;
    if (!m::project_from_json(text, pr, &err, &notes)) {
        std::printf("   DIAG .k2d reader: %s\n", err.c_str());
        return false;
    }
    for (const auto& n : notes) std::printf("   DIAG reader note: %s\n", n.message.c_str());
    const auto issues = katai::io::validate_project(pr);
    for (const auto& is : issues.issues)
        if (is.severity == katai::io::Severity::Error) {
            std::printf("   DIAG validator: %s: %s\n", is.path.c_str(), is.message.c_str());
            return false;
        }
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    if (!M.ok) { std::printf("   DIAG mesh: %s\n", M.message.c_str()); return false; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2) {
        std::printf("   DIAG phases=%zu%s\n", res.size(),
                    res.empty() ? "" : (" msg: " + res[0].message).c_str());
        return false;
    }
    if (!res[0].ok) std::printf("   DIAG initial: %s\n", res[0].message.c_str());
    if (!res[1].ok) std::printf("   DIAG phase: %s\n", res[1].message.c_str());
    lam = res[1].load_factor;
    const auto& mesh = res[1].mesh;
    const auto& S = res[1].stress.stress;
    // Every node strictly inside the block: the platen and the roller are the only places where a
    // smooth-platen idealisation can be violated, so the interior is where the material speaks.
    double sxx_sum = 0.0, syy_sum = 0.0, syy_min = 1e300, syy_max = -1e300;
    int n = 0;
    for (int i = 0; i < mesh.node_count && i < (int)S.size(); ++i) {
        if (mesh.y[i] < 0.15 * kH || mesh.y[i] > 0.85 * kH) continue;
        sxx_sum += S[i][0]; syy_sum += S[i][1];
        if (S[i][1] < syy_min) syy_min = S[i][1];
        if (S[i][1] > syy_max) syy_max = S[i][1];
        ++n;
    }
    if (n == 0) { std::printf("   DIAG no interior nodes\n"); return false; }
    sxx = sxx_sum / n; syy = syy_sum / n;
    spread = (syy_max - syy_min) / std::fabs(syy);
    return true;
}

}  // namespace

int main() {
    // A NOTE ABOUT THE SPECIMEN, because the first version of it did not work and the reason is
    // not this model's. It was weightless -- the cleanest way to make the stress state the load's
    // alone -- and the initial phase would not converge: 0% of the applied load equilibrated after
    // 200 iterations per increment. The same fixture with a very strong MOHR-COULOMB material
    // failed identically, and giving either of them a real unit weight fixed both (lambda = 1.000
    // at gamma = 25). So it is the K0 initial phase that does not take a weightless material, not
    // the rock model, and this test is not the place to chase it. The specimen is therefore heavy
    // and SMALL: at 1 m tall the self-weight reaches 25 kPa, which is 0.24% of the smallest
    // failure stress measured below and well inside the band.
    std::printf("Hoek-Brown through the calculation path\n\n");
    std::printf("   cell press.   sigma_1 plateau   envelope at sigma_3      error   sigma_3 read"
                "   spread   lambda\n");
    double reached[2] = {0.0, 0.0}, envel[2] = {0.0, 0.0};
    const double confs[2] = {1000.0, 8000.0};
    for (int i = 0; i < 2; ++i) {
        // Squeeze well past yield: the elastic strain the deviator needs is (sigma_1-sigma_3)/E,
        // which is 0.19% at the low cell pressure and 0.52% at the high one, so 2% of the
        // specimen's height leaves the block sitting on the plateau in both runs.
        const double squeeze = 0.02 * kH;
        double sxx = 0.0, syy = 0.0, spread = 0.0, lam = 0.0;
        const bool ran = element_test(confs[i], squeeze, sxx, syy, spread, lam);
        check(ran, "the rock specimen ran");
        if (!ran) continue;   // silent-drop-ok: the failed check above is the report
        // The envelope is evaluated at the lateral stress the specimen IS CARRYING, not at the
        // pressure applied to its face: that way the test states what the FE stress point does,
        // and does not depend on the fixture reproducing sigma_3 exactly.
        reached[i] = syy; envel[i] = envelope_sigma1(sxx);
        std::printf("   %10.1f   %15.1f   %17.1f   %+8.4f%%   %12.1f   %6.2f%%   %6.3f\n",
                    -confs[i], reached[i], envel[i], 100.0 * (reached[i] / envel[i] - 1.0), sxx,
                    100.0 * spread, lam);
        check(spread < 0.005, "and it stayed homogeneous, so the reading is the material's");
        check(std::fabs(reached[i] / envel[i] - 1.0) < 0.001,
              "the stress it carries is the Hoek-Brown envelope at that confinement, within 0.1%");
    }
    // ---- AND THE ANALYSIS THIS MODEL CANNOT BE ASKED FOR --------------------------------
    // phi-c reduction divides c and phi; Hoek-Brown reads neither, so every trial would run the
    // rock at full strength, nothing would ever bring it down, and the search would report its
    // cap itself -- "FoS > 3.0", whatever the rock. The refusal is what is checked here,
    // because a wrong number in the safe-looking direction is the worst kind this program can
    // produce. (A strength reduction that is defined for this model has to reformulate the
    // Hoek-Brown yield function itself, and this build does not have it.)
    {
        auto pr = specimen(1000.0, 0.001);
        m::Phase sf; sf.name = "FoS"; sf.type = m::PhaseType::Safety;
        sf.load_active = {1, 0}; sf.disp_active = {0};
        pr.phases.push_back(sf);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool refused = res.size() == 3 && !res[2].ok &&
                             res[2].message.find("Hoek-Brown") != std::string::npos;
        check(refused, "a Safety phase on rock is REFUSED, not answered with an inflated factor");
        if (res.size() == 3) std::printf("      %s\n", res[2].message.c_str());
        check(res.size() == 3 && res[2].fos < 0.0,
              "and no factor of safety is reported at all");
        // The same refusal at the schema, before a mesh is ever built.
        const auto rep = katai::io::validate_project(pr);
        bool schema = false;
        for (const auto& is : rep.issues)
            schema |= is.severity == katai::io::Severity::Error &&
                      is.message.find("Hoek-Brown") != std::string::npos;
        check(schema, "and the .k2d contract refuses it too, before a mesh is built");
    }

    // ---- AND THE DESIGN APPROACH THAT WOULD HAVE FACTORED NOTHING ------------------------
    // Same defect at a different seam, and worse: EC7 DA1-C2 / DA3 divide c' and tan(phi'), so
    // the rock would be solved at its CHARACTERISTIC strength under a report saying the design
    // approach had been applied. A design verification that quietly used unfactored strength is
    // not a conservative approximation, it is a wrong verdict.
    {
        auto pr = specimen(1000.0, 0.001);
        pr.phases[0].design_approach = m::DesignApproach::EC7_DA3;
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool refused = false;
        for (const auto& x : res)
            refused |= !x.ok && x.message.find("Hoek-Brown material") != std::string::npos;
        check(refused, "a material-factored design approach on rock is REFUSED");
        const auto rep = katai::io::validate_project(pr);
        bool schema = false;
        for (const auto& is : rep.issues)
            schema |= is.severity == katai::io::Severity::Error &&
                      is.message.find("no partial factor") != std::string::npos;
        check(schema, "and the .k2d contract names EN 1997-1 in refusing it");
        // The resistance-factored approaches never touch the material, so they still run.
        auto ok_pr = specimen(1000.0, 0.001);
        ok_pr.phases[0].design_approach = m::DesignApproach::EC7_DA2;
        const auto M2 = katai::app::mesh_from_project(ok_pr, 0.5, 6);
        const auto res2 = katai::app::solve_phases(ok_pr, M2.mesh, InitialPhase::K0Procedure);
        check(res2.size() == 2 && res2[1].ok,
              "while a resistance-factored one runs, because it factors no material");
    }

    // ---- AND THE JOINT THAT WOULD HAVE BORROWED A STRENGTH THE ROCK HAS NOT GOT ----------
    // An interface takes c_i = R * c' and phi_i from the material beside it. Beside rock those
    // boxes hold the schema defaults -- 1 kPa and 30 degrees -- which the Hoek-Brown model never
    // reads, so the joint would have been given a Mohr-Coulomb strength unrelated to the rock it
    // is cut into. The remedy the refusal names already exists in the schema, and the second half
    // of this check is that it WORKS: point the interface at a Mohr-Coulomb material and the same
    // model runs.
    {
        auto pr = specimen(1000.0, 0.001);
        m::PlateMaterial pm; pm.name = "Liner"; pm.EA = 6.0e6; pm.EI = 4.0e4; pm.w = 1.0;
        pr.plates.push_back(pm);
        m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Liner";
        w.x1 = 0.5 * kW; w.y1 = 0.0; w.x2 = 0.5 * kW; w.y2 = kH;
        w.material = 0; w.iface_pos = true; w.iface_neg = true;
        pr.structs.push_back(w);
        // Active from the start: a plate that carries interfaces cannot be switched on per phase
        // in this build, and the run says so -- which is not what this case is about.
        pr.initial.struct_active = {1};
        pr.phases[0].struct_active = {1};
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool refused = false;
        for (const auto& x : res)
            refused |= !x.ok && x.message.find("Hoek-Brown model -- and that model has no c'") !=
                                    std::string::npos;
        for (const auto& x : res)
            if (!x.ok) std::printf("      %s\n", x.message.c_str());
        check(refused, "an interface against rock is REFUSED, not given the unused c'/phi' boxes");
        // The remedy: a Mohr-Coulomb material named as the interface's own, which is how a rock
        // JOINT is modelled -- the discontinuity's friction, not the rock mass's envelope.
        m::Material j; j.name = "Joint"; j.model = m::SoilModel::MohrCoulomb;
        j.E = 5.0e6; j.nu = 0.25; j.c = 30.0; j.phi = 35.0;
        j.gamma_unsat = j.gamma_sat = 25.0;
        pr.materials.push_back(j);
        pr.structs[0].iface_material = 1;
        const auto M2 = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res2 = katai::app::solve_phases(pr, M2.mesh, InitialPhase::K0Procedure);
        check(res2.size() == 2 && res2[1].ok,
              "and the remedy the refusal names actually runs the same model");
    }

    // THE CURVE, not a line. A Mohr-Coulomb fit through the first point would predict the second
    // by a straight line; the criterion's own ratio is smaller, and that difference is the whole
    // reason this model exists.
    const double ratio = envel[1] / envel[0];
    const double linear = (confs[1] / confs[0]);
    std::printf("   the envelope's ratio between the two confinements is %.3f; a line through the "
                "origin would give %.3f\n", ratio, linear);
    check(std::fabs(reached[1] / reached[0] - ratio) < 0.001,
          "and the two runs reproduce the CURVE's ratio, which a straight line cannot");

    if (g_failures == 0) {
        std::printf("\nOK: a .k2d file that names Hoek-Brown gets the published rock criterion\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
