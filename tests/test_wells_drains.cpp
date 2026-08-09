// Wells and drains: the two ways of taking water out of the ground from inside the model.
//
// Until now KATAI could only say what happens at the EDGE of a model -- a head, a seepage face,
// a flux. Dewatering happens in the middle: a well pumps a discharge and the ground draws down
// around it; a drain holds a level and takes whatever reaches it. Neither could be drawn, so an
// excavation being dewatered had to be modelled by moving the phreatic surface by hand, which
// states the answer instead of computing it.
//
// verify: KV-FLW-003
//   oracle:   closed_form
//   source:   Dupuit-Forchheimer one-dimensional confined flow (Darcy's law integrated across a strip of aquifer), the classical basis of the well and trench formulae; PLAXIS 2D 2025.1 Reference Manual sec. 5.9.1-5.9.2 for what the two objects mean -- a well prescribes "the discharge of the well in the unit of volume per unit of time, per unit of width" and stops at h_min ("when the groundwater head reduces below the h_min level no further extraction will occur"), a drain prescribes a head where "the pore pressure in all nodes of the drain is reduced such that it is equivalent to the given head", and a NORMAL drain leaves ground that is already drier alone ("pore pressures lower than the equivalent to the given head are not affected by the drain")
//   locator:  a confined strip of aquifer of thickness b and permeability k, held at head H on both sides, with a FULL-DEPTH vertical line at mid-span: as a well pumping Q, and as a drain holding h_d. The flow to a full-depth line is one-dimensional on each side, so each side carries half of it and Darcy gives the answer exactly
//   quantity: drawdown at the line [m], the discharge the line removes [m3/day per m], and the head profile away from it
//   expected: a well of discharge Q draws the line down to H - Q L / (2 k b) and the head varies linearly to H at each boundary; a drain at h_d removes 2 k b (H - h_d)/L and produces the same linear profiles; a well whose rated discharge would pull the head below h_min instead settles at h_min and removes exactly what a drain at that level would -- the well BECOMES a drain at its floor; and a normal drain set above the surrounding head does nothing at all
//   band:     0.5% on the drawdown and on the extracted discharge, as asserted below (measured -0.00% and +0.00%, the 1D confined field being exactly representable; the head profile away from the line is linear to 5e-14 m); the "does nothing" case is compared bit-for-bit against the SAME MESH with the drain switched off for the phase, since a hydraulic line is a mesh constraint and a model drawn without it is a different discretisation

#include <katai/analysis/seepage.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/flow_driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

// A confined aquifer strip: 40 m long, 8 m thick, held at H on both ends, closed top and bottom.
// The line at mid-span is where the well or the drain goes.
constexpr double kL = 40.0, kB = 8.0, kK = 2.0, kH = 10.0;   // m, m, m/day, m
constexpr double kMid = 0.5 * kL;

m::Project aquifer() {
    m::Project pr;
    pr.name = "KV-FLW-003 dewatering";
    pr.x_min = 0; pr.x_max = kL; pr.y_min = 0; pr.y_max = kB;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, kL};
    pr.wy = {kH, kH};

    m::Material sand;
    sand.name = "Aquifer sand";
    sand.model = m::SoilModel::LinearElastic;
    sand.E = 2.0e4; sand.nu = 0.3;
    sand.gamma_unsat = 18.0; sand.gamma_sat = 20.0;
    sand.kx = kK; sand.ky = kK;
    pr.materials.push_back(sand);

    m::SoilPolygon P;
    P.name = "Aquifer";
    P.material = 0;
    P.x = {0, kL, kL, 0};
    P.y = {0, 0, kB, kB};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    // Confined: head H at both ends, closed above and below. The aquifer is full (H > b), so the
    // flow is saturated everywhere and the linear (confined) solution is the exact one.
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Head};
    P.edge_head = {0.0, kH, 0.0, kH};
    pr.polygons.push_back(P);
    return pr;
}

m::HydroLine well(double Q, double h_min) {
    m::HydroLine H;
    H.name = "Pumping well";
    H.kind = m::HydroKind::Well;
    H.behaviour = (int)m::WellBehaviour::Extraction;
    H.x1 = kMid; H.y1 = 0.0; H.x2 = kMid; H.y2 = kB;   // full depth: the flow to it is 1D
    H.q = Q;
    H.h_min = h_min;
    return H;
}

m::HydroLine drain(double head, m::DrainBehaviour behaviour = m::DrainBehaviour::Normal) {
    m::HydroLine H;
    H.name = "Drain";
    H.kind = m::HydroKind::Drain;
    H.behaviour = (int)behaviour;
    H.x1 = kMid; H.y1 = 0.0; H.x2 = kMid; H.y2 = kB;
    H.head = head;
    return H;
}

// Head at the dewatering line (all its nodes carry the same head in a 1D field).
double head_at_line(const katai::app::FlowResult& R, const katai::mesh::Mesh& mesh) {
    double h = 0.0;
    int count = 0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - kMid) < 1e-9) { h += R.head[n]; ++count; }
    return count ? h / count : 0.0;
}
// Worst deviation from the expected linear profile H -> h_line -> H.
double profile_error(const katai::app::FlowResult& R, const katai::mesh::Mesh& mesh, double h_line) {
    double worst = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double d = std::fabs(mesh.x[n] - kMid);
        const double want = h_line + (kH - h_line) * (d / kMid);
        worst = std::fmax(worst, std::fabs(R.head[n] - want));
    }
    return worst;
}

struct Run {
    katai::app::FlowResult R;
    katai::mesh::Mesh mesh;
    bool ok = false;
};

Run run(const m::Project& pr) {
    Run out;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { out.R.message = M.message; return out; }
    out.mesh = M.mesh;
    out.R = katai::app::solve_groundwater_flow(pr, out.mesh);
    out.ok = out.R.ok;
    return out;
}

}  // namespace

int main() {
    std::printf("== a well pumping from a confined aquifer (KV-FLW-003) ==\n");
    std::printf("  strip L = %.0f m, thickness b = %.0f m, k = %.0f m/day, boundaries at h = %.0f m\n",
                kL, kB, kK, kH);

    // Each half of the strip carries Q/2 over the distance L/2: Darcy gives
    // Q/2 = k b (H - h_well) / (L/2), so the drawdown is Q L / (4 k b) ... per HALF. Writing it
    // once, for the whole strip: h_well = H - Q L / (4 k b) is the drawdown of ONE side feeding
    // half the discharge, which is the same thing.
    // Q is chosen so the drawdown stays above the aquifer roof: these are the CONFINED
    // formulae, and a strip pumped below its own roof is a different (unconfined) problem.
    const double Q = 1.6;                                   // m3/day per m out of plane
    const double s_exact = Q * (kL / 2.0) / (2.0 * kK * kB);  // = (Q/2) (L/2) / (k b)
    {
        m::Project pr = aquifer();
        pr.hydros.push_back(well(Q, -50.0));                // h_min far below: the well is free
        check(katai::io::validate_project(pr).ok(), "the pumped aquifer validates");
        const Run r = run(pr);
        check(r.ok, "the pumped aquifer solves");
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); return 1; }
        const double h_line = head_at_line(r.R, r.mesh);
        const double s = kH - h_line;
        std::printf("  Q = %.4g m3/day/m  ->  drawdown %.6f m (exact %.6f, %+.2f%%), "
                    "profile error %.2e m\n",
                    Q, s, s_exact, 100.0 * (s - s_exact) / s_exact, profile_error(r.R, r.mesh, h_line));
        std::printf("  the run reports the wells and drains removed %.4f m3/day/m\n",
                    r.R.hydro_discharge);
        check(std::fabs(s - s_exact) < 0.005 * s_exact,
              "the drawdown is Q L / (4 k b), the Dupuit answer for a full-depth line");
        check(profile_error(r.R, r.mesh, h_line) < 1e-6 * kH,
              "and the head falls off linearly to the boundaries, as 1D confined flow does");
        check(std::fabs(r.R.hydro_discharge - Q) < 0.005 * Q,
              "the well removes exactly the discharge it was asked for");
        check(r.R.hydro_limited == 0, "and it is not limited: h_min was set far below");
    }

    // A drain at a level takes whatever the ground brings to that level: the same 1D field read
    // the other way round.
    std::printf("\n== a drain holding a level ==\n");
    const double h_d = 9.0;   // above the aquifer roof (8 m): the flow stays confined
    const double Q_drain_exact = 2.0 * kK * kB * (kH - h_d) / (kL / 2.0);
    {
        m::Project pr = aquifer();
        pr.hydros.push_back(drain(h_d));
        check(katai::io::validate_project(pr).ok(), "the drained aquifer validates");
        const Run r = run(pr);
        check(r.ok, "the drained aquifer solves");
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); return 1; }
        const double h_line = head_at_line(r.R, r.mesh);
        std::printf("  drain head %.2f m  ->  head at the line %.6f m, removed %.6f m3/day/m "
                    "(exact %.6f, %+.2f%%)\n",
                    h_d, h_line, r.R.hydro_discharge, Q_drain_exact,
                    100.0 * (r.R.hydro_discharge - Q_drain_exact) / Q_drain_exact);
        check(std::fabs(h_line - h_d) < 1e-9,
              "the drain holds its nodes at exactly the head it was given");
        check(std::fabs(r.R.hydro_discharge - Q_drain_exact) < 0.005 * Q_drain_exact,
              "and removes 2 k b (H - h_d) / (L/2), what the two sides deliver to it");
        check(profile_error(r.R, r.mesh, h_line) < 1e-6 * kH,
              "with the same linear profiles either side");
    }

    // The two objects meet at h_min: a well asked for more than the ground can give at its floor
    // level settles AT that level and removes exactly what a drain there would.
    std::printf("\n== a well that hits its floor becomes a drain at h_min ==\n");
    {
        const double h_min = 9.0;
        m::Project pr = aquifer();
        pr.hydros.push_back(well(100.0, h_min));   // far more than this aquifer can supply
        const Run r = run(pr);
        check(r.ok, "the over-pumped aquifer solves");
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); return 1; }
        const double h_line = head_at_line(r.R, r.mesh);
        std::printf("  rated Q = 100 m3/day/m with h_min = %.2f m  ->  head at the well %.6f m, "
                    "removed %.6f m3/day/m (a drain at that level: %.6f)\n",
                    h_min, h_line, r.R.hydro_discharge, Q_drain_exact);
        check(std::fabs(h_line - h_min) < 1e-9,
              "the head at the well stops exactly at h_min, as the manual says it must");
        check(std::fabs(r.R.hydro_discharge - Q_drain_exact) < 0.005 * Q_drain_exact,
              "and the extraction is what the ground can supply at that head, not the rated 100");
        check(r.R.hydro_limited > 0, "the run reports that the well is limited by h_min");
        check(r.R.message.find("h_min") != std::string::npos,
              "and says so in words, so a reader does not take the rated discharge for the answer");
    }

    // A NORMAL drain above the water it sits in does nothing at all -- not approximately, but
    // exactly: the same model without it gives the same field.
    std::printf("\n== a normal drain above the ground water does nothing ==\n");
    {
        // The comparison is the SAME model with the drain switched OFF for the phase, not a model
        // drawn without it: a hydraulic line is a mesh constraint, so leaving it out changes the
        // discretisation, and the two runs would then differ for a reason that has nothing to do
        // with drains. Per-phase activity is exactly the switch this needs.
        m::Project pr = aquifer();
        pr.hydros.push_back(drain(kH + 5.0));       // five metres above the boundary head
        const auto M = katai::app::mesh_from_project(pr);
        check(M.ok, "the model with a high drain meshes");
        m::Phase off, on;
        off.hydro_active = {0};
        on.hydro_active = {1};
        const auto Roff = katai::app::solve_groundwater_flow(pr, M.mesh, &off);
        const auto Ron = katai::app::solve_groundwater_flow(pr, M.mesh, &on);
        check(Roff.ok && Ron.ok, "both runs solve");
        if (Roff.ok && Ron.ok) {
            double worst = 0.0;
            for (int n = 0; n < M.mesh.node_count; ++n)
                worst = std::fmax(worst, std::fabs(Roff.head[n] - Ron.head[n]));
            std::printf("  worst head difference against the same mesh with the drain switched "
                        "off: %.1e m; removed %.3g m3/day/m\n", worst, Ron.hydro_discharge);
            check(worst == 0.0,
                  "a normal drain above the surrounding head leaves the field bit-identical");
            check(Ron.hydro_discharge == 0.0, "and removes nothing");
        }
        // A VACUUM drain at the same level does NOT leave it alone -- that is the difference
        // between the two behaviours, and it is what vacuum consolidation is.
        m::Project pv = aquifer();
        pv.hydros.push_back(drain(kH + 5.0, m::DrainBehaviour::Vacuum));
        const Run vac = run(pv);
        check(vac.ok, "the vacuum-drain model solves");
        if (vac.ok) {
            const double h_line = head_at_line(vac.R, vac.mesh);
            std::printf("  the same level as a VACUUM drain holds the line at %.3f m "
                        "(the ground is pushed up to it)\n", h_line);
            check(std::fabs(h_line - (kH + 5.0)) < 1e-9,
                  "a vacuum drain holds its head in both directions");
        }
    }

    // The refusals: a well that pumps nothing, and a line the mesh cannot see.
    std::printf("\n== what is refused ==\n");
    {
        m::Project pr = aquifer();
        pr.hydros.push_back(well(0.0, 0.0));
        const auto rep = katai::io::validate_project(pr);
        check(!rep.ok(), "a well with no discharge is refused by the input contract");
    }
    {
        m::Project pr = aquifer();
        m::HydroLine H = well(10.0, -50.0);
        H.x1 = H.x2 = kL + 20.0;             // outside the model
        pr.hydros.push_back(H);
        const Run r = run(pr);
        check(!r.ok && r.R.message.find("does not lie on the mesh") != std::string::npos,
              "a well drawn outside the soil is refused rather than pumping nothing: " +
                  r.R.message.substr(0, 60) + "...");
    }

    std::printf(g_failures ? "\nFAILED (%d)\n" : "\nAll checks passed.\n", g_failures);
    return g_failures ? 1 : 0;
}
