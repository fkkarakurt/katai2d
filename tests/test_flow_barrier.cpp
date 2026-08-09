// A wall that the water cannot walk through.
//
// A cut-off wall is drawn to hold a head difference. Until now this program drew it, meshed it,
// and then computed the flow straight through it: the seepage under an excavation was the seepage
// of a model with no wall in it, and the uplift and the inflow that followed belonged to a
// different problem. The deformation half of the program has always split its mesh along a wall
// (that is how soil-structure interaction works here); the flow half never did, so the two sides
// of the wall shared one pore-pressure degree of freedom and were therefore the same water.
//
// verify: KV-FLW-004
//   oracle:   closed_form
//   source:   PLAXIS 2D 2025.1 Reference Manual Table 5-2 and sec. 6.1.7.4, Scientific Manual sec. 3.4 -- an interface is Impermeable ("a full separation of the pore pressure degrees-of-freedom of the interface node pairs", zero cross permeability), Semi-permeable with a HYDRAULIC RESISTANCE d/k in units of time ("to determine d/k, one needs to measure the average discharge q through a wall (per unit of area) for a given head difference dh, so d/k = dh/q"), or Fully permeable. Read against Darcy's law in series: the two halves of a confined aquifer and the screen between them are three resistances the same discharge passes through
//   locator:  a confined aquifer strip of length L, thickness b and permeability k, held at two different heads at its ends, with a full cut-off at mid-span: absent, impermeable, and semi-permeable at two resistances
//   quantity: the discharge through the strip [m3/day per m] and the head on each side of the barrier [m]
//   expected: with no barrier, q = k b dh / L; with an impermeable barrier, q = 0 EXACTLY and each side stands at its own boundary head, so the barrier holds the whole difference; with a semi-permeable barrier, q = dh / (L/(k b) + R/b) -- the soil either side and the screen in series -- which for the values used here is exactly half of the unobstructed discharge at R = 20 days
//   band:     0.5% on the discharges, as asserted below (measured +0.00% unobstructed, -0.10% and -0.15% semi-permeable -- the seam's leakage is integrated with linear segments along a quadratic element edge); the impermeable case is asserted at 1e-9 of the unobstructed discharge (measured 3e-13 relative) and its two side heads against the boundary heads to 1e-12 m

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

// Both boundary heads stay ABOVE the aquifer roof, so the strip is confined everywhere and the
// series-resistance answer is the exact one.
constexpr double kL = 40.0, kB = 8.0, kK = 2.0;
constexpr double kH1 = 12.0, kH2 = 10.0, kDh = kH1 - kH2;
constexpr double kMid = 0.5 * kL;

// barrier: 0 = none, otherwise the cross-permeability setting (1 impermeable, 2 semi-permeable).
m::Project strip(int barrier, double resistance) {
    m::Project pr;
    pr.name = "KV-FLW-004 cut-off";
    pr.x_min = 0; pr.x_max = kL; pr.y_min = 0; pr.y_max = kB;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, kL};
    pr.wy = {kH1, kH1};

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
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Head};
    P.edge_head = {0.0, kH2, 0.0, kH1};
    pr.polygons.push_back(P);

    if (barrier) {
        // The line runs from BELOW the base to ABOVE the surface on purpose. The splitter keeps
        // the end nodes shared -- which is the manual's own rule, "the end points of an interface
        // are always permeable" -- so a screen that stops inside the soil leaks around its ends.
        // Here the ends are outside the mesh, which is what a full cut-off means.
        m::StructElement w;
        w.kind = m::StructKind::Interface;
        w.name = "Cut-off";
        w.x1 = kMid; w.y1 = -1.0; w.x2 = kMid; w.y2 = kB + 1.0;
        w.material = -1;
        w.flow_barrier = barrier;
        w.hydraulic_resistance = resistance;
        pr.structs.push_back(w);
    }
    return pr;
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

// Mean head on the two sides of the barrier line (the near side and, through head_far, the far one).
void heads_at_barrier(const Run& r, double& near_side, double& far_side) {
    double a = 0.0, b = 0.0;
    int count = 0;
    for (int n = 0; n < r.mesh.node_count; ++n) {
        if (std::fabs(r.mesh.x[n] - kMid) > 1e-9) continue;
        a += r.R.head[n];
        b += r.R.head_far.empty() ? r.R.head[n] : r.R.head_far[n];
        ++count;
    }
    near_side = count ? a / count : 0.0;
    far_side = count ? b / count : 0.0;
}

}  // namespace

int main() {
    std::printf("== a confined strip with and without a cut-off (KV-FLW-004) ==\n");
    std::printf("  L = %.0f m, b = %.0f m, k = %.0f m/day, heads %.0f -> %.0f m\n",
                kL, kB, kK, kH1, kH2);

    const double q_open = kK * kB * kDh / kL;          // Darcy, no obstruction
    {
        const Run r = run(strip(0, 0.0));
        check(r.ok, "the unobstructed strip solves");
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); return 1; }
        std::printf("  no barrier            : q = %.6f m3/day/m (exact %.6f, %+.2f%%)\n",
                    r.R.discharge, q_open, 100.0 * (r.R.discharge - q_open) / q_open);
        check(std::fabs(r.R.discharge - q_open) < 0.005 * q_open,
              "the strip carries k b dh / L, as Darcy says");
    }

    // An impermeable screen: no path at all, so no water -- and the whole head difference stands
    // across the line, which is the number a cut-off is designed for.
    {
        const m::Project pr = strip(1, 0.0);
        check(katai::io::validate_project(pr).ok(), "the cut-off model validates");
        const Run r = run(pr);
        check(r.ok, "the cut-off model solves");
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); return 1; }
        double near_side = 0.0, far_side = 0.0;
        heads_at_barrier(r, near_side, far_side);
        std::printf("  impermeable cut-off   : q = %.3e m3/day/m;  heads across the line "
                    "%.9f | %.9f m (boundary heads %.0f | %.0f)\n",
                    r.R.discharge, near_side, far_side, kH1, kH2);
        // Not "== 0": the two sides are exactly decoupled, but the reported discharge is a sum of
        // K h products and carries the solver's round-off. What is asserted is that it is nothing
        // NEXT TO the flow the same strip carries without the wall.
        check(std::fabs(r.R.discharge) < 1e-9 * q_open,
              "no water crosses an impermeable barrier (round-off next to the open strip)");
        const double hi = std::fmax(near_side, far_side), lo = std::fmin(near_side, far_side);
        check(std::fabs(hi - kH1) < 1e-12 && std::fabs(lo - kH2) < 1e-12,
              "and each side stands at its own boundary head: the barrier holds the whole "
              "difference");
        check(std::fabs((hi - lo) - kDh) < 1e-12,
              "which is the head difference the wall was drawn to hold");
    }

    // A semi-permeable screen: the soil either side and the screen are three resistances in
    // series, and the manual defines the middle one as R = d/k = dh / q per unit area.
    for (const double R_res : {20.0, 5.0}) {
        const double q_ex = kDh / (kL / (kK * kB) + R_res / kB);
        const Run r = run(strip(2, R_res));
        check(r.ok, "the semi-permeable model solves at R = " + std::to_string((int)R_res));
        if (!r.ok) { std::printf("      (%s)\n", r.R.message.c_str()); continue; }
        double near_side = 0.0, far_side = 0.0;
        heads_at_barrier(r, near_side, far_side);
        std::printf("  semi-permeable R = %-4.0f: q = %.6f m3/day/m (exact %.6f, %+.2f%%);  "
                    "the screen itself holds %.4f m of the %.0f m\n",
                    R_res, r.R.discharge, q_ex, 100.0 * (r.R.discharge - q_ex) / q_ex,
                    std::fabs(near_side - far_side), kDh);
        check(std::fabs(r.R.discharge - q_ex) < 0.005 * q_ex,
              "q = dh / (L/(k b) + R/b): the soil and the screen in series");
        // The screen's own share of the head difference is q R / b -- the same statement read
        // through the wall instead of through the soil.
        const double dh_screen = r.R.discharge * R_res / kB;
        check(std::fabs(std::fabs(near_side - far_side) - dh_screen) < 0.01 * dh_screen,
              "and the drop across the screen itself is q R / b");
    }

    // The limits, which is where a resistance model is easiest to get wrong: a vanishing
    // resistance must reproduce the unobstructed strip, a huge one the cut-off.
    {
        const Run tiny = run(strip(2, 1e-4));
        const Run huge = run(strip(2, 1e8));
        check(tiny.ok && huge.ok, "both limiting resistances solve");
        if (tiny.ok && huge.ok) {
            std::printf("  R -> 0 : q = %.6f (unobstructed %.6f)    R -> inf : q = %.3e "
                        "(cut-off 0)\n", tiny.R.discharge, q_open, huge.R.discharge);
            check(std::fabs(tiny.R.discharge - q_open) < 0.01 * q_open,
                  "a vanishing resistance is a wall that is not there");
            check(huge.R.discharge < 1e-6 * q_open,
                  "and an enormous one is a cut-off");
        }
    }

    std::printf(g_failures ? "\nFAILED (%d)\n" : "\nAll checks passed.\n", g_failures);
    return g_failures ? 1 : 0;
}
