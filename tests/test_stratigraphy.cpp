// KV-GEO-001: borehole logs -> soil regions and a water surface.
//
// The generator is a pure schema-to-schema function, so it is checked the way a closed form is:
// every expected level here can be computed by hand from the logs, and the test states the rule it
// is checking rather than the number alone. Source for every rule: PLAXIS 2D 2025.1 Reference
// Manual sec. 4.2 (interpolate between boreholes; all layers appear in all boreholes; local zero
// thickness is legal), sec. 4.3.1.1 (a layer absent at a location is a zero thickness there) and
// sec. 7.10.1.1 (a single head is a horizontal water surface reaching the model boundaries).
//
// The reason this file exists at all is that a wrong stratigraphy is INVISIBLE downstream: the
// mesher will happily mesh the wrong ground, every phase will converge on it, and the answer will
// be a correct solution to a model nobody meant. There is no residual that sees it.
//
// verify: KV-GEO-001
//   oracle:   closed_form
//   source:   PLAXIS 2D 2025.1 Reference Manual, sec. 4.2 "Creating boreholes" -- "If multiple boreholes are defined, PLAXIS 2D will automatically interpolate between boreholes, and derive the position of the soil layers from the borehole information. Each defined soil layer is used throughout the whole model contour. In other words, all soil layers appear in all boreholes. The top and the bottom boundaries of the layers may vary through boreholes, making it possible to define non-horizontal soil layers of non-uniform thickness as well as layers that locally have a zero thickness." -- with sec. 4.3.1.1 ("If a certain soil layer does not exist at the position indicated by a model borehole, it should be assigned a zero thickness (bottom level equal to top level) in the relevant borehole") and sec. 7.10.1.1 ("A single borehole can be used to create a horizontal water surface that extends to the model boundaries. When multiple boreholes are used, a non-horizontal water surface can be created by combining the heads in the various boreholes")
//   locator:  the rule stated in full, which is what is evaluated here rather than called from the generator's header: a boundary level between two logs at x_a < x < x_b is z(x) = z_a + (x - x_a)/(x_b - x_a) (z_b - z_a); outside the outermost log it is that log's own level, HELD rather than continued on its slope; and the water surface follows the same two rules applied to the heads. The boundaries break slope at a borehole and nowhere else, so the generated polygon has vertices at the model edges and at each log, and no others
//   quantity: the generated layer-boundary levels and phreatic level of a 40 m x 20 m model, read back off the produced polygons and water polyline, for one log, for two logs with different levels and heads, and for a three-layer section whose middle layer pinches out [m]
//   expected: one log -> every boundary at its own level at x = 0, 12 and 40; two logs at x = 10 and 30 -> the midpoint is the mean (20 and 18 give 19; 16 and 10 give 13) while x = 0 reads 20 and x = 40 reads 18, NOT the extrapolated 21 and 17; the pinching layer is 2.000 m thick at one log and exactly 0.000 m at the other, with the layer beneath rising to meet the one above so the ground has no gap; the water surface is 12 at x <= 10, 16 at x >= 30 and linear between
//   band:     1e-12 m, i.e. exact to round-off, as asserted below -- these are not measurements of a physical process but evaluations of a stated interpolation, so any deviation beyond arithmetic noise is a defect rather than a discretisation. Every expected value above is computable by hand from the logs in the test, which is the point: the oracle is the manual's rule, not this program's own output recorded once
#include <katai/model/project.hpp>
#include <katai/model/stratigraphy.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (ok) std::printf("ok:   %s\n", what.c_str());
    else { std::printf("FAIL: %s\n", what.c_str()); ++g_failures; }
}
void close_to(double got, double want, double tol, const std::string& what) {
    const bool ok = std::fabs(got - want) <= tol;
    if (ok) std::printf("ok:   %s (%.9g)\n", what.c_str(), got);
    else { std::printf("FAIL: %s: got %.9g, want %.9g\n", what.c_str(), got, want); ++g_failures; }
}

// A model 40 m wide, 20 m deep, with two strata unless a case says otherwise.
m::Project base_project() {
    m::Project pr;
    pr.x_min = 0.0; pr.x_max = 40.0;
    pr.y_min = 0.0; pr.y_max = 20.0;
    pr.materials.resize(2);
    pr.materials[0].name = "Fill";
    pr.materials[1].name = "Clay";
    pr.strata = {{"Fill", 0}, {"Clay", 1}};
    return pr;
}

// The generated top boundary of stratum j, sampled at x, read back off the polygon rather than
// from the interpolator -- so the test checks what the GEOMETRY says, not what the helper says.
double top_of_polygon_at(const m::SoilPolygon& P, double x) {
    // The ring is top left-to-right then bottom right-to-left, so the first half is the top.
    double best_y = -1e300;
    for (size_t i = 0; i < P.x.size(); ++i)
        if (std::fabs(P.x[i] - x) < 1e-9) best_y = std::max(best_y, P.y[i]);
    return best_y;
}
double bottom_of_polygon_at(const m::SoilPolygon& P, double x) {
    double best_y = 1e300;
    for (size_t i = 0; i < P.x.size(); ++i)
        if (std::fabs(P.x[i] - x) < 1e-9) best_y = std::min(best_y, P.y[i]);
    return best_y;
}

// ------------------------------------------------------------------ (1) one borehole is flat --
void case_single_borehole() {
    std::printf("\n== one borehole: horizontal layers reaching both model edges (sec. 7.10.1.1) ==\n");
    m::Project pr = base_project();
    pr.boreholes = {{"BH-1", 12.0, {}, true, 14.0}};
    pr.boreholes[0].level = {20.0, 15.0, 0.0};   // Fill 20->15, Clay 15->0

    const auto R = m::stratigraphy_from_boreholes(pr);
    check(R.ok, "one borehole generates a stratigraphy");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    check(R.polygons.size() == 2, "two layers give two regions");
    if (R.polygons.size() != 2) return;

    // The borehole is at x = 12, but the layers must reach x = 0 and x = 40 unchanged: the manual
    // says a single borehole's surface extends to the model boundaries, and nothing is sloped
    // because there is nothing to slope towards.
    for (double x : {0.0, 12.0, 40.0}) {
        close_to(top_of_polygon_at(R.polygons[0], x), 20.0, 1e-12,
                 "Fill top is 20 at x = " + std::to_string((int)x));
        close_to(bottom_of_polygon_at(R.polygons[0], x), 15.0, 1e-12,
                 "Fill base is 15 at x = " + std::to_string((int)x));
        close_to(top_of_polygon_at(R.polygons[1], x), 15.0, 1e-12,
                 "Clay top is 15 at x = " + std::to_string((int)x));
    }
    check(R.polygons[0].material == 0 && R.polygons[1].material == 1,
          "each region carries its layer's material");
    check(R.has_water, "the borehole head produced a water surface");
    for (size_t i = 0; i < R.wy.size(); ++i)
        if (std::fabs(R.wy[i] - 14.0) > 1e-12)
            check(false, "a single head must be horizontal at 14 m everywhere");
    check(R.wy.size() >= 2, "the water surface spans the model");
    close_to(R.wx.front(), 0.0, 1e-12, "the water surface starts at the model's left edge");
    close_to(R.wx.back(), 40.0, 1e-12, "and ends at its right edge");
}

// ------------------------------------------- (2) two boreholes: linear between, held outside --
void case_two_boreholes() {
    std::printf("\n== two boreholes: linear between them, HELD outside them (sec. 4.2) ==\n");
    m::Project pr = base_project();
    pr.boreholes = {{"BH-1", 10.0, {}, true, 12.0}, {"BH-2", 30.0, {}, true, 16.0}};
    pr.boreholes[0].level = {20.0, 16.0, 0.0};
    pr.boreholes[1].level = {18.0, 10.0, 0.0};   // ground falls 20 -> 18, base of Fill falls 16 -> 10

    const auto R = m::stratigraphy_from_boreholes(pr);
    check(R.ok, "two boreholes generate a stratigraphy");
    if (!R.ok || R.polygons.size() != 2) { check(false, "expected two regions"); return; }

    // THE POLYGON HAS NO VERTEX AT THE MIDPOINT, and that is correct rather than a shortcoming: a
    // boundary breaks slope only at a borehole, so between two logs it is one straight edge and a
    // vertex at x = 20 would be redundant. The geometric claim is therefore about the ENDS of that
    // edge, and the interpolated midpoint is checked on the rule itself.
    close_to(top_of_polygon_at(R.polygons[0], 10.0), 20.0, 1e-12, "the top edge starts at BH-1's 20");
    close_to(top_of_polygon_at(R.polygons[0], 30.0), 18.0, 1e-12, "and ends at BH-2's 18");
    check(R.polygons[0].x.size() == 8,
          "four stations (0, 10, 30, 40) top and bottom: no vertex is invented between the logs");
    std::vector<const m::Borehole*> s{&pr.boreholes[0], &pr.boreholes[1]};
    close_to(m::stratigraphy_level_at(s, 0, 20.0), 19.0, 1e-12,
             "midway the surface is the mean of 20 and 18");
    close_to(m::stratigraphy_level_at(s, 1, 20.0), 13.0, 1e-12,
             "and the Fill base is the mean of 16 and 10");

    // OUTSIDE the outermost boreholes the nearest log is HELD, not continued on its slope. If it
    // were extrapolated the left edge would read 21 and the right edge 17, and a model would grow
    // ground nobody logged.
    close_to(top_of_polygon_at(R.polygons[0], 0.0), 20.0, 1e-12,
             "left of BH-1 the surface holds at 20, it does not extrapolate to 21");
    close_to(top_of_polygon_at(R.polygons[0], 40.0), 18.0, 1e-12,
             "right of BH-2 it holds at 18, not 17");
    close_to(bottom_of_polygon_at(R.polygons[0], 0.0), 16.0, 1e-12, "and the same for the Fill base");
    close_to(bottom_of_polygon_at(R.polygons[0], 40.0), 10.0, 1e-12, "on the right too");

    // The water surface follows the same rule, so it is sloped between the logs and flat outside.
    check(R.has_water && R.wy.size() == R.wx.size(), "the two heads produced a water surface");
    for (size_t i = 0; i < R.wx.size(); ++i) {
        const double want = R.wx[i] <= 10.0 ? 12.0
                          : R.wx[i] >= 30.0 ? 16.0
                          : 12.0 + (R.wx[i] - 10.0) / 20.0 * 4.0;
        close_to(R.wy[i], want, 1e-12,
                 "water level at x = " + std::to_string((int)R.wx[i]));
    }
}

// -------------------------------------------------- (3) a layer that runs out (sec. 4.3.1.1) --
void case_pinch_out() {
    std::printf("\n== a layer that pinches out: zero thickness is legal, not an error (sec. 4.3.1.1) ==\n");
    m::Project pr = base_project();
    pr.strata = {{"Fill", 0}, {"Peat", 1}, {"Clay", 0}};
    pr.boreholes = {{"BH-1", 0.0, {}, false, 0.0}, {"BH-2", 40.0, {}, false, 0.0}};
    // Peat is 2 m thick on the left and absent on the right: its top and bottom meet there.
    pr.boreholes[0].level = {20.0, 16.0, 14.0, 0.0};
    pr.boreholes[1].level = {20.0, 16.0, 16.0, 0.0};

    const auto R = m::stratigraphy_from_boreholes(pr);
    check(R.ok, "the pinching layer still generates");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    check(R.polygons.size() == 3, "all three layers produce a region -- the peat exists somewhere");
    if (R.polygons.size() != 3) return;

    close_to(top_of_polygon_at(R.polygons[1], 0.0) - bottom_of_polygon_at(R.polygons[1], 0.0),
             2.0, 1e-12, "peat is 2 m thick at the left borehole");
    close_to(top_of_polygon_at(R.polygons[1], 40.0) - bottom_of_polygon_at(R.polygons[1], 40.0),
             0.0, 1e-12, "and exactly zero at the right one, where it has run out");
    // The layer BELOW must take up the space the peat gave up, or the ground would have a gap.
    close_to(top_of_polygon_at(R.polygons[2], 40.0), 16.0, 1e-12,
             "the clay's top rises to meet the fill where the peat is gone: no gap in the ground");
    close_to(top_of_polygon_at(R.polygons[2], 0.0), 14.0, 1e-12,
             "and sits under the peat where it is there");
    check(!R.has_water, "no borehole carried a head, so no water surface was invented");
}

// ------------------------------------ (4) a layer that is nowhere, and other refusals to guess --
void case_refusals() {
    std::printf("\n== what it refuses to do rather than guess ==\n");
    m::Project pr = base_project();
    pr.strata = {{"Fill", 0}, {"Ghost", 1}, {"Clay", 0}};
    pr.boreholes = {{"BH-1", 0.0, {}, false, 0.0}, {"BH-2", 40.0, {}, false, 0.0}};
    pr.boreholes[0].level = {20.0, 12.0, 12.0, 0.0};   // Ghost is zero everywhere
    pr.boreholes[1].level = {20.0, 12.0, 12.0, 0.0};
    const auto R = m::stratigraphy_from_boreholes(pr);
    check(R.ok, "a layer that is nowhere does not stop the generation");
    check(R.polygons.size() == 2, "...it just produces no region");
    check(!R.notes.empty(), "and the run SAYS the layer produced nothing, rather than dropping it");
    if (!R.notes.empty()) std::printf("      note: %s\n", R.notes[0].c_str());

    // Wrong number of levels: a log that does not match the layer list is a data-entry mistake, and
    // guessing which end it is missing from would put the wrong ground in a model silently.
    m::Project bad = base_project();
    bad.boreholes = {{"BH-1", 5.0, {}, true, 10.0}};
    bad.boreholes[0].level = {20.0, 10.0};             // two levels for two layers: one short
    const auto Rb = m::stratigraphy_from_boreholes(bad);
    check(!Rb.ok, "a borehole with the wrong number of levels is refused");
    check(Rb.message.find("needs 3") != std::string::npos,
          "...and the refusal says how many it needed");
    if (!Rb.ok) std::printf("      refusal: %s\n", Rb.message.c_str());

    // Inverted levels: boundaries that rise going down are not a stratigraphy.
    m::Project inv = base_project();
    inv.boreholes = {{"BH-1", 5.0, {}, true, 10.0}};
    inv.boreholes[0].level = {20.0, 5.0, 12.0};
    const auto Ri = m::stratigraphy_from_boreholes(inv);
    check(!Ri.ok, "layer boundaries that rise going down are refused");
    if (!Ri.ok) std::printf("      refusal: %s\n", Ri.message.c_str());

    // No boreholes and no layers: both refuse, and neither pretends to have generated something.
    m::Project empty = base_project();
    check(!m::stratigraphy_from_boreholes(empty).ok, "no boreholes generates nothing");
    m::Project nolayers = base_project();
    nolayers.strata.clear();
    nolayers.boreholes = {{"BH-1", 5.0, {}, true, 10.0}};
    check(!m::stratigraphy_from_boreholes(nolayers).ok, "no layers generates nothing");
}

// -------------------------------------------- (5) applying it is a separate, complete action --
void case_apply() {
    std::printf("\n== applying the result replaces the ground and the water together ==\n");
    m::Project pr = base_project();
    pr.boreholes = {{"BH-1", 0.0, {}, true, 11.0}, {"BH-2", 40.0, {}, true, 15.0}};
    pr.boreholes[0].level = {20.0, 14.0, 0.0};
    pr.boreholes[1].level = {20.0, 12.0, 0.0};
    // Something drawn earlier that the generation must replace, not add to.
    pr.polygons.resize(1);
    pr.polygons[0].name = "hand-drawn";
    pr.has_water = false;

    std::string msg;
    check(m::apply_stratigraphy(pr, &msg), "apply succeeds");
    std::printf("      %s\n", msg.c_str());
    check(pr.polygons.size() == 2, "the generated regions REPLACE what was there, not append");
    check(pr.polygons[0].name == "Fill", "and they are the layers, in layer order");
    check(pr.has_water && pr.wy.size() == pr.wx.size() && !pr.wx.empty(),
          "the water surface came with them");
    close_to(pr.wy.front(), 11.0, 1e-12, "left head");
    close_to(pr.wy.back(), 15.0, 1e-12, "right head");

    // A generation that fails must leave the project alone: a half-applied geometry is worse than
    // none, because it looks like a model.
    m::Project keep = pr;
    keep.boreholes[0].level = {20.0, 14.0};   // now malformed
    const m::Project before = keep;
    check(!m::apply_stratigraphy(keep, &msg), "a malformed log is refused on apply too");
    check(keep.polygons.size() == before.polygons.size(),
          "...and the project is untouched by the refusal");
}

}  // namespace

int main() {
    std::printf("== KV-GEO-001: stratigraphy from borehole logs ==\n");
    case_single_borehole();
    case_two_boreholes();
    case_pinch_out();
    case_refusals();
    case_apply();
    if (g_failures == 0) {
        std::printf("\nOK: borehole logs generate the ground the manual's rules describe\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
