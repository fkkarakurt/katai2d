#pragma once
// Boreholes -> soil polygons and a water surface.
//
// The one job in this header is the interpolation rule stated in Project's Borehole comment, taken
// from the PLAXIS 2D 2025.1 Reference Manual sec. 4.2 / 4.3.1.1 / 7.10.1.1. It is a pure function
// from schema to schema: no mesh, no engine, no file. That is deliberate -- the generator is the
// part most likely to be wrong in a way a solve cannot detect, so it is the part that has to be
// testable on its own, against levels a person can compute by hand.
//
// WHAT COMES OUT. One polygon per stratum, in the layer order, each a closed ring: the top boundary
// left to right, then the bottom boundary right to left. The vertices sit at STATIONS -- the model
// edges and every borehole x -- because the boundaries are piecewise linear with a break at each
// borehole and nowhere else. A stratum that is everywhere zero-thickness produces no polygon at
// all; one that pinches out locally produces a polygon that touches itself at that station, which
// is what a layer running out actually looks like and what the manual explicitly allows.
//
// WHAT DOES NOT COME OUT: materials are carried, but nothing else about the polygons is invented.
// Edge boundary conditions, flow conditions and coarseness stay at their defaults, because a
// borehole log says where the ground changes and says nothing about how the model is supported.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <katai/model/project.hpp>

namespace katai::model {

// The x positions at which a generated boundary may change slope: the model edges and every
// borehole. Sorted, de-duplicated, and always at least the two edges.
inline std::vector<double> stratigraphy_stations(const Project& pr) {
    std::vector<double> xs{pr.x_min, pr.x_max};
    for (const Borehole& b : pr.boreholes)
        if (b.x > pr.x_min && b.x < pr.x_max) xs.push_back(b.x);
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(),
                         [](double a, double b) { return std::fabs(a - b) < 1e-9; }),
             xs.end());
    return xs;
}

// Level of boundary `k` at position x, by the manual's rule: linear between the two bracketing
// boreholes, and HELD -- not extrapolated -- outside the outermost ones.
//
// `sorted` must be the boreholes ordered by x. Boreholes whose level vector is too short for `k`
// are skipped rather than read out of range: a malformed log must not decide a geometry, and the
// validator refuses it long before this is called.
inline double stratigraphy_level_at(const std::vector<const Borehole*>& sorted, size_t k, double x) {
    std::vector<const Borehole*> usable;
    usable.reserve(sorted.size());
    for (const Borehole* b : sorted)
        if (k < b->level.size()) usable.push_back(b);
    if (usable.empty()) return 0.0;
    if (usable.size() == 1 || x <= usable.front()->x) return usable.front()->level[k];
    if (x >= usable.back()->x) return usable.back()->level[k];
    for (size_t i = 0; i + 1 < usable.size(); ++i) {
        const Borehole* a = usable[i];
        const Borehole* b = usable[i + 1];
        if (x >= a->x && x <= b->x) {
            const double span = b->x - a->x;
            if (span <= 1e-12) return a->level[k];
            const double t = (x - a->x) / span;
            return a->level[k] + t * (b->level[k] - a->level[k]);
        }
    }
    return usable.back()->level[k];
}

// The phreatic level at x, from the boreholes that carry a head. Same rule as a layer boundary --
// sec. 7.10.1.1 combines the heads into one surface, and a single head is horizontal to the edges.
inline double stratigraphy_head_at(const std::vector<const Borehole*>& sorted, double x) {
    std::vector<const Borehole*> h;
    for (const Borehole* b : sorted)
        if (b->has_head) h.push_back(b);
    if (h.empty()) return 0.0;
    if (h.size() == 1 || x <= h.front()->x) return h.front()->head;
    if (x >= h.back()->x) return h.back()->head;
    for (size_t i = 0; i + 1 < h.size(); ++i) {
        if (x >= h[i]->x && x <= h[i + 1]->x) {
            const double span = h[i + 1]->x - h[i]->x;
            if (span <= 1e-12) return h[i]->head;
            const double t = (x - h[i]->x) / span;
            return h[i]->head + t * (h[i + 1]->head - h[i]->head);
        }
    }
    return h.back()->head;
}

struct StratigraphyResult {
    std::vector<SoilPolygon> polygons;   // one per stratum that has thickness somewhere
    bool has_water = false;              // any borehole carried a head
    std::vector<double> wx, wy;          // the generated phreatic polyline
    std::vector<std::string> notes;      // strata that produced nothing, and why
    bool ok = false;
    std::string message;
};

// Generate the geometry a set of borehole logs describes. Pure: `pr` is not modified, and applying
// the result is the caller's separate decision.
inline StratigraphyResult stratigraphy_from_boreholes(const Project& pr) {
    StratigraphyResult R;
    if (pr.boreholes.empty()) {
        R.message = "No boreholes: there is nothing to generate a stratigraphy from.";
        return R;
    }
    if (pr.strata.empty()) {
        R.message = "No soil layers defined: a borehole records where the layers are, so at least "
                    "one layer has to exist first.";
        return R;
    }
    const size_t n_lv = pr.strata.size() + 1;
    for (const Borehole& b : pr.boreholes) {
        if (b.level.size() != n_lv) {
            R.message = "Borehole '" + b.name + "' carries " + std::to_string(b.level.size()) +
                        " level(s) for " + std::to_string(pr.strata.size()) +
                        " layer(s); it needs " + std::to_string(n_lv) +
                        " (a top for every layer, then the base of the lowest).";
            return R;
        }
        for (size_t k = 0; k + 1 < b.level.size(); ++k)
            if (b.level[k] < b.level[k + 1] - 1e-9) {
                R.message = "Borehole '" + b.name + "': layer boundaries must not rise going down "
                            "(level " + std::to_string(k) + " is below level " +
                            std::to_string(k + 1) + "). Equal values are a layer that pinches out "
                            "here, which is allowed.";
                return R;
            }
    }

    std::vector<const Borehole*> sorted;
    sorted.reserve(pr.boreholes.size());
    for (const Borehole& b : pr.boreholes) sorted.push_back(&b);
    std::sort(sorted.begin(), sorted.end(),
              [](const Borehole* a, const Borehole* b) { return a->x < b->x; });

    const std::vector<double> xs = stratigraphy_stations(pr);

    for (size_t j = 0; j < pr.strata.size(); ++j) {
        std::vector<double> top(xs.size()), bot(xs.size());
        double thickest = 0.0;
        for (size_t i = 0; i < xs.size(); ++i) {
            top[i] = stratigraphy_level_at(sorted, j, xs[i]);
            bot[i] = stratigraphy_level_at(sorted, j + 1, xs[i]);
            thickest = std::max(thickest, top[i] - bot[i]);
        }
        // A layer with no thickness ANYWHERE is not a degenerate polygon to be cleaned up later --
        // it is a row the user filled in and then zeroed, and saying so is cheaper than leaving a
        // sliver in the geometry for the mesher to trip over.
        if (thickest <= 1e-9) {
            R.notes.push_back("Layer '" + pr.strata[j].name +
                              "' has zero thickness at every borehole, so it produced no region.");
            continue;
        }
        SoilPolygon P;
        P.name = pr.strata[j].name;
        P.material = pr.strata[j].material;
        const auto push = [&](double x, double y) {
            if (!P.x.empty() && std::fabs(P.x.back() - x) < 1e-9 && std::fabs(P.y.back() - y) < 1e-9)
                return;                                  // never emit the same vertex twice
            P.x.push_back(x);
            P.y.push_back(y);
        };
        for (size_t i = 0; i < xs.size(); ++i) push(xs[i], top[i]);
        for (size_t i = xs.size(); i-- > 0;) push(xs[i], bot[i]);
        // Closing vertex duplicates the first (the ring is implicitly closed).
        if (P.x.size() >= 2 && std::fabs(P.x.front() - P.x.back()) < 1e-9 &&
            std::fabs(P.y.front() - P.y.back()) < 1e-9) {
            P.x.pop_back();
            P.y.pop_back();
        }
        if (P.x.size() < 3) {
            R.notes.push_back("Layer '" + pr.strata[j].name +
                              "' collapsed to fewer than three distinct corners and produced no "
                              "region.");
            continue;
        }
        P.edge_bc.assign(P.x.size(), 0);
        R.polygons.push_back(std::move(P));
    }

    for (const Borehole* b : sorted)
        if (b->has_head) { R.has_water = true; break; }
    if (R.has_water)
        for (double x : xs) {
            R.wx.push_back(x);
            R.wy.push_back(stratigraphy_head_at(sorted, x));
        }

    if (R.polygons.empty()) {
        R.message = "Every layer came out with zero thickness, so no ground was generated.";
        return R;
    }
    R.ok = true;
    R.message = "Generated " + std::to_string(R.polygons.size()) + " region(s) from " +
                std::to_string(pr.boreholes.size()) + " borehole(s).";
    return R;
}

// Replace the project's ground with what the boreholes describe. Separate from the generator above
// so that a caller can show the result before committing to it, and so the test can check the
// generation without a Project to write into.
inline bool apply_stratigraphy(Project& pr, std::string* message = nullptr) {
    const StratigraphyResult R = stratigraphy_from_boreholes(pr);
    if (message) *message = R.message;
    if (!R.ok) return false;
    pr.polygons = R.polygons;
    if (R.has_water) {
        pr.has_water = true;
        pr.wx = R.wx;
        pr.wy = R.wy;
    }
    return true;
}

}  // namespace katai::model
