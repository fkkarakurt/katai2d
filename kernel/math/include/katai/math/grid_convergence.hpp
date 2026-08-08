#pragma once
// Discretisation-error estimation from a systematically refined mesh triplet: the observed
// order of convergence, the Richardson-extrapolated value at h -> 0, and the Grid Convergence
// Index (GCI) that turns the two into a reportable numerical uncertainty band.
//
// WHY THIS IS HERE, AND WHERE IT COMES FROM
// ----------------------------------------
// A benchmark reported at a single mesh density -- "Cox +3.7% on the file's own 0.25 m mesh" --
// leaves a reader unable to separate the MODEL's error from the MESH's. The geotechnical finite
// element codes this program is measured against do not answer that question either: PLAXIS,
// MIDAS GTS NX and GEO5 all provide mesh refinement CONTROLS and qualitative advice ("finer is
// more accurate; balance accuracy against run time"), and PLAXIS's Scientific Manual chapter 9
// governs ITERATIVE convergence -- when a Newton loop may stop -- which is a different question
// from how far a converged answer sits from the exact solution of the same equations. No
// vendor manual in this field defines a discretisation-error estimator.
//
// The procedure implemented here therefore comes from the verification literature, where it is
// standardised:
//
//   [1] P. J. Roache (1994), "Perspective: A Method for Uniform Reporting of Grid Refinement
//       Studies", ASME J. Fluids Eng. 116(3):405-413 -- the Grid Convergence Index and its
//       factor of safety (1.25 when three or more grids establish the observed order, 3.0 when
//       only two grids are available and the order must be assumed).
//   [2] I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad (2008),
//       "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD
//       Applications", ASME J. Fluids Eng. 130(7):078001 -- the journal's editorial policy
//       statement; the five-step procedure and the fixed-point form of the observed order that
//       admits UNEQUAL refinement ratios, which is the case an unstructured mesher produces.
//   [3] ASME V&V 20-2009 (R2016), "Standard for Verification and Validation in Computational
//       Fluid Dynamics and Heat Transfer" -- adopts a GCI-based numerical uncertainty; and
//       ASME V&V 10-2006, "Guide for Verification and Validation in Computational Solid
//       Mechanics", the companion whose scope actually contains this program. Their vocabulary
//       is the one used here: CODE verification asks whether the code solves the equations it
//       claims to; SOLUTION verification -- this file -- estimates the numerical accuracy of
//       one particular calculation.
//
// The equations below were re-derived before they were coded, because the secondary literature
// transcribes them inconsistently (one widely-read restatement prints the extrapolated value as
// r^p (phi1 - phi2)/(r^p - 1), which is not what Richardson extrapolation gives). Writing
// phi_i = phi_exact + C h_i^p and eliminating C:
//
//   eps21 = phi2 - phi1 = C h1^p (r21^p - 1)
//   eps32 = phi3 - phi2 = C h1^p r21^p (r32^p - 1)
//   => ln|eps32/eps21| = p ln r21 + ln|(r32^p - 1)/(r21^p - 1)|
//   => p = [ ln|eps32/eps21| + q(p) ] / ln r21,  q(p) = ln((r21^p - s)/(r32^p - s))
//   => phi_ext = phi1 - C h1^p = (r21^p phi1 - phi2)/(r21^p - 1)
//
// which is [2] exactly, and which tests/test_grid_convergence.cpp pins against manufactured
// data whose exact answer is known in closed form. The sign s = sign(eps32/eps21) is [2]'s
// treatment of oscillatory convergence (after Celik and Karatekin); it is an accepted
// convention rather than a derived result, and the literature says so -- which is why an
// oscillatory triplet is REPORTED as such here instead of being quietly extrapolated.

#include <cmath>
#include <string>

namespace katai::math {

// How the three solutions behave as the mesh is refined, from R = eps21/eps32 ([2]; the same
// classification appears in the applied literature as the "convergence condition"). Only
// MonotonicConvergence justifies quoting an extrapolated value without qualification.
enum class ConvergenceKind {
    MonotonicConvergence,    // 0 < R < 1: each refinement moves the answer less than the last
    MonotonicDivergence,     // R >= 1: refinement moves the answer MORE each time
    OscillatoryConvergence,  // -1 < R < 0: the answer alternates about the limit, narrowing
    OscillatoryDivergence,   // R <= -1: alternates and widens
    Exact                    // eps21 = 0: the quantity does not change with the mesh at all
};

inline const char* convergence_kind_name(ConvergenceKind k) {
    switch (k) {
        case ConvergenceKind::MonotonicConvergence:   return "monotonic convergence";
        case ConvergenceKind::MonotonicDivergence:    return "monotonic divergence";
        case ConvergenceKind::OscillatoryConvergence: return "oscillatory convergence";
        case ConvergenceKind::OscillatoryDivergence:  return "oscillatory divergence";
        case ConvergenceKind::Exact:                  return "mesh-independent";
    }
    return "unknown";
}

// Three solutions of the SAME quantity on three systematically refined meshes, 1 = finest.
// `h` is the representative cell size: for a 2D mesh, sqrt(total area / element count)
// ([2], step 1), which is what mesh_representative_size() below computes.
struct GridTriplet {
    double h1 = 0.0, h2 = 0.0, h3 = 0.0;          // h1 < h2 < h3
    double phi1 = 0.0, phi2 = 0.0, phi3 = 0.0;    // the quantity on each
};

struct ConvergenceEstimate {
    bool ok = false;              // false: the triplet cannot support an estimate; see `message`
    std::string message;          // why not, or a qualification the reader must carry
    ConvergenceKind kind = ConvergenceKind::MonotonicConvergence;
    double r21 = 0.0, r32 = 0.0;  // refinement factors h2/h1, h3/h2
    double p = 0.0;               // observed ("apparent") order of convergence
    double phi_extrapolated = 0.0;// Richardson value at h -> 0
    double e_approx = 0.0;        // e_a^21, approximate relative error   (fraction, not %)
    double e_extrapolated = 0.0;  // e_ext^21, relative to the extrapolated value
    double gci_fine = 0.0;        // GCI^21_fine from the OBSERVED order  (fraction, not %)
    double gci_medium = 0.0;      // GCI^32
    double asymptotic_ratio = 0.0;// GCI^32 / (r21^p GCI^21): -> 1 inside the asymptotic range
    int iterations = 0;           // fixed-point iterations spent on p

    // ---- what to publish -------------------------------------------------------------------
    // The observed order is only meaningful when the triplet is in the asymptotic range. When
    // it is not -- and on an unstructured mesher that rebuilds the mesh from scratch at every
    // density, it often is not -- p comes out far from anything the elements can deliver, and
    // an extrapolation built on it is worse than no extrapolation at all: 1/(r^p - 1) grows
    // without bound as p -> 0, so a tiny mesh-to-mesh wobble is amplified into a large,
    // confident-looking correction. `band` is therefore the number to report, and `band_basis`
    // says how it was obtained.
    bool asymptotic = false;      // observed order accepted: extrapolation may be quoted
    double p_used = 0.0;          // the order `band` was computed with (observed or assumed)
    double safety_used = 0.0;     // the factor of safety `band` was computed with
    double band = 0.0;            // THE numerical uncertainty to publish  (fraction, not %)
    std::string band_basis;       // "observed order" / "assumed order (Fs = 3)" + the reason
};

// Representative cell size of a 2D mesh, [2] step 1: h = sqrt((1/N) sum_i dA_i). Passing the
// total area and the element count keeps this header free of any mesh type.
inline double mesh_representative_size(double total_area, int element_count) {
    return element_count > 0 && total_area > 0.0
               ? std::sqrt(total_area / (double)element_count)
               : 0.0;
}

// The estimate. `safety_factor` is Roache's Fs: 1.25 is the value [1] and [2] prescribe when
// three grids establish the observed order, and it is the default because that is what this
// program's sweeps do; it is exposed because a two-grid study must use 3.0 and say so.
//
// Everything is returned as a FRACTION (0.022 = 2.2%). No exceptions are thrown and no case is
// silently patched: a triplet that cannot carry an estimate comes back ok = false with the
// reason, because a numerical uncertainty that was quietly invented is worse than none.
// The formal order this program's elements can deliver for a recovered field, and the window
// outside which an observed order is evidence that the triplet is NOT in the asymptotic range
// rather than evidence about the discretisation. Both are stated here so that a report can
// quote them, and both are conventions of this program -- see the note on `band` below.
struct OrderPolicy {
    double p_assumed = 2.0;       // used when the observed order is rejected
    double p_min = 0.5;           // below this, 1/(r^p - 1) amplifies mesh noise into "signal"
    double p_max = 3.0;           // above this, no element here can be responsible for it
    double asymptotic_tol = 0.1;  // |asymptotic_ratio - 1| must be within this
    double safety_observed = 1.25;// Roache's Fs when three grids establish the order
    double safety_assumed = 3.0;  // Roache's Fs when the order is assumed rather than measured
};

inline ConvergenceEstimate estimate_grid_convergence(const GridTriplet& g,
                                                     double safety_factor = 1.25) {
    ConvergenceEstimate out;
    if (!(g.h1 > 0.0) || !(g.h2 > g.h1) || !(g.h3 > g.h2)) {
        out.message = "the three meshes must be ordered fine to coarse with positive sizes";
        return out;
    }
    out.r21 = g.h2 / g.h1;
    out.r32 = g.h3 / g.h2;

    const double eps21 = g.phi2 - g.phi1;
    const double eps32 = g.phi3 - g.phi2;
    if (eps21 == 0.0) {
        // The quantity did not move between the two finest meshes. That is not a failure -- it
        // is the strongest possible result -- but the order of convergence is then undefined
        // (0/0), so it is reported as what it is rather than as a number.
        out.ok = true;
        out.kind = ConvergenceKind::Exact;
        out.phi_extrapolated = g.phi1;
        out.message = "the finest two meshes give the same value; no discretisation error is "
                      "measurable and the observed order is undefined";
        return out;
    }

    const double R = eps21 / eps32;
    out.kind = eps32 == 0.0             ? ConvergenceKind::MonotonicDivergence
             : R >= 1.0                 ? ConvergenceKind::MonotonicDivergence
             : R > 0.0                  ? ConvergenceKind::MonotonicConvergence
             : R > -1.0                 ? ConvergenceKind::OscillatoryConvergence
                                        : ConvergenceKind::OscillatoryDivergence;
    if (eps32 == 0.0) {
        out.message = "the two coarsest meshes give the same value; the triplet carries no "
                      "information about the order of convergence";
        return out;
    }

    // Observed order, [2]: p = |ln|eps32/eps21| + q(p)| / ln(r21), solved by fixed-point
    // iteration from q = 0. s = sign(eps32/eps21) is 1 for monotone behaviour and -1 for
    // oscillatory behaviour.
    const double ratio = eps32 / eps21;
    const double s = ratio > 0.0 ? 1.0 : -1.0;
    const double ln_r21 = std::log(out.r21);
    const double ln_ratio = std::log(std::fabs(ratio));
    double p = std::fabs(ln_ratio) / ln_r21;    // q = 0 start, i.e. the equal-ratio closed form
    int it = 0;
    for (; it < 200; ++it) {
        const double dr21 = std::pow(out.r21, p) - s;
        const double dr32 = std::pow(out.r32, p) - s;
        if (!(dr21 > 0.0) || !(dr32 > 0.0)) {   // s = 1 with p -> 0 leaves q undefined
            out.message = "the observed order cannot be evaluated: r^p - s is not positive, "
                          "which happens when the refinement is too small to resolve a trend";
            return out;
        }
        const double p_next = std::fabs(ln_ratio + std::log(dr21 / dr32)) / ln_r21;
        const double dp = std::fabs(p_next - p);
        p = p_next;
        if (dp < 1e-13 * std::fmax(1.0, p)) { ++it; break; }
    }
    if (!(p > 0.0) || !std::isfinite(p)) {
        out.message = "the observed order did not resolve to a positive finite value";
        return out;
    }
    out.p = p;
    out.iterations = it;

    const double r21p = std::pow(out.r21, p);
    const double r32p = std::pow(out.r32, p);
    if (!(std::fabs(r21p - 1.0) > 0.0) || !(std::fabs(r32p - 1.0) > 0.0)) {
        out.message = "r^p = 1: the meshes are not distinguishable at the observed order";
        return out;
    }
    out.phi_extrapolated = (r21p * g.phi1 - g.phi2) / (r21p - 1.0);

    // Relative errors are taken against the fine solution ([2]) -- with a floor, because a
    // quantity whose value is legitimately near zero (a heave, a residual) would otherwise
    // report an infinite relative error and hide the real, small, absolute one.
    const double scale1 = std::fmax(std::fabs(g.phi1), 1e-300);
    const double scale2 = std::fmax(std::fabs(g.phi2), 1e-300);
    out.e_approx = std::fabs(eps21) / scale1;
    out.e_extrapolated = std::fabs((out.phi_extrapolated - g.phi1) /
                                   std::fmax(std::fabs(out.phi_extrapolated), 1e-300));
    out.gci_fine = safety_factor * out.e_approx / std::fabs(r21p - 1.0);
    out.gci_medium = safety_factor * (std::fabs(eps32) / scale2) / std::fabs(r32p - 1.0);
    // Asymptotic-range check: with an exact power law the two indices differ by exactly r^p,
    // so this ratio tends to 1 as the triplet enters the asymptotic range. It is a diagnostic,
    // not a gate -- a value far from 1 means the triplet is too coarse for the extrapolation
    // to be trusted, and that is a statement the report must carry rather than suppress.
    out.asymptotic_ratio = out.gci_fine > 0.0 ? out.gci_medium / (r21p * out.gci_fine) : 0.0;

    out.ok = true;
    // Qualifications the reader must carry with the number. [2] recommends a refinement factor
    // of at least 1.3; below that the observed order is dominated by round-off and by the
    // mesher's own irregularity rather than by the discretisation.
    if (out.r21 < 1.3 || out.r32 < 1.3)
        out.message = "refinement factor below the recommended 1.3: the observed order is not "
                      "reliable and the band should be read as indicative";
    else if (out.kind == ConvergenceKind::OscillatoryConvergence)
        out.message = "oscillatory convergence: the extrapolation follows the accepted "
                      "convention (s = -1) rather than a derived result";
    else if (out.kind != ConvergenceKind::MonotonicConvergence)
        out.message = std::string("the triplet shows ") + convergence_kind_name(out.kind) +
                      ": the extrapolated value is not a limit and must not be quoted as one";
    return out;
}

// Decide WHICH band to publish, and say why. This is the one place in this header that is a
// policy rather than a citation, and it is written out because it has to be auditable.
//
// Roache [1] gives both halves: a factor of safety of 1.25 when three grids ESTABLISH the
// order, and 3.0 when the order is ASSUMED instead of measured. Celik et al. [2] add the
// warning that an observed order far from the formal one indicates the triplet is not in the
// asymptotic range -- their example of a bad sign is an order of 7 for a second-order scheme.
// Neither source states a numeric acceptance window, so the window in OrderPolicy is this
// program's reading, chosen for a reason that can be checked: as p -> 0 the amplification
// 1/(r^p - 1) diverges, so an observed order below p_min turns mesh-to-mesh irregularity into
// a large "correction", and an order above p_max cannot be produced by these elements.
//
// When the observed order is rejected, the band is recomputed with the assumed order and the
// larger factor of safety -- a wider, honestly-labelled band -- and `asymptotic` stays false so
// no caller quotes the extrapolated value as a limit.
inline void apply_order_policy(ConvergenceEstimate& e, const GridTriplet& g,
                               const OrderPolicy& pol = OrderPolicy{}) {
    if (!e.ok) return;
    if (e.kind == ConvergenceKind::Exact) {
        e.asymptotic = true;
        e.p_used = 0.0; e.safety_used = 0.0; e.band = 0.0;
        e.band_basis = "no measurable change between the two finest meshes";
        return;
    }
    const bool monotone = e.kind == ConvergenceKind::MonotonicConvergence;
    const bool p_ok = e.p >= pol.p_min && e.p <= pol.p_max;
    const bool ratio_ok = std::fabs(e.asymptotic_ratio - 1.0) <= pol.asymptotic_tol;
    e.asymptotic = monotone && p_ok && ratio_ok;
    if (e.asymptotic) {
        e.p_used = e.p;
        e.safety_used = pol.safety_observed;
        e.band = e.gci_fine;
        e.band_basis = "observed order (Fs = 1.25)";
        return;
    }
    const double r21p = std::pow(e.r21, pol.p_assumed);
    e.p_used = pol.p_assumed;
    e.safety_used = pol.safety_assumed;
    e.band = std::fabs(r21p - 1.0) > 0.0
                 ? pol.safety_assumed * e.e_approx / std::fabs(r21p - 1.0)
                 : e.e_approx;
    e.band_basis = std::string("assumed order p = ") +
                   (pol.p_assumed == 2.0 ? "2" : "given") + " with Fs = 3, because " +
                   (!monotone   ? std::string("the triplet shows ") + convergence_kind_name(e.kind)
                    : !p_ok     ? "the observed order is outside the plausible window"
                                : "the asymptotic-range check is not met");

    // A floor the evidence sets, applied ONLY outside the asymptotic range. There, the three
    // values are not a hierarchy of errors that the finest one crowns -- they are three answers
    // whose differences are dominated by the mesher's irregularity, and e_a (which sees only
    // the finest pair) can happen to be the smallest of them. A band narrower than the range
    // the family actually produced would claim a reproducibility the data denies.
    //
    // Inside the asymptotic range the same floor would be WRONG: the coarse mesh is expected
    // to be further off, its error is not the fine mesh's uncertainty, and the GCI is the
    // established estimate of the latter. Hence the asymmetry, deliberately.
    const double lo = std::fmin(g.phi1, std::fmin(g.phi2, g.phi3));
    const double hi = std::fmax(g.phi1, std::fmax(g.phi2, g.phi3));
    const double spread = (hi - lo) / std::fmax(std::fabs(g.phi1), 1e-300);
    if (spread > e.band) {
        e.band = spread;
        e.band_basis += "; widened to the spread across the three meshes, which is what they "
                        "actually did";
    }
}

// Convenience: estimate and apply the policy in one call, which is how a report should use it.
inline ConvergenceEstimate grid_convergence_band(const GridTriplet& g,
                                                 const OrderPolicy& pol = OrderPolicy{}) {
    ConvergenceEstimate e = estimate_grid_convergence(g, pol.safety_observed);
    apply_order_policy(e, g, pol);
    return e;
}

}  // namespace katai::math
