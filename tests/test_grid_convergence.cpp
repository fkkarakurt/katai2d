// Solution verification, the instrument itself: does KATAI's implementation of the Grid
// Convergence Index actually solve the equation it claims to?
//
// Every number this program will publish with a numerical-uncertainty band gets that band from
// katai/math/grid_convergence.hpp, so the estimator has to be verified before it is used to
// verify anything else -- and it can be, exactly, because its defining assumption is a closed
// form. Richardson extrapolation assumes
//
//     phi(h) = phi_exact + C h^p
//
// so a triplet MANUFACTURED from that expression has a known exact answer: the estimator must
// return p and phi_exact to round-off, whatever the refinement factors are. That is a stronger
// oracle than reproducing a published worked example, because it tests the equation rather than
// a transcription of a table -- and transcription is the risk here: the secondary literature
// restating Celik et al. disagrees with itself on the extrapolation formula (one widely-read
// version prints r^p (phi1 - phi2)/(r^p - 1), which is not Richardson extrapolation).
//
// The unequal-ratio fixtures are the ones that matter for this program: an unstructured mesher
// asked for three densities delivers three element counts, never three exact ratios, so the
// fixed-point form with q(p) is the form that will actually run.
//
// verify: KV-NUM-004
//   oracle:   closed_form
//   source:   P. J. Roache (1994), "Perspective: A Method for Uniform Reporting of Grid Refinement Studies", ASME J. Fluids Eng. 116(3):405-413 (the Grid Convergence Index and its factor of safety); I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad (2008), "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications", ASME J. Fluids Eng. 130(7):078001 (the procedure and the unequal-ratio form of the observed order); ASME V&V 20-2009 and V&V 10-2006 for the solution-verification framework
//   locator:  phi(h) = phi_exact + C h^p; p = |ln|eps32/eps21| + q(p)| / ln(r21) with q(p) = ln((r21^p - s)/(r32^p - s)) and s = sign(eps32/eps21); phi_ext^21 = (r21^p phi1 - phi2)/(r21^p - 1); e_a^21 = |(phi1 - phi2)/phi1|; GCI_fine^21 = Fs e_a^21/(r21^p - 1) with Fs = 1.25 for three grids (stated in full)
//   quantity: the observed order p [-], the extrapolated value phi_ext [same unit as the quantity], and the fine-grid GCI [-], recovered from triplets manufactured from the closed form above at equal and at unequal refinement factors
//   expected: p = p_true and phi_ext = phi_exact exactly; and, since the data is an exact power law, the identity GCI_fine = Fs |phi1 - phi_exact| / |phi1| holds, i.e. the reported band is exactly the factor of safety times the TRUE relative error of the fine mesh
//   band:     1e-10 relative on p and phi_ext, 1e-12 relative on the GCI identity, as asserted below -- these are round-off tolerances, not physical ones: the manufactured triplet satisfies the estimator's assumption exactly, so any larger deviation is an error in the implementation

#include <katai/math/grid_convergence.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace gc = katai::math;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
void check_close(double got, double want, double rtol, const std::string& what) {
    const double err = std::fabs(got - want) / std::fmax(std::fabs(want), 1e-300);
    std::printf(err <= rtol ? "ok:   %s (%.15g vs %.15g, rel %.2e)\n"
                            : "FAIL: %s (%.15g vs %.15g, rel %.2e)\n",
                what.c_str(), got, want, err);
    if (!(err <= rtol)) ++g_failures;
}

// A triplet that satisfies phi(h) = phi_exact + C h^p exactly.
gc::GridTriplet manufacture(double phi_exact, double C, double p, double h1, double h2,
                            double h3) {
    gc::GridTriplet g;
    g.h1 = h1; g.h2 = h2; g.h3 = h3;
    g.phi1 = phi_exact + C * std::pow(h1, p);
    g.phi2 = phi_exact + C * std::pow(h2, p);
    g.phi3 = phi_exact + C * std::pow(h3, p);
    return g;
}

// ---------------------------------------------------------------------------------------------
// 1. The case this program will actually meet: three meshes from an unstructured mesher, so the
//    refinement factors differ. Only the q(p) fixed point solves this; the equal-ratio closed
//    form would return a wrong order, which is exactly the failure this fixture would catch.
void case_unequal_ratios() {
    std::printf("\n== exact power law, unequal refinement factors ==\n");
    const double phi_exact = 7.5, C = 0.8, p_true = 2.0;
    const gc::GridTriplet g = manufacture(phi_exact, C, p_true, 0.5, 0.8, 1.4);  // r21=1.6, r32=1.75
    const gc::ConvergenceEstimate e = gc::estimate_grid_convergence(g);
    std::printf("      r21 = %.6f, r32 = %.6f, %d fixed-point iterations\n", e.r21, e.r32,
                e.iterations);
    check(e.ok, "the triplet yields an estimate");
    check(e.kind == gc::ConvergenceKind::MonotonicConvergence, "classified as monotonic convergence");
    check_close(e.p, p_true, 1e-10, "observed order recovers p");
    check_close(e.phi_extrapolated, phi_exact, 1e-10, "extrapolation recovers the exact value");

    // The identity that makes the band meaningful: with exact power-law data the fine mesh's
    // TRUE relative error is |phi1 - phi_exact|/|phi1|, and the GCI is exactly Fs times it.
    const double true_rel_error = std::fabs(g.phi1 - phi_exact) / std::fabs(g.phi1);
    check_close(e.gci_fine, 1.25 * true_rel_error, 1e-12,
                "GCI_fine == 1.25 x the TRUE relative error of the fine mesh");
    check(e.gci_fine > true_rel_error, "the band contains the true error, with margin");
    std::printf("      GCI_fine = %.4f%%, true error = %.4f%%, asymptotic ratio = %.6f\n",
                100.0 * e.gci_fine, 100.0 * true_rel_error, e.asymptotic_ratio);
}

// 2. Equal factors must reduce to the textbook closed form p = ln(eps32/eps21)/ln(r): the
//    general implementation may not drift from the special case everyone checks by hand.
void case_equal_ratios_closed_form() {
    std::printf("\n== exact power law, equal refinement factors ==\n");
    const double phi_exact = -0.0421, C = 1.3e-3, p_true = 1.7, r = 1.5;
    const gc::GridTriplet g = manufacture(phi_exact, C, p_true, 0.4, 0.4 * r, 0.4 * r * r);
    const gc::ConvergenceEstimate e = gc::estimate_grid_convergence(g);
    const double closed_form = std::log(std::fabs((g.phi3 - g.phi2) / (g.phi2 - g.phi1))) /
                               std::log(r);
    check(e.ok, "the triplet yields an estimate");
    check_close(e.p, p_true, 1e-10, "observed order recovers p");
    check_close(e.p, closed_form, 1e-12, "and equals the equal-ratio closed form exactly");
    check_close(e.phi_extrapolated, phi_exact, 1e-10,
                "extrapolation recovers the exact value (negative quantity)");
    // The asymptotic-range check is a diagnostic, and its exact value here is phi1/phi2 -- it
    // tends to 1 only as the solutions themselves converge, which is what "asymptotic" means.
    check_close(e.asymptotic_ratio, g.phi1 / g.phi2, 1e-12,
                "asymptotic ratio equals phi1/phi2 for an exact power law");
}

// 3. A non-integer order, and one below the nominal order of the elements. The estimator must
//    report what the data says, not what the formulation promises: a singularity, a coarse
//    mesh or a non-associated flow rule all degrade the observed order, and hiding that behind
//    an assumed p = 2 is how a discretisation error gets under-reported.
void case_non_integer_order() {
    std::printf("\n== observed order is measured, not assumed ==\n");
    const double phi_exact = 233.9, C = -12.0, p_true = 1.37;
    const gc::GridTriplet g = manufacture(phi_exact, C, p_true, 0.25, 0.375, 0.5625);
    const gc::ConvergenceEstimate e = gc::estimate_grid_convergence(g);
    check(e.ok, "the triplet yields an estimate");
    check_close(e.p, p_true, 1e-10, "observed order recovers a non-integer p");
    check_close(e.phi_extrapolated, phi_exact, 1e-10,
                "extrapolation recovers the exact value from below");
}

// 4. Oscillatory data. The sign convention s = -1 is an accepted patch, not a theorem, so the
//    contract here is the CLASSIFICATION and the caveat -- not a precise order.
void case_oscillatory() {
    std::printf("\n== oscillatory convergence is named, not extrapolated silently ==\n");
    const double phi_exact = 2.0, C = 0.05, p = 2.0;
    gc::GridTriplet g;
    g.h1 = 0.5; g.h2 = 0.75; g.h3 = 1.125;
    g.phi1 = phi_exact - C * std::pow(g.h1, p);
    g.phi2 = phi_exact + C * std::pow(g.h2, p);
    g.phi3 = phi_exact - C * std::pow(g.h3, p);
    const gc::ConvergenceEstimate e = gc::estimate_grid_convergence(g);
    check(e.kind == gc::ConvergenceKind::OscillatoryConvergence, "classified as oscillatory");
    check(!e.message.empty(), "and carries the qualification with it");
    std::printf("      message: %s\n", e.message.c_str());
}

// 5. Divergence under refinement. This is the answer to "our answer got worse on a finer mesh",
//    and the one case where an extrapolated value must never be quoted as a limit.
void case_divergence() {
    std::printf("\n== divergence under refinement is refused as a limit ==\n");
    gc::GridTriplet g;
    g.h1 = 0.5; g.h2 = 0.75; g.h3 = 1.125;
    g.phi1 = 10.0; g.phi2 = 10.4; g.phi3 = 10.5;   // eps21 = 0.4 > eps32 = 0.1 -> R = 4
    const gc::ConvergenceEstimate e = gc::estimate_grid_convergence(g);
    check(e.kind == gc::ConvergenceKind::MonotonicDivergence, "classified as monotonic divergence");
    check(e.message.find("not a limit") != std::string::npos,
          "the message forbids quoting the extrapolated value");
    std::printf("      message: %s\n", e.message.c_str());
}

// 6. Degenerate and abusive input: no exceptions, no invented numbers.
void case_degenerate() {
    std::printf("\n== degenerate triplets are refused, not patched ==\n");
    gc::GridTriplet same;
    same.h1 = 0.5; same.h2 = 0.75; same.h3 = 1.125;
    same.phi1 = same.phi2 = same.phi3 = 4.2;
    const gc::ConvergenceEstimate e0 = gc::estimate_grid_convergence(same);
    check(e0.ok && e0.kind == gc::ConvergenceKind::Exact,
          "identical solutions: mesh-independent, and said so");
    check_close(e0.phi_extrapolated, 4.2, 1e-15, "and the extrapolated value is that value");
    check(e0.gci_fine == 0.0, "with a zero band");

    gc::GridTriplet unordered = manufacture(1.0, 0.1, 2.0, 0.8, 0.5, 1.4);   // h2 < h1
    const gc::ConvergenceEstimate e1 = gc::estimate_grid_convergence(unordered);
    check(!e1.ok, "meshes out of order: refused");
    check(!e1.message.empty(), "with a reason");

    // Refinement below the 1.3 that Celik et al. recommend: the estimate is produced, and it
    // arrives with the warning that makes it quotable only as indicative.
    const gc::GridTriplet timid = manufacture(5.0, 0.2, 2.0, 1.0, 1.1, 1.21);
    const gc::ConvergenceEstimate e2 = gc::estimate_grid_convergence(timid);
    check(e2.ok, "timid refinement still estimates");
    check(e2.message.find("1.3") != std::string::npos,
          "and states that the refinement factor is below the recommended 1.3");
}

// 7. The two-grid case has to use Roache's larger factor of safety, and the caller has to be
//    able to say so: the parameter exists precisely so that a wider band is a stated choice.
void case_safety_factor() {
    std::printf("\n== the factor of safety is the caller's declaration ==\n");
    const gc::GridTriplet g = manufacture(7.5, 0.8, 2.0, 0.5, 0.8, 1.4);
    const gc::ConvergenceEstimate a = gc::estimate_grid_convergence(g, 1.25);
    const gc::ConvergenceEstimate b = gc::estimate_grid_convergence(g, 3.0);
    check_close(b.gci_fine / a.gci_fine, 3.0 / 1.25, 1e-14,
                "the band scales exactly with the factor of safety");
    check_close(a.p, b.p, 1e-15, "and the observed order does not depend on it");
}

// 8. The representative cell size of [2], step 1 -- the only place a mesh enters this header.
void case_representative_size() {
    std::printf("\n== representative cell size ==\n");
    // 800 m2 of soil in 3200 elements: sqrt(800/3200) = 0.5 m.
    check_close(gc::mesh_representative_size(800.0, 3200), 0.5, 1e-15,
                "h = sqrt(total area / element count)");
    check(gc::mesh_representative_size(0.0, 10) == 0.0, "no area: no size, rather than a guess");
    check(gc::mesh_representative_size(10.0, 0) == 0.0, "no elements: same");
}

}  // namespace

int main() {
    std::printf("== grid convergence index (KV-NUM-004) ==\n");
    case_unequal_ratios();
    case_equal_ratios_closed_form();
    case_non_integer_order();
    case_oscillatory();
    case_divergence();
    case_degenerate();
    case_safety_factor();
    case_representative_size();
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
