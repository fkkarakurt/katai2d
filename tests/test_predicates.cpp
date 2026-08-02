// Verification of the robust orientation predicate (P1.4a).
//
// Three properties matter for a crash-proof mesher: (i) correct sign on clear
// configurations, (ii) exactly zero on degenerate (collinear) input, and -- the
// reason an adaptive predicate exists at all -- (iii) the correct sign even when
// the points are so nearly collinear that a naive double determinant rounds to
// the wrong sign.
#include <katai/mesh/geom_predicates.hpp>

#include <cmath>
#include <cstdio>

using katai::mesh::incircle;
using katai::mesh::orient2d;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
int sign(double v) { return (v > 0.0) - (v < 0.0); }

void test_clear_orientation() {
    // Counter-clockwise triangle -> positive; reversed -> negative.
    check(orient2d(0, 0, 1, 0, 0, 1) > 0.0, "CCW triangle is positive");
    check(orient2d(0, 0, 0, 1, 1, 0) < 0.0, "CW triangle is negative");
    // Sign is invariant under cyclic permutation, negated under a swap.
    const double s = orient2d(0, 0, 4, 1, 1, 5);
    check(sign(orient2d(4, 1, 1, 5, 0, 0)) == sign(s), "cyclic permutation keeps sign");
    check(sign(orient2d(4, 1, 0, 0, 1, 5)) == -sign(s), "single swap negates sign");
}

void test_exact_collinear() {
    // Points on a line must return EXACTLY zero (integer coords, slope 1/2).
    check(orient2d(0, 0, 4, 2, 2, 1) == 0.0, "collinear -> exactly 0 (a)");
    check(orient2d(0, 0, 6, 3, 2, 1) == 0.0, "collinear -> exactly 0 (b)");
    check(orient2d(-3, -3, 5, 5, 100, 100) == 0.0, "collinear on y=x -> 0");
}

void test_near_degenerate() {
    // For a = (0,0), b = (1,1) the exact determinant equals (cy - cx). Place c one
    // representable step above / below the line y = x: the robust predicate must
    // report the correct sign of a single-ulp perturbation.
    const double base = 0.3;
    const double above = std::nextafter(base, 1.0);   // base + 1 ulp
    const double below = std::nextafter(base, -1.0);  // base - 1 ulp

    check(orient2d(0, 0, 1, 1, base, above) > 0.0,
          "one ulp above y=x -> strictly positive");
    check(orient2d(0, 0, 1, 1, base, below) < 0.0,
          "one ulp below y=x -> strictly negative");
    check(orient2d(0, 0, 1, 1, base, base) == 0.0, "exactly on y=x -> 0");

    // A skew, large-magnitude case where the two determinant halves nearly cancel.
    const double x = 12.34, y = 56.78;
    check(orient2d(x, y, x + 1.0, y + 1.0, x + 2.0, std::nextafter(y + 2.0, 1e9)) > 0.0,
          "near-collinear skew line: correct positive sign");
}

void test_incircle() {
    // CCW reference triangle; circumcircle centred at (0.5,0.5), radius sqrt(0.5).
    const double ax = 0, ay = 0, bx = 1, by = 0, cx = 0, cy = 1;
    check(incircle(ax, ay, bx, by, cx, cy, 0.5, 0.5) > 0.0,
          "circumcentre is inside the circumcircle");
    check(incircle(ax, ay, bx, by, cx, cy, 2.0, 2.0) < 0.0,
          "far point is outside the circumcircle");
    check(incircle(ax, ay, bx, by, cx, cy, 1.0, 1.0) == 0.0 ||
              std::fabs(incircle(ax, ay, bx, by, cx, cy, 1.0, 1.0)) < 1e-9,
          "point on the circumcircle gives ~0");

    // Four points exactly on the unit circle are cocircular.
    check(std::fabs(incircle(1, 0, 0, 1, -1, 0, 0, -1)) < 1e-9,
          "cocircular points -> ~0");

    // Near-cocircular: nudge d radially in/out by one representable step.
    const double r_in = std::nextafter(1.0, 0.0);   // just inside
    const double r_out = std::nextafter(1.0, 2.0);  // just outside
    check(incircle(1, 0, 0, 1, -1, 0, 0, -r_in) > 0.0,
          "one ulp inside the circle -> positive");
    check(incircle(1, 0, 0, 1, -1, 0, 0, -r_out) < 0.0,
          "one ulp outside the circle -> negative");
}

} // namespace

int main() {
    test_clear_orientation();
    test_exact_collinear();
    test_near_degenerate();
    test_incircle();
    if (g_failures == 0) {
        std::printf("OK: robust orientation predicate verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
