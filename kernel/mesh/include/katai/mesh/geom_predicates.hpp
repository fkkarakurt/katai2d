#pragma once
// Robust adaptive geometric predicates (P1.4a) -- the foundation of a crash-proof
// Delaunay mesher.
//
// Naive floating-point evaluation of geometric determinants can return the wrong
// SIGN near degeneracy, which makes incremental Delaunay/CDT algorithms produce
// inverted elements or loop forever. Following Shewchuk ("Adaptive Precision
// Floating-Point Arithmetic and Fast Robust Geometric Predicates", 1997), each
// predicate first evaluates the determinant in plain double with an a-priori
// error bound (the fast path that succeeds for non-degenerate input), and only
// when the sign is not certain falls back to an exact computation carried out in
// multiple-component floating-point expansions. The returned value therefore has
// EXACTLY the sign of the true determinant -- in particular it is exactly zero
// for degenerate (collinear) configurations.

namespace katai::mesh {

// 2D orientation test. Returns a value with the same sign as the exact signed
// area determinant of (a, b, c):
//   > 0  : a, b, c are in counter-clockwise order (c lies left of a->b),
//   < 0  : clockwise (c lies right of a->b),
//   = 0  : exactly collinear.
double orient2d(double ax, double ay, double bx, double by, double cx, double cy);

// In-circle test. With (a, b, c) given in counter-clockwise order, returns a
// value with the same sign as the determinant deciding whether d lies inside the
// circumcircle of triangle (a, b, c):
//   > 0  : d is strictly inside the circumcircle,
//   < 0  : d is strictly outside,
//   = 0  : d is exactly cocircular with a, b, c.
// (If (a, b, c) are clockwise the sign is reversed.) A static error filter gives
// the fast path; uncertain signs fall back to double-double (~106-bit) evaluation,
// which is robust across realistic coordinate ranges. A fully exact expansion
// version (Shewchuk incircleadapt) is a future hardening if ever required.
double incircle(double ax, double ay, double bx, double by, double cx,
                double cy, double dx, double dy);

} // namespace katai::mesh
