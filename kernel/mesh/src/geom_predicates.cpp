#include <katai/mesh/geom_predicates.hpp>

#include <cmath>

// Adaptive, exact-sign geometric predicates after Shewchuk (1997). The exact
// fallback uses nonoverlapping floating-point "expansions": a value is an
// increasing-magnitude sequence of doubles whose exact sum it represents. The
// error-free transforms below (TwoSum/TwoDiff/TwoProduct) are exact in IEEE-754
// double; TwoProduct relies on a fused multiply-add (std::fma), so the product
// error is captured without loss.

namespace katai::mesh {
namespace {

// Machine epsilon (half ulp at 1.0) and Shewchuk's orientation error bound.
constexpr double kEpsilon = 1.1102230246251565e-16;          // 2^-53
constexpr double kCcwErrboundA = (3.0 + 16.0 * kEpsilon) * kEpsilon;
constexpr double kIncErrboundA = (10.0 + 96.0 * kEpsilon) * kEpsilon;

// Error-free transforms (a op b == x + y exactly).
inline void two_sum(double a, double b, double& x, double& y) {
    x = a + b;
    const double bvirt = x - a;
    const double avirt = x - bvirt;
    const double bround = b - bvirt;
    const double around = a - avirt;
    y = around + bround;
}

inline void two_diff(double a, double b, double& x, double& y) {
    x = a - b;
    const double bvirt = a - x;
    const double avirt = x + bvirt;
    const double bround = bvirt - b;
    const double around = a - avirt;
    y = around + bround;
}

inline void two_product(double a, double b, double& x, double& y) {
    x = a * b;
    y = std::fma(a, b, -x);  // exact product error
}

// (a1 + a0) - b  ->  x2 + x1 + x0
inline void two_one_diff(double a1, double a0, double b, double& x2, double& x1,
                         double& x0) {
    double i;
    two_diff(a0, b, i, x0);
    two_sum(a1, i, x2, x1);
}

// (a1 + a0) - (b1 + b0)  ->  x3 + x2 + x1 + x0
inline void two_two_diff(double a1, double a0, double b1, double b0, double& x3,
                         double& x2, double& x1, double& x0) {
    double j, t;
    two_one_diff(a1, a0, b0, j, t, x0);
    two_one_diff(j, t, b1, x3, x2, x1);
}

// Sum of two nonoverlapping, increasing-magnitude expansions e and f, with zero
// components eliminated (Shewchuk, fast_expansion_sum_zeroelim). Returns length.
int expansion_sum(int elen, const double* e, int flen, const double* f,
                  double* h) {
    double q, qnew, hh;
    int eindex = 0, findex = 0, hindex = 0;
    double enow = e[0], fnow = f[0];

    if ((fnow > enow) == (fnow > -enow)) {
        q = enow;
        enow = (++eindex < elen) ? e[eindex] : 0.0;
    } else {
        q = fnow;
        fnow = (++findex < flen) ? f[findex] : 0.0;
    }
    while ((eindex < elen) && (findex < flen)) {
        if ((fnow > enow) == (fnow > -enow)) {
            two_sum(q, enow, qnew, hh);
            enow = (++eindex < elen) ? e[eindex] : 0.0;
        } else {
            two_sum(q, fnow, qnew, hh);
            fnow = (++findex < flen) ? f[findex] : 0.0;
        }
        q = qnew;
        if (hh != 0.0) h[hindex++] = hh;
    }
    while (eindex < elen) {
        two_sum(q, enow, qnew, hh);
        enow = (++eindex < elen) ? e[eindex] : 0.0;
        q = qnew;
        if (hh != 0.0) h[hindex++] = hh;
    }
    while (findex < flen) {
        two_sum(q, fnow, qnew, hh);
        fnow = (++findex < flen) ? f[findex] : 0.0;
        q = qnew;
        if (hh != 0.0) h[hindex++] = hh;
    }
    if ((q != 0.0) || (hindex == 0)) h[hindex++] = q;
    return hindex;
}

// Exact orientation determinant: the sign of the most significant component of
// the expansion equals the sign of the true determinant.
double orient2d_exact(double ax, double ay, double bx, double by, double cx,
                      double cy) {
    double axby1, axby0, axcy1, axcy0, bxcy1, bxcy0, bxay1, bxay0, cxay1, cxay0,
        cxby1, cxby0;
    double aterms[4], bterms[4], cterms[4];

    // det = ax*by - ax*cy + bx*cy - bx*ay + cx*ay - cx*by  (exact).
    two_product(ax, by, axby1, axby0);
    two_product(ax, cy, axcy1, axcy0);
    two_two_diff(axby1, axby0, axcy1, axcy0, aterms[3], aterms[2], aterms[1],
                 aterms[0]);

    two_product(bx, cy, bxcy1, bxcy0);
    two_product(bx, ay, bxay1, bxay0);
    two_two_diff(bxcy1, bxcy0, bxay1, bxay0, bterms[3], bterms[2], bterms[1],
                 bterms[0]);

    two_product(cx, ay, cxay1, cxay0);
    two_product(cx, by, cxby1, cxby0);
    two_two_diff(cxay1, cxay0, cxby1, cxby0, cterms[3], cterms[2], cterms[1],
                 cterms[0]);

    double ab[8];
    const int ablen = expansion_sum(4, aterms, 4, bterms, ab);
    double det[12];
    const int detlen = expansion_sum(ablen, ab, 4, cterms, det);
    return det[detlen - 1];  // most significant component -> exact sign
}

// --- Double-double arithmetic (~106-bit) for the in-circle exact fallback ----
// A DD value (hi, lo) represents the unevaluated sum hi + lo with |lo| <= ulp(hi)/2.
struct DD {
    double hi, lo;
};

inline DD dd_diff(double a, double b) {  // exact difference of two doubles
    DD r;
    two_diff(a, b, r.hi, r.lo);
    return r;
}

inline DD dd_add(DD a, DD b) {
    double s, e;
    two_sum(a.hi, b.hi, s, e);
    e += a.lo + b.lo;
    DD r;
    two_sum(s, e, r.hi, r.lo);
    return r;
}

inline DD dd_sub(DD a, DD b) { return dd_add(a, {-b.hi, -b.lo}); }

inline DD dd_mul(DD a, DD b) {
    double p, e;
    two_product(a.hi, b.hi, p, e);
    e += a.hi * b.lo + a.lo * b.hi;
    DD r;
    two_sum(p, e, r.hi, r.lo);
    return r;
}

// In-circle determinant recomputed in double-double; the differences are exact
// and every product/sum keeps ~106 bits, resolving near-cocircular signs.
double incircle_dd(double ax, double ay, double bx, double by, double cx,
                   double cy, double dx, double dy) {
    const DD adx = dd_diff(ax, dx), ady = dd_diff(ay, dy);
    const DD bdx = dd_diff(bx, dx), bdy = dd_diff(by, dy);
    const DD cdx = dd_diff(cx, dx), cdy = dd_diff(cy, dy);

    const DD alift = dd_add(dd_mul(adx, adx), dd_mul(ady, ady));
    const DD blift = dd_add(dd_mul(bdx, bdx), dd_mul(bdy, bdy));
    const DD clift = dd_add(dd_mul(cdx, cdx), dd_mul(cdy, cdy));

    const DD bc = dd_sub(dd_mul(bdx, cdy), dd_mul(cdx, bdy));
    const DD ca = dd_sub(dd_mul(cdx, ady), dd_mul(adx, cdy));
    const DD ab = dd_sub(dd_mul(adx, bdy), dd_mul(bdx, ady));

    const DD det = dd_add(dd_add(dd_mul(alift, bc), dd_mul(blift, ca)),
                          dd_mul(clift, ab));
    return det.hi + det.lo;  // sign carried faithfully
}

} // namespace

double orient2d(double ax, double ay, double bx, double by, double cx,
                double cy) {
    const double detleft = (ax - cx) * (by - cy);
    const double detright = (ay - cy) * (bx - cx);
    const double det = detleft - detright;

    // Static error filter: when the two halves have opposite signs (or one is
    // zero) there is no cancellation and the sign is certain immediately.
    double detsum;
    if (detleft > 0.0) {
        if (detright <= 0.0) return det;
        detsum = detleft + detright;
    } else if (detleft < 0.0) {
        if (detright >= 0.0) return det;
        detsum = -detleft - detright;
    } else {
        return det;
    }
    const double errbound = kCcwErrboundA * detsum;
    if (det >= errbound || -det >= errbound) return det;

    return orient2d_exact(ax, ay, bx, by, cx, cy);  // sign not certain -> exact
}

double incircle(double ax, double ay, double bx, double by, double cx,
                double cy, double dx, double dy) {
    const double adx = ax - dx, ady = ay - dy;
    const double bdx = bx - dx, bdy = by - dy;
    const double cdx = cx - dx, cdy = cy - dy;

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;

    const double alift = adx * adx + ady * ady;
    const double blift = bdx * bdx + bdy * bdy;
    const double clift = cdx * cdx + cdy * cdy;

    const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) +
                       clift * (adxbdy - bdxady);

    const double permanent =
        (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
        (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
        (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
    const double errbound = kIncErrboundA * permanent;
    if (det > errbound || -det > errbound) return det;

    return incircle_dd(ax, ay, bx, by, cx, cy, dx, dy);  // uncertain -> dd
}

} // namespace katai::mesh
