// Verification of the incremental Delaunay triangulation (P1.4b).
//
// Checks the three properties that define a correct triangulation: every element
// is positively oriented, the elements exactly tile the convex hull (areas sum to
// the hull area -- no gaps or overlaps), and the Delaunay empty-circumcircle
// property holds (no vertex lies inside any triangle's circumcircle). Robustness
// is exercised on a cocircular grid and a large random cloud.
#include <katai/mesh/delaunay.hpp>
#include <katai/mesh/geom_predicates.hpp>

#include <cstdint>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

using katai::mesh::constrained_delaunay;
using katai::mesh::delaunay_triangulate;
using katai::mesh::incircle;
using katai::mesh::orient2d;
using katai::mesh::quality_mesh;
using katai::mesh::Triangulation;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// Twice the signed area of triangle (a,b,c).
double area2(const std::vector<double>& x, const std::vector<double>& y,
             const std::array<int, 3>& t) {
    return orient2d(x[t[0]], y[t[0]], x[t[1]], y[t[1]], x[t[2]], y[t[2]]);
}

// Twice the area of the convex hull (Andrew's monotone chain, via orient2d).
double hull_area2(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = static_cast<int>(x.size());
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return x[a] < x[b] || (x[a] == x[b] && y[a] < y[b]);
    });
    std::vector<int> h(2 * n);
    int k = 0;
    auto cross = [&](int o, int a, int b) {
        return orient2d(x[o], y[o], x[a], y[a], x[b], y[b]);
    };
    for (int i = 0; i < n; ++i) {  // lower hull
        while (k >= 2 && cross(h[k - 2], h[k - 1], idx[i]) <= 0) --k;
        h[k++] = idx[i];
    }
    const int lower = k + 1;
    for (int i = n - 2; i >= 0; --i) {  // upper hull
        while (k >= lower && cross(h[k - 2], h[k - 1], idx[i]) <= 0) --k;
        h[k++] = idx[i];
    }
    h.resize(k - 1);  // drop the repeated first vertex
    double a2 = 0.0;
    for (std::size_t i = 1; i + 1 < h.size(); ++i)
        a2 += orient2d(x[h[0]], y[h[0]], x[h[i]], y[h[i]], x[h[i + 1]], y[h[i + 1]]);
    return a2;
}

bool all_ccw(const Triangulation& tg, const std::vector<double>& x,
             const std::vector<double>& y) {
    for (const auto& t : tg.triangles)
        if (area2(x, y, t) <= 0.0) return false;
    return true;
}

double total_area2(const Triangulation& tg, const std::vector<double>& x,
                   const std::vector<double>& y) {
    double s = 0.0;
    for (const auto& t : tg.triangles) s += area2(x, y, t);
    return s;
}

// Deterministic pseudo-random points in the unit square.
void random_cloud(int n, std::vector<double>& x, std::vector<double>& y) {
    std::uint64_t s = 0x9e3779b97f4a7c15ULL;
    auto next = [&] {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return (s >> 11) * (1.0 / 9007199254740992.0);  // [0,1)
    };
    x.resize(n); y.resize(n);
    for (int i = 0; i < n; ++i) { x[i] = next(); y[i] = next(); }
}

void test_square() {
    const std::vector<double> x = {0, 1, 1, 0}, y = {0, 0, 1, 1};
    const Triangulation tg = delaunay_triangulate(x, y);
    check(tg.triangles.size() == 2, "unit square -> 2 triangles");
    check(all_ccw(tg, x, y), "square: triangles CCW");
    check(std::fabs(total_area2(tg, x, y) - 2.0) < 1e-12, "square: area = 1");
}

void test_random_delaunay() {
    std::vector<double> x, y;
    random_cloud(400, x, y);
    const Triangulation tg = delaunay_triangulate(x, y);
    check(!tg.triangles.empty(), "random: produced triangles");
    check(all_ccw(tg, x, y), "random: all triangles CCW");
    check(std::fabs(total_area2(tg, x, y) - hull_area2(x, y)) < 1e-9,
          "random: triangles tile the convex hull");

    // Delaunay empty-circumcircle property.
    bool empty = true;
    const int nt = static_cast<int>(tg.triangles.size());
    for (int ti = 0; ti < nt && empty; ++ti) {
        const auto& t = tg.triangles[ti];
        for (int q = 0; q < tg.point_count; ++q) {
            if (q == t[0] || q == t[1] || q == t[2]) continue;
            if (incircle(x[t[0]], y[t[0]], x[t[1]], y[t[1]], x[t[2]], y[t[2]],
                         x[q], y[q]) > 1e-9) {
                empty = false;
                break;
            }
        }
    }
    check(empty, "random: empty-circumcircle (Delaunay) property holds");
}

void test_degenerate_grid() {
    // A regular grid is maximally cocircular -- a stress test for the predicates.
    std::vector<double> x, y;
    const int m = 8;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i) { x.push_back(i); y.push_back(j); }
    const Triangulation tg = delaunay_triangulate(x, y);
    check(all_ccw(tg, x, y), "grid: all triangles CCW");
    const double expected = 2.0 * (m - 1) * (m - 1);  // 2 * hull area
    check(std::fabs(total_area2(tg, x, y) - expected) < 1e-9,
          "grid: tiles the square (no gaps/overlaps)");
}

void test_large_robustness() {
    std::vector<double> x, y;
    random_cloud(3000, x, y);
    const Triangulation tg = delaunay_triangulate(x, y);
    check(all_ccw(tg, x, y), "3000 pts: all triangles CCW");
    check(std::fabs(total_area2(tg, x, y) - hull_area2(x, y)) < 1e-8,
          "3000 pts: tiles the convex hull");
}

bool has_edge(const Triangulation& tg, int a, int b) {
    for (const auto& t : tg.triangles) {
        bool ha = false, hb = false;
        for (int k = 0; k < 3; ++k) { ha |= t[k] == a; hb |= t[k] == b; }
        if (ha && hb) return true;
    }
    return false;
}

void test_constrained_diagonal() {
    // A general convex quad has one Delaunay diagonal; constraining the OTHER one
    // must force it into the triangulation via edge flipping.
    const std::vector<double> x = {0, 4, 5, 1}, y = {0, 0, 3, 2};
    const Triangulation un = delaunay_triangulate(x, y);
    const int da = has_edge(un, 0, 2) ? 1 : 0;  // absent diagonal endpoints
    const int db = has_edge(un, 0, 2) ? 3 : 2;
    check(!has_edge(un, da, db), "the other diagonal is initially absent");

    const Triangulation cd = constrained_delaunay(x, y, {{{da, db}}});
    check(has_edge(cd, da, db), "constrained diagonal is recovered");
    check(all_ccw(cd, x, y), "constrained quad: triangles CCW");
    check(std::fabs(total_area2(cd, x, y) - hull_area2(x, y)) < 1e-12,
          "constrained quad: still tiles the hull");
}

void test_constrained_long_edge() {
    // A 6x6 grid with a long constraint that crosses many triangles. The chosen
    // endpoints (0,0)-(5,4) have no grid vertex in the segment interior.
    std::vector<double> x, y;
    const int m = 6;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i) { x.push_back(i); y.push_back(j); }
    auto idx = [&](int i, int j) { return j * m + i; };
    const int a = idx(0, 0), b = idx(5, 4);

    const Triangulation cd = constrained_delaunay(x, y, {{{a, b}}});
    check(has_edge(cd, a, b), "long constraint crossing many triangles recovered");
    check(all_ccw(cd, x, y), "constrained grid: triangles CCW");
    check(std::fabs(total_area2(cd, x, y) - 2.0 * (m - 1) * (m - 1)) < 1e-9,
          "constrained grid: tiles the square (no gaps/overlaps)");
}

double tri_min_angle_deg(const Triangulation& tg, const std::array<int, 3>& t) {
    const double pi = 3.14159265358979323846;
    auto angle_at = [&](int a, int b, int c) {
        const double abx = tg.x[b] - tg.x[a], aby = tg.y[b] - tg.y[a];
        const double acx = tg.x[c] - tg.x[a], acy = tg.y[c] - tg.y[a];
        const double dot = abx * acx + aby * acy;
        const double mag = std::sqrt((abx * abx + aby * aby) * (acx * acx + acy * acy));
        double cosv = dot / mag;
        cosv = std::max(-1.0, std::min(1.0, cosv));
        return std::acos(cosv) * 180.0 / pi;
    };
    return std::min({angle_at(t[0], t[1], t[2]), angle_at(t[1], t[2], t[0]),
                     angle_at(t[2], t[0], t[1])});
}

void test_quality_refinement() {
    // Square domain with its four boundary edges as constraints; refine to a
    // minimum angle and a maximum element area.
    const std::vector<double> x = {0, 4, 4, 0}, y = {0, 0, 4, 4};
    const std::vector<std::array<int, 2>> segs = {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}};
    const double min_angle = 20.0, max_area = 0.25;
    const Triangulation tg = quality_mesh(x, y, segs, min_angle, max_area);

    check(tg.triangles.size() >= 64, "refinement produced a fine mesh");
    check(tg.x.size() > 4, "Steiner points were added");

    double worst = 180.0, area2_sum = 0.0;
    bool ccw = true, within_area = true;
    for (const auto& t : tg.triangles) {
        const double o = orient2d(tg.x[t[0]], tg.y[t[0]], tg.x[t[1]], tg.y[t[1]],
                                  tg.x[t[2]], tg.y[t[2]]);
        if (o <= 0.0) ccw = false;
        area2_sum += o;
        if (0.5 * o > max_area + 1e-9) within_area = false;
        worst = std::min(worst, tri_min_angle_deg(tg, t));
    }
    check(ccw, "refined: all triangles CCW");
    check(within_area, "refined: every triangle within max_area");
    check(std::fabs(area2_sum - 2.0 * 16.0) < 1e-7, "refined: tiles the square");
    std::printf("  [refine] triangles=%zu vertices=%zu min_angle=%.2f deg\n",
                tg.triangles.size(), tg.x.size(), worst);
    check(worst >= min_angle - 0.5, "refined: minimum angle bound achieved");
}

// Non-convex (L-shaped) domain: the re-entrant corner makes the convex hull
// strictly larger than the polygon. The quality mesher must tile the POLYGON, not
// the hull -- triangles filling the concave notch (inside the hull, outside the
// domain) must be neither produced nor emitted. This guards the slope-toe class
// of bug where phantom exterior elements break the downstream FE solve.
void test_nonconvex_domain() {
    // CCW L-shape; re-entrant vertex at (2,2). Area = 4*4 - 2*2 = 12 (hull is 14).
    const std::vector<double> x = {0, 4, 4, 2, 2, 0}, y = {0, 0, 2, 2, 4, 4};
    const std::vector<std::array<int, 2>> segs = {
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 4}}, {{4, 5}}, {{5, 0}}};
    const Triangulation tg = quality_mesh(x, y, segs, 20.0, 0.5);

    bool ccw = true;
    double area2_sum = 0.0;
    for (const auto& t : tg.triangles) {
        const double o = orient2d(tg.x[t[0]], tg.y[t[0]], tg.x[t[1]], tg.y[t[1]],
                                  tg.x[t[2]], tg.y[t[2]]);
        if (o <= 0.0) ccw = false;
        area2_sum += o;
    }
    check(ccw, "L-shape: all triangles CCW");
    // Tiles the polygon (area 12), not the convex hull (area 14).
    check(std::fabs(0.5 * area2_sum - 12.0) < 1e-7,
          "L-shape: mesh tiles the polygon, not the convex hull");
    // No vertex escapes the polygon bounding box (no runaway circumcentres).
    bool in_box = true;
    for (std::size_t i = 0; i < tg.x.size(); ++i)
        if (tg.x[i] < -1e-9 || tg.x[i] > 4 + 1e-9 ||
            tg.y[i] < -1e-9 || tg.y[i] > 4 + 1e-9)
            in_box = false;
    check(in_box, "L-shape: no vertex outside the domain bounding box");
    std::printf("  [non-convex] triangles=%zu vertices=%zu area=%.4f (poly=12)\n",
                tg.triangles.size(), tg.x.size(), 0.5 * area2_sum);
}

// Acute-corner domain: a thin wedge whose apex angle (~21.8 deg) is far below the
// 20 deg quality bound's comfortable regime. Plain midpoint subsegment splitting
// fails to terminate here (cascading mutual encroachment between the two legs);
// concentric-shell splitting must make it terminate, tile the wedge, and achieve
// the angle bound everywhere EXCEPT the elements incident to the sharp apex (no
// algorithm can beat the input angle there).
void test_acute_corner_domain() {
    // Right-triangle wedge; apex angle at vertex 0 = atan(4/10) ~= 21.8 deg.
    const std::vector<double> x = {0, 10, 10}, y = {0, 0, 4};
    const std::vector<std::array<int, 2>> segs = {{{0, 1}}, {{1, 2}}, {{2, 0}}};
    const Triangulation tg = quality_mesh(x, y, segs, 20.0, 1.0);

    check(!tg.triangles.empty(), "acute: refinement terminated with a mesh");
    bool ccw = true;
    double area2_sum = 0.0, worst_off_apex = 180.0;
    for (const auto& t : tg.triangles) {
        const double o = orient2d(tg.x[t[0]], tg.y[t[0]], tg.x[t[1]], tg.y[t[1]],
                                  tg.x[t[2]], tg.y[t[2]]);
        if (o <= 0.0) ccw = false;
        area2_sum += o;
        // The apex is input vertex 0 (extract preserves input indices).
        if (t[0] != 0 && t[1] != 0 && t[2] != 0)
            worst_off_apex = std::min(worst_off_apex, tri_min_angle_deg(tg, t));
    }
    check(ccw, "acute: all triangles CCW");
    check(std::fabs(0.5 * area2_sum - 20.0) < 1e-7, "acute: tiles the wedge (area 20)");
    std::printf("  [acute] triangles=%zu vertices=%zu min_angle(off-apex)=%.2f deg\n",
                tg.triangles.size(), tg.x.size(), worst_off_apex);
    check(worst_off_apex >= 20.0 - 0.5,
          "acute: angle bound achieved away from the sharp apex");
}

} // namespace

int main() {
    test_square();
    test_random_delaunay();
    test_degenerate_grid();
    test_large_robustness();
    test_constrained_diagonal();
    test_constrained_long_edge();
    test_quality_refinement();
    test_nonconvex_domain();
    test_acute_corner_domain();
    if (g_failures == 0) {
        std::printf("OK: incremental Delaunay triangulation verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
