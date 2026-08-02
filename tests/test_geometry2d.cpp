// Planar arrangement face extraction (GUI geometry -> soil clusters). Verifies that a soup of line
// segments yields the correct closed regions: a box = 1 face; a box split by a line = 2 faces; etc.
#include <katai/geometry/geometry2d.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::geometry::V2;
using katai::geometry::arrangement_faces;
using Seg = std::array<V2, 2>;

namespace {
int g_fail = 0;
void check(bool ok, const char* what) { if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_fail; } }

double area(const std::vector<V2>& p) {
    double a = 0; for (size_t k = 0; k < p.size(); ++k) { const V2& u = p[k]; const V2& v = p[(k + 1) % p.size()]; a += u.x * v.y - v.x * u.y; }
    return 0.5 * std::fabs(a);
}
std::vector<Seg> box(double x0, double y0, double x1, double y1) {
    return {Seg{V2{x0, y0}, V2{x1, y0}}, Seg{V2{x1, y0}, V2{x1, y1}},
            Seg{V2{x1, y1}, V2{x0, y1}}, Seg{V2{x0, y1}, V2{x0, y0}}};
}
double total_area(const std::vector<std::vector<V2>>& f) { double s = 0; for (auto& p : f) s += area(p); return s; }
}  // namespace

int main() {
    {  // a single box -> exactly one face of the box area
        auto f = arrangement_faces(box(0, 0, 40, 20));
        std::printf("  box: %zu face(s), area=%.1f (exp 800)\n", f.size(), total_area(f));
        check(f.size() == 1, "box -> 1 face");
        check(std::fabs(total_area(f) - 800.0) < 1e-6, "box face area = 800");
    }
    {  // box split by a horizontal mid-line -> two faces
        auto s = box(0, 0, 40, 20);
        s.push_back(Seg{V2{0, 10}, V2{40, 10}});
        auto f = arrangement_faces(s);
        std::printf("  box + mid line: %zu face(s), area=%.1f (exp 800)\n", f.size(), total_area(f));
        check(f.size() == 2, "box + horizontal line -> 2 faces");
        check(std::fabs(total_area(f) - 800.0) < 1e-6, "two faces cover the box");
    }
    {  // box split by a diagonal -> two triangles
        auto s = box(0, 0, 10, 10);
        s.push_back(Seg{V2{0, 0}, V2{10, 10}});
        auto f = arrangement_faces(s);
        std::printf("  box + diagonal: %zu face(s), area=%.1f (exp 100)\n", f.size(), total_area(f));
        check(f.size() == 2, "box + diagonal -> 2 faces");
        check(std::fabs(total_area(f) - 100.0) < 1e-6, "two triangles cover the box");
    }
    {  // a crossing inside a box -> four faces
        auto s = box(0, 0, 10, 10);
        s.push_back(Seg{V2{0, 5}, V2{10, 5}});
        s.push_back(Seg{V2{5, 0}, V2{5, 10}});
        auto f = arrangement_faces(s);
        std::printf("  box + cross: %zu face(s), area=%.1f (exp 100)\n", f.size(), total_area(f));
        check(f.size() == 4, "box + cross -> 4 faces");
        check(std::fabs(total_area(f) - 100.0) < 1e-6, "four faces cover the box");
    }
    if (g_fail == 0) { std::printf("OK: planar arrangement faces verified\n"); return 0; }
    std::fprintf(stderr, "%d check(s) failed\n", g_fail);
    return 1;
}
