// A line between two refined corners is refined along its length. The automatic refinement at an angle
// that is not a multiple of 90 degrees multiplies the local element size of the geometry; read as a size
// at the corner point alone, it gave an embankment one fine patch at each toe and elements of the full
// size along the rest of its base, where the slip surfaces of its safety analysis run. The size now runs
// along the line, from the size at one end to the size at the other.
//
// verify: none -- a mesh sizing rule, held by the sizes it asks for on the line between two corners, by the far field it leaves alone and by the right-angled model it leaves bit-identical; what it does to the factors of safety of the embankments is measured in the changelog
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

m::Project embankment(double elem_size, bool corners) {
    m::Project pr;
    pr.name = "embankment";
    pr.x_min = 0; pr.x_max = 50; pr.y_min = -6; pr.y_max = 4;
    m::Material s; s.name = "Soil"; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3;
    pr.materials.push_back(s);
    const auto poly = [&](std::vector<double> x, std::vector<double> y) {
        m::SoilPolygon p; p.x = std::move(x); p.y = std::move(y); p.material = 0;
        pr.polygons.push_back(std::move(p));
    };
    poly({0, 50, 50, 36, 14, 0}, {-6, -6, 0, 0, 0, 0});
    poly({14, 36, 24, 22}, {0, 0, 4, 4});
    pr.mesh.elem_size = elem_size;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = true;
    pr.mesh.corner_refine = corners;
    return pr;
}

m::Project block(double elem_size, bool corners) {
    m::Project pr;
    pr.name = "block";
    pr.x_min = 0; pr.x_max = 40; pr.y_min = 0; pr.y_max = 20;
    m::Material s; s.name = "Soil"; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3;
    pr.materials.push_back(s);
    m::SoilPolygon p; p.x = {0, 40, 40, 0}; p.y = {0, 0, 20, 20}; p.material = 0;
    pr.polygons.push_back(std::move(p));
    pr.mesh.elem_size = elem_size;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = true;
    pr.mesh.corner_refine = corners;
    return pr;
}

// The mean length of the element edges lying on y = 0 between x0 and x1, and of those inside a window.
double edges_on_base(const katai::mesh::Mesh& mesh, double x0, double x1) {
    double sum = 0.0;
    int n = 0;
    for (int e = 0; e < mesh.element_count; ++e)
        for (int c = 0; c < 3; ++c) {
            const int a = mesh.node_of(e, c), b = mesh.node_of(e, (c + 1) % 3);
            if (std::fabs(mesh.y[a]) > 1e-9 || std::fabs(mesh.y[b]) > 1e-9) continue;
            if (std::min(mesh.x[a], mesh.x[b]) < x0 || std::max(mesh.x[a], mesh.x[b]) > x1) continue;
            sum += std::fabs(mesh.x[a] - mesh.x[b]);
            ++n;
        }
    return n ? sum / n : 0.0;
}
}  // namespace

int main() {
    std::printf("== a line between refined corners ==\n");
    constexpr double kSize = 2.0;
    const auto R = katai::app::mesh_from_project(embankment(kSize, true));
    const auto P = katai::app::mesh_from_project(embankment(kSize, false));
    check(R.ok && P.ok, "both mesh");
    if (!R.ok || !P.ok) return 1;
    // The corners of this embankment are 26.6 and 18.4 degrees: factors 0.42 and 0.30 on the size, so the middle
    // of the base asks for about 0.36 of it.
    const double mid_r = edges_on_base(R.mesh, 20.0, 30.0), mid_p = edges_on_base(P.mesh, 20.0, 30.0);
    std::printf("      elements %d against %d; edge at the middle of the base %.2f against %.2f m\n",
                R.mesh.element_count, P.mesh.element_count, mid_r, mid_p);
    check(mid_r < 0.55 * kSize, "the middle of the base is refined");
    check(mid_p > 0.8 * kSize, "and was not without the rule");
    check(mid_r > 0.2 * kSize, "but no finer than its two corners ask");
    // Off the model's own lines the size is what it was.
    const double far_r = edges_on_base(R.mesh, 0.0, 5.0), far_p = edges_on_base(P.mesh, 0.0, 5.0);
    std::printf("      edge far along the surface %.2f against %.2f m\n", far_r, far_p);
    check(far_r > 0.6 * kSize, "the far end of the surface keeps about its size");
    check(R.mesh.element_count < 4 * P.mesh.element_count, "and the mesh is not refined everywhere");
    // Right angles ask for nothing.
    const auto A = katai::app::mesh_from_project(block(kSize, true));
    const auto B = katai::app::mesh_from_project(block(kSize, false));
    check(A.ok && B.ok && A.mesh.node_count == B.mesh.node_count && A.mesh.element_count == B.mesh.element_count,
          "a model of right angles is the same mesh with the rule on");
    bool same = A.ok && B.ok && A.mesh.node_count == B.mesh.node_count;
    for (int n = 0; same && n < A.mesh.node_count; ++n)
        same = A.mesh.x[n] == B.mesh.x[n] && A.mesh.y[n] == B.mesh.y[n];
    check(same, "node for node");
    if (g_failures) {
        std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
        return 1;
    }
    std::printf("OK: a line between refined corners verified\n");
    return 0;
}
