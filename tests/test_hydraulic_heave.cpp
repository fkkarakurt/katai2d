// Hydraulic heave below exposed ground (katai/analysis/hydraulic_heave.hpp), checked against the
// one head field whose answer is closed form: steady upward flow at a uniform gradient.
//
//   problem:  a 10 m wide, 6 m deep ground (y = 0..6) under free water standing 1 m above its
//             surface, flowing upwards at a uniform gradient a: h(y) = 7 + a (6 - y). Then
//             u(z) = gamma_w (1 + (1 + a) z) at depth z, the hydrostatic part u_h = gamma_w (1 + z)
//             and the excess gamma_w a z, against a submerged weight (gamma_sat - gamma_w) z: the
//             factor is gamma' / (gamma_w a) at EVERY depth -- the critical gradient over the
//             gradient. A linear head is reproduced exactly by both element types and the unit
//             weight is constant, so nothing but round-off separates the check from it.
//   expected: (a) a = 0.3, gamma_sat = 20: F = 10/3; (b) a = 1.25: F = 0.8, heave; (c) a = 0:
//             no excess, nothing reported; (d) the left half dug 2 m deep, the dug-out soil
//             gamma_sat = 22: the smallest factor is in the excavation, F = 1/a at y_top = 4,
//             because the undug columns carry the heavier layer on top (their least factor,
//             at the foot of the column, is (2*12 + 4*10) / (6*10 a) = 1.07/a); (e) a dry
//             surface (water 1 m below it) is not a column of water to lift: nothing reported;
//             (f) the water at the surface itself, h(y) = 6 + a (6 - y) -- the seepage face of a
//             floor pumped dry -- the same gamma' / (gamma_w a).
//   band:     1e-9 relative -- round-off.
//
// verify: KV-FLW-005
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley -- failure by heave: the critical hydraulic gradient i_c = gamma'/gamma_w, at which upward seepage removes the vertical effective stress
//   locator:  steady upward flow at a uniform gradient a below free water standing 1 m on the surface, h(y) = 7 + a (6 - y): u(z) = gamma_w (1 + (1 + a) z), u_h(z) = gamma_w (1 + z), sigma_v(z) = gamma_w + gamma_sat z, so (sigma_v - u_h) / (u - u_h) = (gamma_sat - gamma_w) / (gamma_w a) at every depth z (stated in full); the same with the water at the surface (a seepage face), and a dug half whose floor carries none of the heavier layer above it
//   quantity: the least ratio of submerged weight to excess pore pressure below the exposed surface [-]
//   expected: 10/3 at a = 0.3 (gamma_sat = 20, gamma_w = 10), 0.8 at a = 1.25, 2 below the floor dug at a = 0.5; nothing reported for hydrostatic water or a dry surface
//   band:     1e-9 relative, as asserted below -- a linear head is reproduced exactly by 6- and 15-noded elements and the unit weight is constant along each column, so only round-off remains

#include <katai/analysis/hydraulic_heave.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kGw = 10.0;

Eigen::VectorXd linear_head(const katai::mesh::Mesh& m, double top, double a) {
    Eigen::VectorXd h(m.node_count);
    for (int n = 0; n < m.node_count; ++n) h(n) = top + a * (6.0 - m.y[n]);
    return h;
}

void run(const char* label, katai::mesh::Mesh (*gen)(const katai::geometry::RectangularDomain&, int, int)) {
    std::printf("\n%s\n", label);
    const katai::geometry::RectangularDomain dom{0.0, 0.0, 10.0, 6.0, 0};
    katai::mesh::Mesh m = gen(dom, 10, 6);
    const std::vector<double> g{18.0}, gs{20.0};
    const std::vector<char> all;

    // (a)
    auto r = katai::core::hydraulic_heave(m, all, g, gs, linear_head(m, 7.0, 0.3), kGw);
    std::printf("      a = 0.3:  F = %.12f over %d columns (closed form %.12f)\n", r.factor,
                r.columns, 10.0 / 3.0);
    check(r.found && std::fabs(r.factor - 10.0 / 3.0) <= 1e-9 * 10.0 / 3.0 && r.columns == 10,
          "uniform upward gradient: gamma'/(gamma_w i), every top column measured");
    // (b)
    r = katai::core::hydraulic_heave(m, all, g, gs, linear_head(m, 7.0, 1.25), kGw);
    std::printf("      a = 1.25: F = %.12f\n", r.factor);
    check(r.found && std::fabs(r.factor - 0.8) <= 1e-9, "past the critical gradient: F = 0.8 < 1");
    // (c)
    r = katai::core::hydraulic_heave(m, all, g, gs, linear_head(m, 7.0, 0.0), kGw);
    check(!r.found, "hydrostatic water: no excess, nothing reported");
    // (e)
    r = katai::core::hydraulic_heave(m, all, g, gs, linear_head(m, 5.0, 0.3), kGw);
    check(!r.found, "a dry surface: nothing reported");
    // (f)
    r = katai::core::hydraulic_heave(m, all, g, gs, linear_head(m, 6.0, 0.3), kGw);
    std::printf("      seepage face: F = %.12f\n", r.factor);
    check(r.found && std::fabs(r.factor - 10.0 / 3.0) <= 1e-9 * 10.0 / 3.0,
          "the water at the surface itself (a floor pumped dry): the same factor");

    // (d) the left half dug 2 m deep; the top 2 m are a heavier soil
    for (int e = 0; e < m.element_count; ++e) {
        double cy = 0.0;
        for (int k = 0; k < 3; ++k) cy += m.y[m.node_of(e, k)] / 3.0;
        if (cy > 4.0) m.element_material[e] = 1;
    }
    std::vector<char> act(m.element_count, 1);
    for (int e = 0; e < m.element_count; ++e) {
        double cx = 0.0, cy = 0.0;
        for (int k = 0; k < 3; ++k) { cx += m.x[m.node_of(e, k)] / 3.0; cy += m.y[m.node_of(e, k)] / 3.0; }
        if (cx < 5.0 && cy > 4.0) act[e] = 0;
    }
    const std::vector<double> g2{18.0, 20.0}, gs2{20.0, 22.0};
    r = katai::core::hydraulic_heave(m, act, g2, gs2, linear_head(m, 7.0, 0.5), kGw);
    std::printf("      dug:      F = %.12f at x = %.2f, y_top = %.2f, depth %.2f\n", r.factor, r.x,
                r.y_top, r.depth);
    check(r.found && std::fabs(r.factor - 2.0) <= 1e-9 * 2.0 && std::fabs(r.y_top - 4.0) < 1e-12 &&
              r.x < 5.0,
          "an excavation: the least factor is below the dug surface, in the soil left there");
}

}  // namespace

int main() {
    run("6-noded elements", katai::mesh::generate_structured_tri6);
    run("15-noded elements", katai::mesh::generate_structured_tri15);
    if (g_failures) { std::printf("\n%d check(s) failed\n", g_failures); return 1; }
    std::printf("\nall checks passed\n");
    return 0;
}
