// Generalised (layered + water-aware) K0 initial stress: compute_k0_initial_stress_layered builds
// sigma'_v at each Gauss point from the VERTICAL integral of the effective unit weight above it,
// and sigma'_h = K0 sigma'_v. This is the foundation for a PLAXIS-style K0 procedure on a real GUI
// model (multiple soil layers, a water table, per-material K0).
//
//   (A) Uniform single layer, flat surface: sigma'_v = -gamma'(H - y) reproduced exactly (the
//       midpoint rule is exact for a constant integrand), sigma'_h = K0 sigma'_v.
//   (B) Two horizontal layers (gamma1 on top, gamma2 below the interface at h1): the overburden
//       is piecewise linear; checked at every Gauss point (fine integration -> sub-kPa error).
#include <katai/analysis/initial_stress.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::GaussState;
using katai::core::K0LayeredOptions;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;
namespace tri6 = katai::core::tri6;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// Gauss-point y at element e, gauss index g (matches solver ordering e*kGaussCount + g).
double gp_xy(const Mesh& mesh, int e, int g, double& x) {
    const auto gauss = tri6::gauss_points();
    const auto N = tri6::shape_functions(gauss[g].xi, gauss[g].eta);
    double y = 0.0; x = 0.0;
    for (int k = 0; k < 6; ++k) {
        const int n = mesh.node_of(e, k);
        x += N(k) * mesh.x[n]; y += N(k) * mesh.y[n];
    }
    return y;
}

void test_uniform_layer() {
    constexpr double W = 6.0, H = 10.0, gamma = 18.0, K0 = 0.5;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 8);

    K0LayeredOptions opt;
    opt.eff_unit_weight = [&](double, double) { return gamma; };
    opt.ground_surface = [&](double) { return H; };
    opt.k0 = {K0};
    const auto st = katai::core::compute_k0_initial_stress_layered(mesh, opt);

    double max_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e)
        for (int g = 0; g < 3; ++g) {
            double x; const double y = gp_xy(mesh, e, g, x);
            const GaussState& s = st[(size_t)e * 3 + g];
            const double sv = -gamma * (H - y), sh = K0 * sv;
            max_err = std::fmax(max_err, std::fabs(s.stress(1) - sv));
            max_err = std::fmax(max_err, std::fabs(s.stress(0) - sh));
            max_err = std::fmax(max_err, std::fabs(s.stress_zz - sh));
        }
    std::printf("  uniform layer: max|err|=%.3e\n", max_err);
    check(max_err < 1e-9 * gamma * H, "uniform K0 field matches -gamma(H-y), K0 sigma_v");
}

void test_two_layers() {
    constexpr double W = 6.0, H = 10.0, h1 = 4.0;   // interface at y = h1
    constexpr double g1 = 16.0, g2 = 20.0, K0 = 0.45;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 10);

    K0LayeredOptions opt;
    opt.eff_unit_weight = [&](double, double y) { return y >= h1 ? g1 : g2; };
    opt.ground_surface = [&](double) { return H; };
    opt.k0 = {K0};
    opt.integration_steps = 800;
    const auto st = katai::core::compute_k0_initial_stress_layered(mesh, opt);

    double max_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e)
        for (int g = 0; g < 3; ++g) {
            double x; const double y = gp_xy(mesh, e, g, x);
            const double sv = y >= h1 ? -g1 * (H - y)
                                      : -(g1 * (H - h1) + g2 * (h1 - y));
            max_err = std::fmax(max_err, std::fabs(st[(size_t)e * 3 + g].stress(1) - sv));
        }
    std::printf("  two layers: max|sigma'_v err|=%.3e\n", max_err);
    check(max_err < 0.1, "layered sigma'_v matches piecewise overburden (< 0.1 kPa)");
}

} // namespace

int main() {
    test_uniform_layer();
    test_two_layers();
    if (g_failures == 0) {
        std::printf("OK: layered K0 initial stress verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
