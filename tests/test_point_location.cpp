// Point location + inverse isoparametric mapping (embedded beam foundation, Faz A.4). Verifies
// physical->local Newton mapping and element location: (a) round-trip phys->local->phys to
// round-off; (b) shape-function partition of unity; (c) exact interpolation of a linear field;
// (d) a DISTORTED (curved-edge) tri6 element so the Newton iteration is genuinely exercised.
#include <katai/fem/elements/point_location.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::Tri6Element;
using katai::core::Tri15Element;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;
namespace ploc = katai::core::ploc;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

template <class E>
void test_mesh(const Mesh& mesh, const char* tag) {
    constexpr int N = E::kNodeCount;
    const double a = 0.37, b = -0.21, c = 1.5;  // linear field f = a x + b y + c
    double max_rt = 0.0, max_pu = 0.0, max_lin = 0.0;
    int located = 0;
    const double pts[][2] = {{2.3, 1.7}, {5.5, 3.1}, {7.9, 0.6}, {1.1, 3.8}, {6.0, 2.0}};
    for (auto& q : pts) {
        const auto loc = ploc::locate_point(mesh, q[0], q[1]);
        if (!loc.found) continue;
        ++located;
        typename E::NodeCoords X;
        for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(loc.element, k)]; X(k, 1) = mesh.y[mesh.node_of(loc.element, k)]; }
        const auto sh = E::shape_functions(loc.xi, loc.eta);
        const double xi_ = X.col(0).dot(sh), yi_ = X.col(1).dot(sh);
        max_rt = std::fmax(max_rt, std::hypot(xi_ - q[0], yi_ - q[1]));   // round-trip
        max_pu = std::fmax(max_pu, std::fabs(sh.sum() - 1.0));            // partition of unity
        double f = 0.0;
        for (int k = 0; k < N; ++k) f += sh(k) * (a * X(k, 0) + b * X(k, 1) + c);
        max_lin = std::fmax(max_lin, std::fabs(f - (a * q[0] + b * q[1] + c)));  // linear exactness
    }
    std::printf("  %s: located %d/5  round-trip=%.2e  partition=%.2e  linear=%.2e\n",
                tag, located, max_rt, max_pu, max_lin);
    check(located == 5, "all interior points located");
    check(max_rt < 1e-8, "round-trip phys->local->phys to round-off");
    check(max_pu < 1e-12, "shape functions sum to 1");
    check(max_lin < 1e-9, "linear field interpolated exactly");
}

void test_distorted_tri6() {
    // Single tri6 with a midside node moved off-centre -> curved edge -> non-affine map.
    Tri6Element::NodeCoords X;
    X << 0.0, 0.0,   // c1
         4.0, 0.0,   // c2
         0.0, 3.0,   // c3
         2.3, 0.4,   // mid 1-2 (off the straight midpoint 2,0)
         2.4, 1.8,   // mid 2-3 (off 2,1.5)
         -0.3, 1.5;  // mid 3-1 (off 0,1.5)
    // Pick a local (xi,eta), map to physical, then invert -> must recover (xi,eta).
    const double targets[][2] = {{0.25, 0.25}, {0.1, 0.6}, {0.6, 0.1}, {0.34, 0.33}};
    double max_err = 0.0;
    for (auto& t : targets) {
        const auto sh = katai::core::tri6::shape_functions(t[0], t[1]);
        const double px = X.col(0).dot(sh), py = X.col(1).dot(sh);
        const auto lc = ploc::physical_to_local<Tri6Element>(X, px, py);
        check(lc.converged, "distorted-element inverse map converged");
        max_err = std::fmax(max_err, std::hypot(lc.xi - t[0], lc.eta - t[1]));
    }
    std::printf("  distorted (curved) tri6: max |local error| = %.2e\n", max_err);
    check(max_err < 1e-10, "inverse map exact on a curved element (Newton)");
}

} // namespace

int main() {
    const RectangularDomain domain{0.0, 0.0, 10.0, 4.0, 0};
    const Mesh m6 = katai::mesh::generate_structured_tri6(domain, 10, 4);
    const Mesh m15 = katai::mesh::generate_structured_tri15(domain, 5, 2);
    test_mesh<Tri6Element>(m6, "tri6");
    test_mesh<Tri15Element>(m15, "tri15");
    test_distorted_tri6();
    if (g_failures == 0) {
        std::printf("OK: point location + inverse isoparametric mapping verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
