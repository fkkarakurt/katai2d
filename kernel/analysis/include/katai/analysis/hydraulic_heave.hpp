#pragma once
// HYDRAULIC HEAVE below exposed ground, read from a steady seepage head field.
//
// Water flowing upwards through the ground loses head, and the head it loses is a force on the
// soil skeleton pointing up. Where that force reaches the submerged weight of the soil above a
// point, the vertical effective stress there is gone and the ground lifts: hydraulic heave, the
// failure an excavation dug below the water table next to a cut-off wall is checked for. The
// deformation analysis does not report it by itself -- a soil with cohesion can carry a zero or
// tensile vertical effective stress over a thin layer and converge -- so it is checked here,
// column by column, in the stress form of the hydraulic-heave limit state: at every depth z below
// the exposed surface the total vertical stress must exceed the pore pressure,
//
//     sigma_v(z) = u_0 + integral_0^z gamma dz'   >=   u(z),
//
// with u_0 the free water standing on the surface (it loads both sides) and gamma the saturated
// unit weight where the ground is wet, the unsaturated one where it is not. Neither the friction
// on the sides of the column nor any strength of the soil is counted, as in the classical check.
//
// The ratio reported is the one that separates the water from the soil:
//
//     F(z) = (sigma_v(z) - u_h(z)) / (u(z) - u_h(z)),     u_h(z) = u_0 + gamma_w z,
//
// the submerged weight of the column over the water pressure in excess of hydrostatic, which is
// gamma' / (gamma_w i_avg) -- the critical gradient over the gradient averaged over the depth z.
// F >= 1 is exactly sigma_v >= u, so F < 1 is heave whichever form is read; F is defined where
// there is an excess (upward flow) and the water stands at or above the top of the column. The smallest F over every
// depth of every column is what is reported, with where it was found: near the surface it is the
// local exit gradient, deeper it is the gradient averaged over the whole column above.
//
// Columns start at the midpoint of every exposed edge of the active ground that faces upwards
// (outward normal within 45 degrees of vertical) and walk down through active ground in steps of
// a quarter of that edge's length, until they leave the active ground or the mesh.

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <katai/fem/elements/element_traits.hpp>
#include <katai/fem/elements/point_location.hpp>
#include <katai/mesh/boundary_extraction.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

struct HeaveCheck {
    bool found = false;     // some column has upward flow below a wet exposed surface
    double factor = 0.0;    // the smallest F over every depth of every such column
    double x = 0.0;         // the column it was found in
    double y_top = 0.0;     // ... the exposed surface it starts from
    double depth = 0.0;     // ... and the depth below that surface
    int columns = 0;        // columns that had an excess to measure
};

namespace detail_heave {
template <class E>
HeaveCheck check(const mesh::Mesh& mesh, const std::vector<char>& active,
                 const std::vector<double>& gamma, const std::vector<double>& gamma_sat,
                 const Eigen::VectorXd& head, double gamma_w) {
    constexpr int N = E::kNodeCount;
    HeaveCheck out;
    const auto on = [&](int e) { return active.empty() || active[e]; };
    const auto coords = [&](int e) {
        typename E::NodeCoords X;
        for (int k = 0; k < N; ++k) {
            X(k, 0) = mesh.x[mesh.node_of(e, k)];
            X(k, 1) = mesh.y[mesh.node_of(e, k)];
        }
        return X;
    };
    // The element holding (px, py): the one the column was last in first, which is where a
    // column walking down usually still is.
    int last = -1;
    const auto locate = [&](double px, double py, double& xi, double& eta) {
        if (last >= 0) {
            const auto lc = ploc::physical_to_local<E>(coords(last), px, py);
            if (lc.converged && ploc::inside_ref_triangle(lc.xi, lc.eta, 1e-7)) {
                xi = lc.xi; eta = lc.eta;
                return last;
            }
        }
        // Only the active ground is searched: a point on the exposed surface also lies, to the
        // tolerance, in the excavated element above it.
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!on(e)) continue;
            double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
            for (int k = 0; k < 3; ++k) {
                const int n = mesh.node_of(e, k);
                xmin = std::min(xmin, mesh.x[n]); xmax = std::max(xmax, mesh.x[n]);
                ymin = std::min(ymin, mesh.y[n]); ymax = std::max(ymax, mesh.y[n]);
            }
            const double pad = 0.05 * (xmax - xmin + ymax - ymin);
            if (px < xmin - pad || px > xmax + pad || py < ymin - pad || py > ymax + pad) continue;
            const auto lc = ploc::physical_to_local<E>(coords(e), px, py);
            if (lc.converged && ploc::inside_ref_triangle(lc.xi, lc.eta, 1e-7)) {
                xi = lc.xi; eta = lc.eta;
                return last = e;
            }
        }
        return -1;
    };
    struct Sample { double h = 0.0, u = 0.0, gamma = 0.0; bool ok = false; };
    const auto sample = [&](double px, double py) {
        Sample s;
        double xi = 0.0, eta = 0.0;
        const int e = locate(px, py, xi, eta);
        if (e < 0) return s;
        const auto Nv = E::shape_functions(xi, eta);
        double h = 0.0;
        for (int k = 0; k < N; ++k) h += Nv(k) * head(mesh.node_of(e, k));
        s.h = h;
        s.u = gamma_w * std::max(0.0, h - py);
        const int m = mesh.element_material[e];
        if (m < 0 || m >= (int)gamma.size()) return s;
        s.gamma = s.u > 0.0 && m < (int)gamma_sat.size() ? gamma_sat[m] : gamma[m];
        s.ok = true;
        return s;
    };
    for (const auto& ed : mesh::extract_exposed_edges(mesh, active)) {
        if (ed.ny < std::sqrt(0.5)) continue;   // not an upward-facing surface
        const int a = ed.node[0], b = ed.node[ed.npe - 1];
        const double len = std::hypot(mesh.x[b] - mesh.x[a], mesh.y[b] - mesh.y[a]);
        if (len <= 0.0) continue;
        const double x0 = 0.5 * (mesh.x[a] + mesh.x[b]), y0 = 0.5 * (mesh.y[a] + mesh.y[b]);
        const double dz = 0.25 * len;
        last = -1;
        // The free water on the surface: the pore pressure just inside it, where the column's
        // depths are counted from.
        const double ys = y0 - 1e-6 * len;
        const Sample s0 = sample(x0, ys);
        // A surface above the water is not a column of water to lift; one at it -- the seepage
        // face of an excavation pumped dry to its floor -- is, and is the usual case.
        if (!s0.ok || s0.h < ys - 1e-3 * len) continue;
        const double u0 = s0.u;
        double sigma = u0, z = 0.0, g_prev = s0.gamma;
        bool counted = false;
        for (int k = 1; k <= 4000; ++k) {
            z = k * dz;
            const Sample s = sample(x0, ys - z);
            if (!s.ok) break;                     // out of the active ground, or of the mesh
            sigma += 0.5 * (g_prev + s.gamma) * dz;
            g_prev = s.gamma;
            const double u_h = u0 + gamma_w * z;
            const double excess = s.u - u_h;
            if (excess <= 1e-9 * std::max(1.0, s.u)) continue;   // no upward flow here
            const double F = (sigma - u_h) / excess;
            if (!counted) { ++out.columns; counted = true; }
            if (!out.found || F < out.factor) {
                out.found = true;
                out.factor = F;
                out.x = x0;
                out.y_top = y0;
                out.depth = z;
            }
        }
    }
    return out;
}
}  // namespace detail_heave

inline HeaveCheck hydraulic_heave(const mesh::Mesh& mesh, const std::vector<char>& active,
                                  const std::vector<double>& gamma,
                                  const std::vector<double>& gamma_sat,
                                  const Eigen::VectorXd& head, double gamma_w) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail_heave::check<Tri15Element>(mesh, active, gamma, gamma_sat, head, gamma_w);
    return detail_heave::check<Tri6Element>(mesh, active, gamma, gamma_sat, head, gamma_w);
}

}  // namespace katai::core
