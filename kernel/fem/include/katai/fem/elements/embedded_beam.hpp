#pragma once
// Embedded beam (pile row) — the skin (shaft) interaction (PLAXIS 2D Sci.Man §7.5, Phase
// A.4). The beam crosses the soil mesh in an ARBITRARY orientation; the skin interaction
// ties the beam nodes to the soil AT the beam's location (interpolated with N_s, a
// mesh-nonconforming "virtual node"). At one Newton-Cotes point:
//   Δu_rel = u_b − u_s = N_b·v_b − N_s·v_s ,  t_skin = T_skin·Δu_rel ,
//   K_skin = ∫ N_rel^T T_skin N_rel dS (4 blocks; Eq 7-57/58).  T_skin = k_a t̂t̂ᵀ + k_n n̂n̂ᵀ
// (beam axis). Newton-Cotes nodal integration (Table 7-6; same as our interface → node
// pairs decouple). This header produces the skin INTEGRATION POINTS (via point_location);
// assembly is at the caller (test/solver). Math: docs/references/embedded-beam-formulation.md.

#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/plate.hpp>          // plate::shape (beam shape functions)
#include <katai/fem/elements/point_location.hpp> // locate_point
#include <katai/mesh/mesh.hpp>

namespace katai::core::ebeam {

// Skin integration point: precomputed coupling data (fast assembly).
struct SkinPoint {
    std::array<int, 3> beam_node;   // beam element nodes (indices into the beam node list) [A,B,mid]
    Eigen::Vector3d Nb;             // beam shape functions N_b(ξ)
    int soil_elem = -1;             // the soil element containing the point
    double xi_s = 0.0, eta_s = 0.0; // soil local coordinates (N_s = E::shape_functions(xi_s,eta_s))
    double wJ = 0.0;                // Newton-Cotes weight × beam segment Jacobian
    Eigen::Matrix2d T;              // global skin stiffness 2×2 (k_a t̂t̂ᵀ + k_n n̂n̂ᵀ) — elastic
    Eigen::Vector2d tang;           // unit tangent t̂ (beam axis) — for the axial plasticity
    double k_a = 0.0, k_n = 0.0;    // axial / lateral skin stiffness
    double t_max = -1.0;            // AXIAL skin capacity (force/length); ≤0 ⇒ unbounded (elastic)
    bool ok = false;               // located in the soil?
};

// The beam is split into 3-node elements (each beam_elem is {A,B,mid} beam-node indices).
// beam_x/beam_y: beam node coordinates. k_axial/k_lateral: skin stiffness per unit length.
// At each element's 3 Newton-Cotes points (ξ=−1,0,+1, w=1/3,4/3,1/3) it locates the point
// in the mesh and builds T_skin from the beam axis.
inline std::vector<SkinPoint> build_skin_points(const mesh::Mesh& mesh,
                                                const std::vector<double>& beam_x,
                                                const std::vector<double>& beam_y,
                                                const std::vector<std::array<int, 3>>& beam_elem,
                                                double k_axial, double k_lateral,
                                                double t_max = -1.0) {
    std::vector<SkinPoint> pts;
    const double xi_nc[3] = {-1.0, 0.0, 1.0};
    const double w_nc[3] = {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0};
    for (const auto& el : beam_elem) {
        Eigen::Vector3d X, Y;
        for (int k = 0; k < 3; ++k) { X(k) = beam_x[el[k]]; Y(k) = beam_y[el[k]]; }
        for (int q = 0; q < 3; ++q) {
            const Eigen::Vector3d N = plate::shape(xi_nc[q]);
            const Eigen::Vector3d dN = plate::shape_deriv_xi(xi_nc[q]);
            const double dx = dN.dot(X), dy = dN.dot(Y);
            const double J = std::sqrt(dx * dx + dy * dy);       // ds/dξ
            const double px = N.dot(X), py = N.dot(Y);
            SkinPoint sp;
            sp.beam_node = el; sp.Nb = N; sp.wJ = w_nc[q] * J;
            const double c = dx / J, s = dy / J;                 // unit tangent t̂=(c,s); n̂=(−s,c)
            Eigen::Vector2d t(c, s), n(-s, c);
            sp.tang = t; sp.k_a = k_axial; sp.k_n = k_lateral; sp.t_max = t_max;
            sp.T = k_axial * (t * t.transpose()) + k_lateral * (n * n.transpose());
            const auto loc = ploc::locate_point(mesh, px, py);
            sp.soil_elem = loc.element; sp.xi_s = loc.xi; sp.eta_s = loc.eta; sp.ok = loc.found;
            pts.push_back(sp);
        }
    }
    return pts;
}

// An embedded-beam package that can be handed to the solver: beam node coordinates +
// INDEPENDENT extra DOFs (translations ux,uy + rotation φ), 3-node beam elements (the
// plate core), skin points (pre-located). Foot (tip) interaction: a single point at the
// beam tip (toe, beam node 0) — axial spring D_foot + capacity F_max. The same N_s
// coupling as the soil (Sci.Man §7.5 Eq 7-60…64).
struct FootPoint {
    int beam_node = -1;             // beam tip node (toe)
    int soil_elem = -1; double xi_s = 0.0, eta_s = 0.0;
    Eigen::Vector2d tang{0.0, 1.0}; // beam axis (axial direction)
    double D_foot = 0.0;            // tip axial spring stiffness
    double f_max = -1.0;            // tip capacity (force); ≤0 ⇒ unbounded
    bool ok = false;
};

struct EmbeddedBeam {
    std::vector<double> node_x, node_y;        // beam node coordinates
    std::vector<int> dof_x, dof_y, dof_phi;    // global extra DOFs per beam node
    std::vector<std::array<int, 3>> elements;  // beam node indices [A,B,mid]
    plate::PlateProps props;
    std::vector<SkinPoint> skin;
    FootPoint foot;                            // toe (beam node 0); D_foot≤0 ⇒ no foot
};

// `px,py`: the beam node polyline (2*ne+1 nodes, bottom to top). Adds extra DOFs to the
// DofMap (BEFORE finalize). k_axial/k_lateral: skin stiffness per unit length. props: beam
// EA/EI/ν. The skin points are located in the mesh (point_location). Returns: the package
// to put into the solver's Structures.embedded_beams.
inline EmbeddedBeam build_embedded_beam(const mesh::Mesh& mesh, DofMap& dofs,
                                        const std::vector<double>& px, const std::vector<double>& py,
                                        const plate::PlateProps& props,
                                        double k_axial, double k_lateral, double t_max = -1.0,
                                        double D_foot = 0.0, double f_max = -1.0) {
    EmbeddedBeam b;
    b.node_x = px; b.node_y = py; b.props = props;
    const int nbn = static_cast<int>(px.size());
    for (int i = 0; i < nbn; ++i) {
        b.dof_x.push_back(dofs.add_extra_dof());
        b.dof_y.push_back(dofs.add_extra_dof());
        b.dof_phi.push_back(dofs.add_extra_dof());
    }
    for (int e = 0; 2 * e + 2 < nbn; ++e) b.elements.push_back({2 * e, 2 * e + 2, 2 * e + 1});
    b.skin = build_skin_points(mesh, px, py, b.elements, k_axial, k_lateral, t_max);
    if (D_foot > 0.0) {  // foot at the toe (beam node 0, first node of the bottom-to-top polyline)
        b.foot.beam_node = 0; b.foot.D_foot = D_foot; b.foot.f_max = f_max;
        const double dx = px[1] - px[0], dy = py[1] - py[0], L = std::sqrt(dx * dx + dy * dy);
        b.foot.tang = Eigen::Vector2d(dx / L, dy / L);
        const auto loc = ploc::locate_point(mesh, px[0], py[0]);
        b.foot.soil_elem = loc.element; b.foot.xi_s = loc.xi; b.foot.eta_s = loc.eta; b.foot.ok = loc.found;
    }
    return b;
}

// PLAXIS default interface-stiffness factors for an embedded beam row
// (Sluis 2012; PLAXIS 2D Reference sec 6.6): ISF_RS = ISF_RN =
// 2.5 (L_spacing/D_eq)^-0.75 for the skin springs and ISF_KF =
// 25 (L_spacing/D_eq)^-0.75 for the foot, each times the adjacent soil's
// shear modulus G (the foot also times the equivalent radius D_eq/2).
inline void default_interface_stiffness(double spacing, double diameter, double G,
                                        double& k_axial, double& k_lateral,
                                        double& foot_stiffness) {
    const double isf = 2.5 * std::pow(spacing / diameter, -0.75);
    k_axial = isf * G;
    k_lateral = isf * G;
    foot_stiffness = 25.0 * std::pow(spacing / diameter, -0.75) * G * (0.5 * diameter);
}

}  // namespace katai::core::ebeam
