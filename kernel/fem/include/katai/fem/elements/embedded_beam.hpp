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

// How the beam's CONNECTION POINT is attached (PLAXIS 2D Ref. Man. sec 5.6.3). PLAXIS also
// offers Rigid (to a plate, rotation coupled too); that is out of scope here and declared in
// docs/references/embedded-beam-formulation.md sec 8.
enum class Connection {
    // "The displacement at the connection point of the beam is directly coupled with the
    // displacement of the element in which the connection point is located ... they undergo
    // exactly the same displacement, but not necessarily in the same rotation." PLAXIS's
    // DEFAULT when no structure shares the point.
    Hinged = 0,
    // "Not directly coupled with the soil element in which the beam top is located, but the
    // interaction through the interface elements is still present."
    Free = 1,
};

struct EmbeddedBeam {
    std::vector<double> node_x, node_y;        // beam node coordinates
    std::vector<int> dof_x, dof_y, dof_phi;    // global extra DOFs per beam node
    std::vector<std::array<int, 3>> elements;  // beam node indices [A,B,mid]
    plate::PlateProps props;
    std::vector<SkinPoint> skin;
    FootPoint foot;                            // toe (beam node 0); D_foot≤0 ⇒ no foot
    // Hinged connection: the beam node at the connection point carries NO translation DOFs of
    // its own (dof_x/dof_y are −1 there) — they ARE the tied mesh node's. The constraint
    // u_b,top = N_s·v_s is therefore satisfied exactly, as a degree-of-freedom identity, because
    // the mesher carries the connection point as a vertex, so N_s collapses to a unit vector.
    // No interpolation, no Lagrange multiplier, no penalty stiffness. The shaft stays
    // mesh-nonconforming — only this one point is a node. See the formulation document sec 7.1.
    Connection connection = Connection::Free;
    int conn_beam_node = -1;                   // beam node tied (−1 = none / Free)
    int conn_mesh_node = -1;                   // the mesh node it is tied to
};

// The global DOF of beam node `k`'s translation component `comp` (0 = x, 1 = y): the beam's own
// extra DOF, or — at a hinged connection point, where it has none — the tied mesh node's. The
// exact counterpart of plate_trans_eq, and for the same reason.
inline int trans_gdof(const EmbeddedBeam& b, int k, int comp, const DofMap& dofs) {
    const int own = comp == 0 ? b.dof_x[k] : b.dof_y[k];
    if (own >= 0) return own;
    if (b.conn_mesh_node < 0) return -1;   // no DOF and nothing to tie to: carries nothing
    return dofs.global_dof(b.conn_mesh_node, comp);
}
inline int trans_eq(const EmbeddedBeam& b, int k, int comp, const DofMap& dofs) {
    const int g = trans_gdof(b, k, comp, dofs);
    return g >= 0 ? dofs.equation(g) : -1;
}

// `px,py`: the beam node polyline (2*ne+1 nodes, bottom to top). Adds extra DOFs to the
// DofMap (BEFORE finalize). k_axial/k_lateral: skin stiffness per unit length. props: beam
// EA/EI/ν. The skin points are located in the mesh (point_location). Returns: the package
// to put into the solver's Structures.embedded_beams.
// `conn_beam_node`/`conn_mesh_node`: a HINGED connection point (>= 0 both) ties that beam node's
// translations to that mesh node — the node the mesher carries at the connection point — so they
// get NO extra DOFs of their own. The rotation stays the beam's own: hinged couples displacement,
// "but not necessarily the same rotation" (Ref. Man. sec 5.6.3). Pass −1/−1 for a Free connection,
// which allocates exactly as before and is bit-identical to the pre-connection engine.
inline EmbeddedBeam build_embedded_beam(const mesh::Mesh& mesh, DofMap& dofs,
                                        const std::vector<double>& px, const std::vector<double>& py,
                                        const plate::PlateProps& props,
                                        double k_axial, double k_lateral, double t_max = -1.0,
                                        double D_foot = 0.0, double f_max = -1.0,
                                        int conn_beam_node = -1, int conn_mesh_node = -1) {
    EmbeddedBeam b;
    b.node_x = px; b.node_y = py; b.props = props;
    const int nbn = static_cast<int>(px.size());
    const bool hinged = conn_beam_node >= 0 && conn_beam_node < nbn && conn_mesh_node >= 0;
    if (hinged) {
        b.connection = Connection::Hinged;
        b.conn_beam_node = conn_beam_node;
        b.conn_mesh_node = conn_mesh_node;
    }
    for (int i = 0; i < nbn; ++i) {
        const bool tied = hinged && i == conn_beam_node;
        b.dof_x.push_back(tied ? -1 : dofs.add_extra_dof());
        b.dof_y.push_back(tied ? -1 : dofs.add_extra_dof());
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

// PLAXIS default interface stiffness for an embedded beam row (Sluis 2012; PLAXIS 2D
// Reference Manual sec 6.6.4). TWO equations, and the second one is the reason this
// function takes the spacing twice over:
//
//   Eq 6-66  the dimensionless FACTORS, from the spacing-to-diameter ratio:
//              ISF_RS = ISF_RN = 2.5 (L_spacing/D)^-0.75 ,  ISF_KF = 25 (L_spacing/D)^-0.75
//            where D is the real diameter (width for a square pile).
//   Eq 6-65  the STIFFNESSES themselves, tied to the surrounding soil's shear modulus:
//              R_S = ISF_RS G/L_spacing ,  R_N = ISF_RN G/L_spacing ,
//              K_F = ISF_KF G R_eq/L_spacing .
//
// The division by L_spacing in Eq 6-65 is what smears one pile of a row over a metre of
// wall, exactly as the driver divides EA, EI, the pile weight and both capacities. It is
// also what makes the units come out: a plane-strain skin traction is a force per metre of
// beam per metre of wall, so T_skin = t/u is kN/m^3, and ISF*G alone would be kN/m^2.
//
// R_eq is NOT half the diameter. Eq 6-67 defines the equivalent diameter from the section
// itself, D_eq = sqrt(12 EI/EA), which for a solid circular pile is 0.866 D rather than D
// -- the ratio EI/EA is the same whether the per-pile or the per-metre-of-wall values are
// passed, so the caller may use either.
//
// Validity, stated by the manual and not by us: these defaults were derived for BORED piles
// loaded statically in the AXIAL direction, in Hardening Soil with small-strain stiffness,
// with the phreatic level at the ground surface. Away from those conditions they are a
// starting point, not a prediction.
inline void default_interface_stiffness(double spacing, double diameter, double G,
                                        double EA, double EI,
                                        double& k_axial, double& k_lateral,
                                        double& foot_stiffness) {
    const double ratio = std::pow(spacing / diameter, -0.75);
    const double D_eq = (EA > 0.0 && EI > 0.0) ? std::sqrt(12.0 * EI / EA) : diameter;
    k_axial = 2.5 * ratio * G / spacing;
    k_lateral = 2.5 * ratio * G / spacing;
    foot_stiffness = 25.0 * ratio * G * (0.5 * D_eq) / spacing;
}

}  // namespace katai::core::ebeam
