#pragma once
// Structural internal-force OUTPUTS (PLAXIS Output counterpart) — the DISTRIBUTION of N
// (axial), Q (shear), M (bending moment) + location along plate/anchor/geogrid elements,
// from the converged global solution.
//
// Pure POST-PROCESSING: mirrors the solver's constitutive/return mapping EXACTLY (the
// anchor/geogrid committed plastic state comes from NewtonResult). NOT tied to MKL — only
// `disp` (NewtonResult.displacement, full global-DOF space, fixed DOFs = prescribed/0) +
// DofMap. The GENERAL + DISTRIBUTED counterpart of wall_force_envelope, which returned
// only max |M|,|Q|,|N|.
//
// PLAXIS Reference Manual "Output → Structures → Bending moments M / Shear forces Q / Axial forces N":
// structural forces along the element, per unit width (kN·m/m, kN/m). Sign convention as in
// plate.hpp: N tension-positive, M = EI·κ, Q = kGA'·γ.
// Math: docs/references/structural-plate-formulation.md §9.

#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/nonlinear_solver.hpp>  // PlateElement(5), AnchorElement, GeogridElement, InterfaceElement(5)
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/embedded_beam.hpp>  // ebeam::EmbeddedBeam (pile row)
#include <katai/fem/elements/geogrid.hpp>
#include <katai/fem/elements/interface.hpp>      // iface Coulomb return (interface slip/tau/sigma_n)
#include <katai/fem/elements/plate.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Internal force + physical location at one station along an element. Q=M=0 for
// anchor/geogrid.
struct ForceStation {
    double s = 0.0;            // arc length along the element chain (from the first node, cumulative)
    double x = 0.0, y = 0.0;   // physical location (undeformed)
    double N = 0.0, Q = 0.0, M = 0.0;  // axial (tension+), shear, bending moment (per unit width)
    // NOTE: ux,uy LAST — keeps existing positional {s,x,y,N,Q,M} aggregate-inits (results_io/tests) intact.
    double ux = 0.0, uy = 0.0; // displacement at the station (interpolated with the SAME shape fns
                               // as the location) — for drawing the diagram on the deformed mesh
                               // (GUI overlay; NOT persisted)
};

// Max |N|,|Q|,|M| envelope of a force diagram (output compatible with wall_force_envelope).
struct ForceEnvelope { double max_abs_N = 0.0, max_abs_Q = 0.0, max_abs_M = 0.0; };
inline ForceEnvelope force_envelope(const std::vector<ForceStation>& diag) {
    ForceEnvelope e;
    for (const auto& st : diag) {
        e.max_abs_N = std::max(e.max_abs_N, std::fabs(st.N));
        e.max_abs_Q = std::max(e.max_abs_Q, std::fabs(st.Q));
        e.max_abs_M = std::max(e.max_abs_M, std::fabs(st.M));
    }
    return e;
}

namespace detail {
// Global translational DOF of a plate element: trans_dof[idx]≥0 if given (embedded wall,
// independent); otherwise shared with the mesh node → global_dof(node, comp) (bonded
// plate-in-soil).
inline int plate_trans_gdof(int trans, int node, int comp, const DofMap& dofs) {
    return trans >= 0 ? trans : dofs.global_dof(node, comp);
}

// Lagrange interpolation: evaluate at xi the polynomial through the NR points (xs, ys).
template <int NR>
inline double lagrange_eval(const std::array<double, NR>& xs, const std::array<double, NR>& ys,
                            double xi) {
    double v = 0.0;
    for (int i = 0; i < NR; ++i) {
        double Li = 1.0;
        for (int j = 0; j < NR; ++j)
            if (j != i) Li *= (xi - xs[j]) / (xs[i] - xs[j]);
        v += Li * ys[i];
    }
    return v;
}

// BARLOW shear recovery (3-node): with selective reduced integration the shear strain is
// superconvergent only at the REDUCED Gauss points (ξ=±1/√3); it oscillates spuriously at
// nodes/ends (Zienkiewicz&Taylor; Prathap "optimal stress points"). Sample Q at those 2
// points, interpolate linearly.
inline std::array<double, 2> plate_reduced_shear(const plate::NodeCoords& X,
                                                 const plate::PlateProps& p, const plate::Dof& u) {
    constexpr double g2 = 0.5773502691896257;  // 1/√3
    return {plate::forces(X, p, u, -g2).Q, plate::forces(X, p, u, g2).Q};
}
inline constexpr std::array<double, 2> kPlateShearXi{-0.5773502691896257, 0.5773502691896257};

// 5-node: shear is optimal at the 4-point Gauss rule (Barlow); interpolate with a cubic Lagrange.
inline std::array<double, 4> plate5_reduced_shear(const plate::NodeCoords5& X,
                                                  const plate::PlateProps& p, const plate::Dof5& u) {
    constexpr double a = 0.3399810435848563, b = 0.8611363115940526;
    return {plate::forces5(X, p, u, -b).Q, plate::forces5(X, p, u, -a).Q,
            plate::forces5(X, p, u, a).Q, plate::forces5(X, p, u, b).Q};
}
inline constexpr std::array<double, 4> kPlate5ShearXi{
    -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
}  // namespace detail

// N,Q,M diagram along a 3-node plate (wall/beam) chain. `plates` is consecutive (ends
// sharing geometry; e.g. WallBuild.plates bottom to top). Samples each element at
// `stations_per_elem` equally spaced ξ ∈ [−1,+1], accumulating arc length from the physical
// locations. Works both for an embedded wall (independent trans_dof) and bonded
// plate-in-soil (mesh-node sharing, trans_dof=−1).
// `plastic_state`: the chain's committed M-N hinge state ([ε_p,κ_p]×3 Gauss × element; the
// slice of NewtonResult.plate_plastic belonging to this chain). NON-EMPTY + props.plastic()
// ⇒ N,M are expanded to the stations by Lagrange from the return-mapped values at the GAUSS
// (stress) points (PLAXIS: "extrapolation of the values at the stress points"; a station
// value may slightly exceed the cap — PLAXIS does not check the node either). EMPTY ⇒
// elastic recovery (D6b: a linear dynamic envelope that solves the plate ELASTICALLY must
// not be silently clipped by a capped report).
inline std::vector<ForceStation> plate_force_diagram(
    const std::vector<PlateElement>& plates, const mesh::Mesh& mesh, const DofMap& dofs,
    const Eigen::VectorXd& disp, const std::vector<double>& plastic_state = {},
    int stations_per_elem = 5) {
    std::vector<ForceStation> out;
    if (stations_per_elem < 2) stations_per_elem = 2;
    const bool has_ps = plastic_state.size() == plates.size() * plate::kPlasticStateSize;
    double s_acc = 0.0;
    bool have_prev = false;
    double px = 0.0, py = 0.0;
    for (size_t ei = 0; ei < plates.size(); ++ei) {
        const auto& pe = plates[ei];
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]]; }
        plate::Dof u = plate::Dof::Zero();
        for (int k = 0; k < 3; ++k) {
            u(3 * k + 0) = disp[detail::plate_trans_gdof(pe.trans_dof[2 * k + 0], pe.nodes[k], 0, dofs)];
            u(3 * k + 1) = disp[detail::plate_trans_gdof(pe.trans_dof[2 * k + 1], pe.nodes[k], 1, dofs)];
            u(3 * k + 2) = disp[pe.rot_dof[k]];
        }
        const bool ep_pl = has_ps && pe.props.plastic();
        std::array<double, 3> Ng{}, Mg{};
        if (ep_pl)
            plate::gauss_forces_plastic(X, pe.props, u,
                                        plastic_state.data() + ei * plate::kPlasticStateSize, Ng, Mg);
        const auto Qr = detail::plate_reduced_shear(X, pe.props, u);  // Barlow shear samples
        for (int q = 0; q < stations_per_elem; ++q) {
            const double xi = -1.0 + 2.0 * q / (stations_per_elem - 1);
            const Eigen::Vector3d N = plate::shape(xi);
            double x = 0.0, y = 0.0, sx = 0.0, sy = 0.0;
            for (int k = 0; k < 3; ++k) { x += N(k) * X(k, 0); y += N(k) * X(k, 1);
                                          sx += N(k) * u(3 * k + 0); sy += N(k) * u(3 * k + 1); }
            if (have_prev) s_acc += std::hypot(x - px, y - py);
            const auto f = plate::forces(X, pe.props, u, xi);
            const double Nst = ep_pl ? detail::lagrange_eval<3>(plate::kBendGaussXi, Ng, xi) : f.N;
            const double Mst = ep_pl ? detail::lagrange_eval<3>(plate::kBendGaussXi, Mg, xi) : f.M;
            const double Q = detail::lagrange_eval<2>(detail::kPlateShearXi, Qr, xi);  // optimal Q
            out.push_back({s_acc, x, y, Nst, Q, Mst, sx, sy});
            px = x; py = y; have_prev = true;
        }
    }
    return out;
}

// Same diagram for a 5-node quartic plate (tri15 edge) chain.
inline std::vector<ForceStation> plate_force_diagram(
    const std::vector<PlateElement5>& plates, const mesh::Mesh& mesh, const DofMap& dofs,
    const Eigen::VectorXd& disp, const std::vector<double>& plastic_state = {},
    int stations_per_elem = 9) {
    std::vector<ForceStation> out;
    if (stations_per_elem < 2) stations_per_elem = 2;
    const bool has_ps = plastic_state.size() == plates.size() * plate::kPlasticStateSize5;
    double s_acc = 0.0;
    bool have_prev = false;
    double px = 0.0, py = 0.0;
    for (size_t ei = 0; ei < plates.size(); ++ei) {
        const auto& pe = plates[ei];
        plate::NodeCoords5 X;
        for (int k = 0; k < 5; ++k) { X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]]; }
        plate::Dof5 u = plate::Dof5::Zero();
        for (int k = 0; k < 5; ++k) {
            u(3 * k + 0) = disp[detail::plate_trans_gdof(pe.trans_dof[2 * k + 0], pe.nodes[k], 0, dofs)];
            u(3 * k + 1) = disp[detail::plate_trans_gdof(pe.trans_dof[2 * k + 1], pe.nodes[k], 1, dofs)];
            u(3 * k + 2) = disp[pe.rot_dof[k]];
        }
        const bool ep_pl = has_ps && pe.props.plastic();
        std::array<double, 5> Ng{}, Mg{};
        if (ep_pl)
            plate::gauss_forces_plastic5(X, pe.props, u,
                                         plastic_state.data() + ei * plate::kPlasticStateSize5, Ng, Mg);
        const auto Qr = detail::plate5_reduced_shear(X, pe.props, u);  // Barlow shear samples
        for (int q = 0; q < stations_per_elem; ++q) {
            const double xi = -1.0 + 2.0 * q / (stations_per_elem - 1);
            const Eigen::Matrix<double, 5, 1> N = plate::detail::shape5(xi);
            double x = 0.0, y = 0.0, sx = 0.0, sy = 0.0;
            for (int k = 0; k < 5; ++k) { x += N(k) * X(k, 0); y += N(k) * X(k, 1);
                                          sx += N(k) * u(3 * k + 0); sy += N(k) * u(3 * k + 1); }
            if (have_prev) s_acc += std::hypot(x - px, y - py);
            const auto f = plate::forces5(X, pe.props, u, xi);
            const double Nst = ep_pl ? detail::lagrange_eval<5>(plate::kBendGaussXi5, Ng, xi) : f.N;
            const double Mst = ep_pl ? detail::lagrange_eval<5>(plate::kBendGaussXi5, Mg, xi) : f.M;
            const double Q = detail::lagrange_eval<4>(detail::kPlate5ShearXi, Qr, xi);  // optimal Q
            out.push_back({s_acc, x, y, Nst, Q, Mst, sx, sy});
            px = x; py = y; have_prev = true;
        }
    }
    return out;
}

// Anchor AXIAL force — mirrors the solver (nonlinear_solver.cpp anchor loop) EXACTLY:
//   U = Σ g_i·u  (g = ±unit direction),  N = N0 + kk·(U − U_p),  N ∈ [−Fmax_comp, +Fmax_tens],
//   N0 = the lock-off prestress (0 for a slack anchor).
// kk = EA/L (L≤0 ⇒ geometric distance). U_p = committed plastic elongation
// (NewtonResult.anchor_plastic[i]; 0 for elastic/purely elastic anchors). Fixed DOFs are
// excluded from U, as in the solver. yielded: at capacity.
// `elastic=true` → the cap is NOT applied (raw N = kk·U): for the LINEAR dynamic envelope —
// that system solves the anchor elastically; a capped report would silently CLIP at Fmax
// and imply "yield" in a non-yielding analysis (D6b rule: the post-processor's constitutive
// law must MATCH the solver's).
struct AnchorForce { double N = 0.0; bool yielded = false; };
inline AnchorForce anchor_force(const AnchorElement& an, const mesh::Mesh& mesh, const DofMap& dofs,
                                const Eigen::VectorXd& disp, double Up = 0.0, bool elastic = false) {
    const Eigen::Vector2d Xa(mesh.x[an.node_a], mesh.y[an.node_a]);
    const Eigen::Vector2d Xb = an.node_b >= 0
        ? Eigen::Vector2d(mesh.x[an.node_b], mesh.y[an.node_b]) : an.fixed_point;
    const Eigen::Vector2d dvec = Xb - Xa;
    const double Lgeom = dvec.norm();
    if (Lgeom < 1e-30) return {};
    const Eigen::Vector2d dir = dvec / Lgeom;
    const double kk = an.EA / (an.L > 0.0 ? an.L : Lgeom);
    const int gdof[4] = {dofs.global_dof(an.node_a, 0), dofs.global_dof(an.node_a, 1),
                         an.node_b >= 0 ? dofs.global_dof(an.node_b, 0) : -1,
                         an.node_b >= 0 ? dofs.global_dof(an.node_b, 1) : -1};
    const double g[4] = {-dir(0), -dir(1), dir(0), dir(1)};
    double U = 0.0;
    for (int i = 0; i < 4; ++i)
        if (gdof[i] >= 0 && !dofs.is_fixed(gdof[i])) U += g[i] * disp[gdof[i]];
    double N = an.prestress + kk * (U - Up);
    AnchorForce r;
    if (elastic) { r.N = N; return r; }   // linear dynamic branch: uncapped elastic N (same as solver)
    const double Ft = an.Fmax_tens, Fc = an.Fmax_comp;
    if (Ft > 0.0 && N > Ft)        { r.N = Ft;  r.yielded = true; }
    else if (Fc > 0.0 && N < -Fc)  { r.N = -Fc; r.yielded = true; }
    else                           { r.N = N; }
    return r;
}

// Geogrid axial-force diagram — mirrors the constitutive/return mapping
// (geogrid::axial_return) EXACTLY: N = clamp(EA(ε − ε_p), 0, N_p); compression cut-off
// (slack) + N_p yield. The plastic state is defined at the 2 Gauss points (like PLAXIS
// stress points) → the diagram is reported at those 2 points. `ep_committed` = this
// geogrid's 2 committed plastic ε (NewtonResult.geogrid_plastic[2*gi .. 2*gi+1]);
// elastic/pure tension-only → {0,0}.
// `elastic=true` → the return mapping is NOT applied (raw N = EA·ε, in compression too):
// for the LINEAR dynamic envelope — that system solves the geogrid fully-EA elastic
// ("geogrids stay elastic in compression"); a cut-off/capped report would silently alter
// the force the solver produced (same D6b rule as `elastic` in anchor_force).
inline std::vector<ForceStation> geogrid_force_diagram(
    const GeogridElement& ge, const mesh::Mesh& mesh, const DofMap& dofs,
    const Eigen::VectorXd& disp, const std::array<double, 2>& ep_committed = {0.0, 0.0},
    bool elastic = false) {
    geogrid::NodeCoords X;
    geogrid::Dof u = geogrid::Dof::Zero();
    for (int k = 0; k < 3; ++k) {
        X(k, 0) = mesh.x[ge.nodes[k]]; X(k, 1) = mesh.y[ge.nodes[k]];
        u(2 * k + 0) = disp[dofs.global_dof(ge.nodes[k], 0)];
        u(2 * k + 1) = disp[dofs.global_dof(ge.nodes[k], 1)];
    }
    const Eigen::Vector2d A(X(0, 0), X(0, 1));
    const auto gxi = geogrid::gauss_xi();
    std::vector<ForceStation> out;
    for (int q = 0; q < geogrid::kGaussCount; ++q) {
        const double xi = gxi[q];
        const auto kin = geogrid::axial_kin(X, xi);
        const double eps = (kin.Be * u)(0);
        geogrid::AxialReturn ret;
        if (elastic) ret = {ge.props.EA * eps, ep_committed[q], ge.props.EA};
        else         ret = geogrid::axial_return(ge.props, eps, ep_committed[q]);
        // Location: quadratic shape (geogrid node order A=ξ−1, B=ξ+1, middle=ξ0 — same as plate).
        const Eigen::Vector3d Nsh = plate::shape(xi);
        double x = 0.0, y = 0.0, sx = 0.0, sy = 0.0;
        for (int k = 0; k < 3; ++k) { x += Nsh(k) * X(k, 0); y += Nsh(k) * X(k, 1);
                                      sx += Nsh(k) * u(2 * k + 0); sy += Nsh(k) * u(2 * k + 1); }
        ForceStation st;
        st.x = x; st.y = y; st.ux = sx; st.uy = sy;
        st.s = std::hypot(x - A(0), y - A(1));
        st.N = ret.N;
        out.push_back(st);
    }
    return out;
}

// Embedded beam (pile row) N,Q,M diagram along the pile (PLAXIS Output -> Structures, embedded
// beam forces). The beam is a chain of 3-node Timoshenko elements on its OWN independent DOFs
// (beam.dof_x/dof_y/dof_phi are global DOF indices into `disp`) with its OWN node coordinates
// (beam.node_x/node_y, NOT mesh nodes) -> it reuses the validated plate::forces kernel + Barlow
// shear, exactly like plate_force_diagram but reading the beam's geometry/DOFs instead of the mesh.
// `disp` is the converged full global-DOF solution (NewtonResult.displacement).
inline std::vector<ForceStation> embedded_beam_force_diagram(
    const ebeam::EmbeddedBeam& beam, const Eigen::VectorXd& disp, int stations_per_elem = 5) {
    std::vector<ForceStation> out;
    if (stations_per_elem < 2) stations_per_elem = 2;
    double s_acc = 0.0;
    bool have_prev = false;
    double px = 0.0, py = 0.0;
    for (const auto& el : beam.elements) {
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = beam.node_x[el[k]]; X(k, 1) = beam.node_y[el[k]]; }
        plate::Dof u = plate::Dof::Zero();
        for (int k = 0; k < 3; ++k) {
            u(3 * k + 0) = disp[beam.dof_x[el[k]]];
            u(3 * k + 1) = disp[beam.dof_y[el[k]]];
            u(3 * k + 2) = disp[beam.dof_phi[el[k]]];
        }
        const auto Qr = detail::plate_reduced_shear(X, beam.props, u);  // Barlow shear samples
        for (int q = 0; q < stations_per_elem; ++q) {
            const double xi = -1.0 + 2.0 * q / (stations_per_elem - 1);
            const Eigen::Vector3d N = plate::shape(xi);
            double x = 0.0, y = 0.0, sx = 0.0, sy = 0.0;
            for (int k = 0; k < 3; ++k) { x += N(k) * X(k, 0); y += N(k) * X(k, 1);
                                          sx += N(k) * u(3 * k + 0); sy += N(k) * u(3 * k + 1); }
            if (have_prev) s_acc += std::hypot(x - px, y - py);
            const auto f = plate::forces(X, beam.props, u, xi);
            const double Q = detail::lagrange_eval<2>(detail::kPlateShearXi, Qr, xi);  // optimal Q
            out.push_back({s_acc, x, y, f.N, Q, f.M, sx, sy});
            px = x; py = y; have_prev = true;
        }
    }
    return out;
}

// --- INTERFACE (soil-structure Coulomb joint) output (PLAXIS Output -> Interfaces) --------------
// Along an interface, at the Newton-Cotes node pairs: shear stress tau, normal (effective)
// stress sigma_n (tension-pos; compression < 0), relative slip = du_s, relative opening
// gap = du_n, and whether at the Coulomb limit (slipping). Mirrors the solver's
// coulomb_return EXACTLY (converged disp + committed plastic slip from NewtonResult).
// Math: docs/references/interface-formulation.md.
struct InterfaceStation {
    double s = 0.0;                 // arc length along the interface
    double x = 0.0, y = 0.0;        // physical location (undeformed)
    double tau = 0.0, sigma_n = 0.0;  // shear / normal effective stress [kPa]
    double slip = 0.0, gap = 0.0;     // relative slip du_s / opening du_n [m]
    bool slipping = false;            // reached the Coulomb limit (plastic slip)
    // Seismic (dynamic) phase: Coulomb demand/capacity ratio |τ| / τ_max. The dynamic
    // interface is solved ELASTICALLY (cannot slip) → this ratio says what the solution
    // could NEVER discover on its own: if >1, the real interface WOULD have slipped and the
    // linear result is on the unsafe side there. Meaningful only when the static state is
    // superposed (τ_max depends on the total σ_n). 0 = not computed. NEW FIELD LAST
    // (results_io is positional).
    double utilisation = 0.0;
};

namespace detail {
template <int NN, class Iface, class NodeC>
std::vector<InterfaceStation> interface_diagram_impl(
    const std::vector<Iface>& ifaces, std::size_t begin, std::size_t end, const mesh::Mesh& mesh,
    const DofMap& dofs, const Eigen::VectorXd& disp, const std::vector<double>& slip_committed,
    const std::array<iface::NCPoint, NN>& ncp, iface::EdgeFrame (*frame)(const NodeC&, double),
    bool elastic = false) {
    std::vector<InterfaceStation> out;
    double s_acc = 0.0, px = 0.0, py = 0.0; bool have_prev = false;
    for (std::size_t ii = begin; ii < end && ii < ifaces.size(); ++ii) {
        const auto& ie = ifaces[ii];
        NodeC Xe;
        for (int k = 0; k < NN; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
        for (int q = 0; q < NN; ++q) {
            const int nd = ncp[q].node;
            const auto fr = frame(Xe, ncp[q].xi);
            const double c = fr.c, s = fr.s;
            const double a[4] = {-c, -s, c, s}, b[4] = {s, -c, -s, c};  // du_s / du_n vs [soil x,y; struct x,y]
            const double u[4] = {disp[dofs.global_dof(ie.soil_nodes[nd], 0)],
                                 disp[dofs.global_dof(ie.soil_nodes[nd], 1)],
                                 disp[ie.struct_dof[2 * nd + 0]],
                                 disp[ie.struct_dof[2 * nd + 1]]};
            double du_s = 0.0, du_n = 0.0;
            for (int i = 0; i < 4; ++i) { du_s += a[i] * u[i]; du_n += b[i] * u[i]; }
            const std::size_t si = ii * NN + q;
            const double slip_c = si < slip_committed.size() ? slip_committed[si] : 0.0;
            const double x = mesh.x[ie.soil_nodes[nd]], y = mesh.y[ie.soil_nodes[nd]];
            if (have_prev) s_acc += std::hypot(x - px, y - py);
            InterfaceStation st;
            st.s = s_acc; st.x = x; st.y = y;
            st.slip = du_s; st.gap = du_n;
            if (elastic) {
                // ELASTIC report: τ = k_s·Δu_s, σ_n = k_n·Δu_n — NO Coulomb cap, NO tension
                // cut-off, NO sigma_n0. For the dynamic (seismic) phase: that system solves
                // the interface ELASTICALLY (k_n,k_s) and starts from zero (not at rest) →
                // the reported stress must be in EQUILIBRIUM with the computed
                // displacements. Applying the Coulomb return here would clip τ at capacity
                // and σ_n0 would add a static quantity; both are things the solver never
                // SAW = silently inconsistent output.
                st.tau = ie.props.ks * du_s;
                st.sigma_n = ie.props.kn * du_n;
                st.slipping = false;
            } else {
                const auto ret = iface::coulomb_return(ie.props, du_s, du_n, slip_c, ie.sigma_n0[q]);
                st.tau = ret.tau; st.sigma_n = ret.sigma_n; st.slipping = (ret.Ds == 0.0);
            }
            out.push_back(st); px = x; py = y; have_prev = true;
        }
    }
    return out;
}
}  // namespace detail

// Stations along a 3-node (tri6) interface chain [begin,end). elastic=true → the Coulomb
// return is NOT applied (τ=k_s·Δu_s, σ_n=k_n·Δu_n): only for a system that solves the
// interface ELASTICALLY (dynamic/seismic phase).
inline std::vector<InterfaceStation> interface_force_diagram(
    const std::vector<InterfaceElement>& ifaces, std::size_t begin, std::size_t end,
    const mesh::Mesh& mesh, const DofMap& dofs, const Eigen::VectorXd& disp,
    const std::vector<double>& slip_committed, bool elastic = false) {
    return detail::interface_diagram_impl<iface::kNodeCount, InterfaceElement, iface::NodeCoords>(
        ifaces, begin, end, mesh, dofs, disp, slip_committed, iface::nc_points(), &iface::edge_frame,
        elastic);
}
// 5-node (tri15) interface chain.
inline std::vector<InterfaceStation> interface_force_diagram(
    const std::vector<InterfaceElement5>& ifaces, std::size_t begin, std::size_t end,
    const mesh::Mesh& mesh, const DofMap& dofs, const Eigen::VectorXd& disp,
    const std::vector<double>& slip_committed, bool elastic = false) {
    return detail::interface_diagram_impl<iface::kNodeCount5, InterfaceElement5, iface::NodeCoords5>(
        ifaces, begin, end, mesh, dofs, disp, slip_committed, iface::nc_points5(), &iface::edge_frame5,
        elastic);
}

}  // namespace katai::core
