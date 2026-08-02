#pragma once
// Embedded retaining WALL (barrier) builder — split_mesh_at_wall + plate + interfaces on
// both sides. A wall is a flexural barrier separating the soil: active behind (retained),
// passive in front (excavated).
//   - the mesh is SPLIT along the wall line with split_mesh_at_wall (right/left detached,
//     continuous below the toe).
//   - the wall's INDEPENDENT translational+rotational extra DOFs (add_extra_dof) → detached
//     from the mesh (like a PLAXIS plate).
//   - the plate element lives on the wall DOFs (geometry from the right mesh nodes); two
//     interfaces (right soil↔wall, left soil↔wall) tie the wall to both soils (slip +
//     separation possible).
//   - the TOE (embedment tip, lowest below-seam node) is SHARED: wall translational DOF =
//     toe mesh DOF → the wall stays tied to the embedment soil (from below). The seam nodes
//     above are detached.
// Math: structural-plate-formulation.md (§18 plate) + interface-formulation.md.

#include <algorithm>
#include <vector>

#include <katai/analysis/initial_stress.hpp>       // K0Options
#include <katai/analysis/nonlinear_solver.hpp>     // PlateElement, InterfaceElement
#include <katai/analysis/staged_construction.hpp>  // SeamPair, split_mesh_at_wall
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

struct WallBuild {
    std::vector<PlateElement> plates;
    std::vector<InterfaceElement> interfaces;  // right + left (2 per plate element)
    // Wall positions (bottom to top; pos0 = TOE, shared): for recovery.
    std::vector<double> y;          // position elevation
    std::vector<int> dof_x, dof_y;  // wall translational global DOFs (toe: mesh DOF; seam: extra DOF)
    std::vector<int> dof_phi;       // wall rotational extra DOFs
    std::vector<int> node_r;        // right (retained) mesh node (toe: shared)
};

// `seam` = split_mesh_at_wall output; `toe_node` = the lowest SHARED mesh node on the wall
// line (on the y_toe corner row, not split). dofs BEFORE finalize. pp/ip: wall + interface
// properties.
inline WallBuild build_embedded_wall(const mesh::Mesh& mesh, std::vector<SeamPair> seam,
                                     int toe_node, DofMap& dofs,
                                     const plate::PlateProps& pp,
                                     const iface::InterfaceProps& ip) {
    std::sort(seam.begin(), seam.end(),
              [](const SeamPair& a, const SeamPair& b) { return a.y < b.y; });
    WallBuild w;
    // Position list: toe (shared, corner) + seam (bottom to top). seam[0]=mid, seam[1]=corner,
    // ... (y_toe on a corner row) → the positions cycle corner,mid,corner,...
    std::vector<int> nr, nl, dx, dy, dphi;  // position → right node, left node, DOFs
    nr.push_back(toe_node); nl.push_back(toe_node);            // TOE shared
    dx.push_back(dofs.global_dof(toe_node, 0));               // wall translation = toe mesh DOF (tied)
    dy.push_back(dofs.global_dof(toe_node, 1));
    dphi.push_back(dofs.add_extra_dof());
    for (const auto& s : seam) {
        nr.push_back(s.right); nl.push_back(s.left);
        dx.push_back(dofs.add_extra_dof());                   // INDEPENDENT wall translation
        dy.push_back(dofs.add_extra_dof());
        dphi.push_back(dofs.add_extra_dof());
    }
    const int npos = static_cast<int>(nr.size());
    for (int i = 0; i < npos; ++i) {
        w.y.push_back(mesh.y[nr[i]]); w.dof_x.push_back(dx[i]); w.dof_y.push_back(dy[i]);
        w.dof_phi.push_back(dphi[i]); w.node_r.push_back(nr[i]);
    }
    // Plate + interface elements: consecutive triples [A=2e(corner), B=2e+2(corner), mid=2e+1].
    for (int e = 0; 2 * e + 2 < npos; ++e) {
        const int a = 2 * e, m = 2 * e + 1, b = 2 * e + 2;
        PlateElement pe;
        pe.nodes = {nr[a], nr[b], nr[m]};                     // geometry: right mesh nodes
        pe.rot_dof = {dphi[a], dphi[b], dphi[m]};
        pe.trans_dof = {dx[a], dy[a], dx[b], dy[b], dx[m], dy[m]};  // independent wall translation
        pe.props = pp;
        w.plates.push_back(pe);
        const std::array<int, 6> sdof = {dx[a], dy[a], dx[b], dy[b], dx[m], dy[m]};
        InterfaceElement ir;  // right (retained) soil ↔ wall — normal outward from the right soil (−x)
        ir.soil_nodes = {nr[a], nr[b], nr[m]}; ir.struct_dof = sdof; ir.props = ip;
        w.interfaces.push_back(ir);
        // Left (excavated) soil ↔ wall. Node order REVERSED (A↔B) → the edge tangent, hence
        // the normal, flips → the normal points OUTWARD from the left soil (+x) (outward-
        // normal convention symmetric with the right). Thus compression is σ_n<0 on both
        // sides (tension cut-off correct; K0 σ_n0 is seeded by a single formula, no per-side
        // sign needed). struct_dof also swaps A↔B → node pairing preserved.
        InterfaceElement il;
        il.soil_nodes = {nl[b], nl[a], nl[m]};
        il.struct_dof = {dx[b], dy[b], dx[a], dy[a], dx[m], dy[m]};
        il.props = ip;
        w.interfaces.push_back(il);
    }
    return w;
}

// SEEDS the wall interfaces with the K0 geostatic horizontal stress (σ_n0 = K0·σ'_v at each
// Newton-Cotes point, by depth). In the split mesh the K0 horizontal stress is discontinuous
// across the wall line; σ_n0 carries it at Δu_n=0 → the "wished-in-place" wall is in
// self-equilibrium at u=0, NO spurious installation movement. Since build_embedded_wall
// flips the left interface normal, compression is σ_n<0 on both sides (tension cut-off
// correct, single formula). In the excavation phase the staged-release load applies only
// the excavation imbalance. (PLAXIS K0 procedure → interface initial stress;
// interface-formulation.md §6.) Call AFTER build_embedded_wall, BEFORE the solve.
inline void seed_interface_k0(WallBuild& w, const mesh::Mesh& mesh, const K0Options& opt) {
    const auto ncp = iface::nc_points();
    for (auto& ie : w.interfaces)
        for (int q = 0; q < 3; ++q) {
            const double y = mesh.y[ie.soil_nodes[ncp[q].node]];
            const double sv = -opt.unit_weight * (opt.surface_elevation - y);  // σ'_v (compression<0)
            ie.sigma_n0[q] = opt.k0 * sv;                                      // σ_n0 = K0·σ'_v
        }
}

// Wall (plate elements) internal-force envelope: gathers each plate element's 9-DOF vector
// from the global solution `disp` (translations=trans_dof, rotations=rot_dof global
// indices), computes N,Q,M at a few ξ (plate::forces, EI·κ). Returns: max |M| (bending
// moment, kN·m/m), max |Q|, max |N|. For wall deflection/moment verification (PLAXIS sheet-
// pile benchmark). `disp` = NewtonResult.displacement (total_dofs).
struct WallForceEnvelope { double max_abs_M = 0.0, max_abs_Q = 0.0, max_abs_N = 0.0; };
inline WallForceEnvelope wall_force_envelope(const WallBuild& w, const mesh::Mesh& mesh,
                                             const Eigen::VectorXd& disp) {
    WallForceEnvelope env;
    const std::array<double, 5> xis{-1.0, -0.5, 0.0, 0.5, 1.0};
    for (const auto& pe : w.plates) {
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]]; }
        plate::Dof u = plate::Dof::Zero();
        for (int k = 0; k < 3; ++k) {
            u(3 * k + 0) = disp[pe.trans_dof[2 * k + 0]];
            u(3 * k + 1) = disp[pe.trans_dof[2 * k + 1]];
            u(3 * k + 2) = disp[pe.rot_dof[k]];
        }
        for (double xi : xis) {
            const auto f = plate::forces(X, pe.props, u, xi);
            env.max_abs_M = std::max(env.max_abs_M, std::fabs(f.M));
            env.max_abs_Q = std::max(env.max_abs_Q, std::fabs(f.Q));
            env.max_abs_N = std::max(env.max_abs_N, std::fabs(f.N));
        }
    }
    return env;
}

// ===========================================================================================
// tri15 (5-node structural) WALL — counterpart of the PLAXIS 15-node soil + 5-node
// plate/interface. Same architecture as build_embedded_wall; the wall-line nodes are
// quarter-spaced (tri15), the plate5/interface5 elements are consecutive groups of 5
// [4e..4e+4] (one macro-element height).
// ===========================================================================================
struct WallBuild5 {
    std::vector<PlateElement5> plates;
    std::vector<InterfaceElement5> interfaces;     // right + left (2 per plate)
    std::vector<double> y;
    std::vector<int> dof_x, dof_y, dof_phi, node_r;
};

inline WallBuild5 build_embedded_wall5(const mesh::Mesh& mesh, std::vector<SeamPair> seam,
                                       int toe_node, DofMap& dofs,
                                       const plate::PlateProps& pp,
                                       const iface::InterfaceProps& ip) {
    std::sort(seam.begin(), seam.end(),
              [](const SeamPair& a, const SeamPair& b) { return a.y < b.y; });
    WallBuild5 w;
    std::vector<int> nr, nl, dx, dy, dphi;
    nr.push_back(toe_node); nl.push_back(toe_node);          // TOE shared (bottom to top)
    dx.push_back(dofs.global_dof(toe_node, 0));
    dy.push_back(dofs.global_dof(toe_node, 1));
    dphi.push_back(dofs.add_extra_dof());
    for (const auto& s : seam) {
        nr.push_back(s.right); nl.push_back(s.left);
        dx.push_back(dofs.add_extra_dof());
        dy.push_back(dofs.add_extra_dof());
        dphi.push_back(dofs.add_extra_dof());
    }
    const int npos = static_cast<int>(nr.size());
    for (int i = 0; i < npos; ++i) {
        w.y.push_back(mesh.y[nr[i]]); w.dof_x.push_back(dx[i]); w.dof_y.push_back(dy[i]);
        w.dof_phi.push_back(dphi[i]); w.node_r.push_back(nr[i]);
    }
    // plate5/interface5: consecutive groups of 5 [4e, 4e+1, 4e+2, 4e+3, 4e+4] (ends shared).
    for (int e = 0; 4 * e + 4 < npos; ++e) {
        std::array<int, 5> p;  // position indices (bottom to top = natural ξ=−1..+1)
        for (int k = 0; k < 5; ++k) p[k] = 4 * e + k;
        PlateElement5 pe;
        for (int k = 0; k < 5; ++k) {
            pe.nodes[k] = nr[p[k]];
            pe.rot_dof[k] = dphi[p[k]];
            pe.trans_dof[2 * k + 0] = dx[p[k]];
            pe.trans_dof[2 * k + 1] = dy[p[k]];
        }
        pe.props = pp;
        w.plates.push_back(pe);
        InterfaceElement5 ir;  // right soil ↔ wall (normal outward from the right soil)
        for (int k = 0; k < 5; ++k) {
            ir.soil_nodes[k] = nr[p[k]];
            ir.struct_dof[2 * k + 0] = dx[p[k]];
            ir.struct_dof[2 * k + 1] = dy[p[k]];
        }
        ir.props = ip;
        w.interfaces.push_back(ir);
        InterfaceElement5 il;  // left soil ↔ wall: node order REVERSED → normal outward from the left soil
        for (int k = 0; k < 5; ++k) {
            const int pk = p[4 - k];                       // reversed order
            il.soil_nodes[k] = nl[pk];
            il.struct_dof[2 * k + 0] = dx[pk];
            il.struct_dof[2 * k + 1] = dy[pk];
        }
        il.props = ip;
        w.interfaces.push_back(il);
    }
    return w;
}

inline void seed_interface_k0(WallBuild5& w, const mesh::Mesh& mesh, const K0Options& opt) {
    const auto ncp = iface::nc_points5();
    for (auto& ie : w.interfaces)
        for (int q = 0; q < 5; ++q) {
            const double y = mesh.y[ie.soil_nodes[ncp[q].node]];
            const double sv = -opt.unit_weight * (opt.surface_elevation - y);
            ie.sigma_n0[q] = opt.k0 * sv;
        }
}

inline WallForceEnvelope wall_force_envelope(const WallBuild5& w, const mesh::Mesh& mesh,
                                             const Eigen::VectorXd& disp) {
    WallForceEnvelope env;
    const std::array<double, 9> xis{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    for (const auto& pe : w.plates) {
        plate::NodeCoords5 X;
        for (int k = 0; k < 5; ++k) { X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]]; }
        plate::Dof5 u = plate::Dof5::Zero();
        for (int k = 0; k < 5; ++k) {
            u(3 * k + 0) = disp[pe.trans_dof[2 * k + 0]];
            u(3 * k + 1) = disp[pe.trans_dof[2 * k + 1]];
            u(3 * k + 2) = disp[pe.rot_dof[k]];
        }
        for (double xi : xis) {
            const auto f = plate::forces5(X, pe.props, u, xi);
            env.max_abs_M = std::max(env.max_abs_M, std::fabs(f.M));
            env.max_abs_Q = std::max(env.max_abs_Q, std::fabs(f.Q));
            env.max_abs_N = std::max(env.max_abs_N, std::fabs(f.N));
        }
    }
    return env;
}

}  // namespace katai::core
