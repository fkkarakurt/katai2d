#pragma once
// Wished-in-place interface prestress baseline (Stage B9). An embedded wall or
// a standalone Coulomb joint is seeded with the geostatic lateral earth
// pressure sigma_n0 = (K0 nx^2 + ny^2) sigma'_v; the K0 baseline B must carry
// the matching internal force, or residual(0) != 0 and the wall relaxes the
// geostatic pressure into a spurious "installation" deflection before any
// excavation -- biasing displacements and N/Q/M. The 3-node (tri6) and 5-node
// (tri15) joints carry the SAME seed, integrated at their Newton-Cotes points.
//
// Callers on the parent-carry path must SKIP this: structural_internal_force
// at the parent state subsumes the sigma_n0 terms (coulomb_return adds
// sigma_n0 internally), and adding both would double-count the wished-in-place
// lateral earth pressure.

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>  // Structures
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>     // iface:: NC points + edge frames
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Add every interface's sigma_n0 internal force to the baseline vector B
// (equation space).
inline void add_interface_sigma_n0_baseline(const Structures& structures, const mesh::Mesh& mesh,
                                            const DofMap& dofs, Eigen::VectorXd& B) {
    const auto ncp = iface::nc_points();
    for (const auto& ie : structures.interfaces) {   // 3-node (tri6) interface sigma_n0
        iface::NodeCoords Xe;
        for (int k = 0; k < 3; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
        for (int q = 0; q < 3; ++q) {
            const int nd = ncp[q].node;
            const auto fr = iface::edge_frame(Xe, ncp[q].xi);
            const double bcoef[4] = {fr.s, -fr.c, -fr.s, fr.c};
            const double wJ = ncp[q].w * fr.J;
            const int idx[4] = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                dofs.equation(ie.struct_dof[2 * nd + 0]),
                                dofs.equation(ie.struct_dof[2 * nd + 1])};
            for (int i = 0; i < 4; ++i) if (idx[i] >= 0) B(idx[i]) += wJ * bcoef[i] * ie.sigma_n0[q];
        }
    }
    const auto ncp5 = iface::nc_points5();
    for (const auto& ie : structures.interfaces5) {  // 5-node (tri15): same seed, 5-pt NC
        iface::NodeCoords5 Xe;
        for (int k = 0; k < 5; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
        for (int q = 0; q < 5; ++q) {
            const int nd = ncp5[q].node;
            const auto fr = iface::edge_frame5(Xe, ncp5[q].xi);
            const double bcoef[4] = {fr.s, -fr.c, -fr.s, fr.c};
            const double wJ = ncp5[q].w * fr.J;
            const int idx[4] = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                dofs.equation(ie.struct_dof[2 * nd + 0]),
                                dofs.equation(ie.struct_dof[2 * nd + 1])};
            for (int i = 0; i < 4; ++i) if (idx[i] >= 0) B(idx[i]) += wJ * bcoef[i] * ie.sigma_n0[q];
        }
    }
}

}  // namespace katai::core
