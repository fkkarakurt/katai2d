#pragma once
// GENERAL interface builder — arbitrarily oriented, in-soil (standalone) OR along any
// structural line. `split_mesh_at_segment` (staged_construction.hpp) PAIRS the mesh nodes
// along the interface line: the original nodes stay tied to the POSITIVE side of the
// segment, the twins to the NEGATIVE side. This builder makes one `InterfaceElement`
// (tri6, 3 nodes) / `InterfaceElement5` (tri15, 5 nodes) per interface edge;
// soil_nodes = the positive-side (original) nodes, struct_dof = the MESH DOFs of the
// NEGATIVE-side (twin) nodes (not extra DOFs — the twins are real mesh nodes). Thus the
// interface ties the two soil sides with k_n,k_s + Coulomb (slip + separation). Difference
// from the embedded wall (embedded_wall.hpp): there the "structure side" is a plate on
// independent extra DOFs; here it is the soil of the other side (or it can attach to a
// plate). In the split mesh the K0 horizontal stress is discontinuous →
// `seed_soil_interface_k0` seeds σ_n0 (orientation-aware) (wished-in-place equilibrium, no
// spurious installation movement). Math: interface-formulation.md §6; staged-construction
// (split).

#include <array>
#include <cmath>
#include <vector>

#include <katai/analysis/initial_stress.hpp>       // K0Options
#include <katai/analysis/nonlinear_solver.hpp>     // Structures, InterfaceElement(5)
#include <katai/analysis/staged_construction.hpp>  // SegSeam, split_mesh_at_segment
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// An interface's [begin,end) slice inside `structures` (for K0 seeding / ownership).
struct InterfaceRange { int order; std::size_t begin; std::size_t end; };

// Builds soil-soil interface elements from the seam (sorted by s). order=6 →
// InterfaceElement (3-node edge: corner,mid,corner), order=15 → InterfaceElement5 (5-node:
// corner,q,mid,q,corner). struct_dof = the mesh DOFs of the twin nodes (from the
// post-split dofs; global_dof is stable before/after dofs.finalize).
inline InterfaceRange build_soil_interface(const mesh::Mesh& mesh, const std::vector<SegSeam>& seam,
                                           const DofMap& dofs, const iface::InterfaceProps& ip,
                                           int order, Structures& structures) {
    if (order == 15) {
        const std::size_t begin = structures.interfaces5.size();
        for (std::size_t e = 0; 4 * e + 4 < seam.size(); ++e) {
            const SegSeam* sp = &seam[4 * e];     // 5 nodes [ξ=−1,−0.5,0,0.5,1]
            InterfaceElement5 ie;
            for (int k = 0; k < 5; ++k) {
                ie.soil_nodes[k] = sp[k].orig;
                ie.struct_dof[2 * k + 0] = dofs.global_dof(sp[k].dup, 0);
                ie.struct_dof[2 * k + 1] = dofs.global_dof(sp[k].dup, 1);
            }
            ie.props = ip;
            structures.interfaces5.push_back(ie);
        }
        return {15, begin, structures.interfaces5.size()};
    }
    const std::size_t begin = structures.interfaces.size();
    for (std::size_t e = 0; 2 * e + 2 < seam.size(); ++e) {
        const SegSeam& A = seam[2 * e], &M = seam[2 * e + 1], &B = seam[2 * e + 2];
        InterfaceElement ie;
        ie.soil_nodes = {A.orig, B.orig, M.orig};   // [A, B, mid] (iface NodeCoords order)
        ie.struct_dof = {dofs.global_dof(A.dup, 0), dofs.global_dof(A.dup, 1),
                         dofs.global_dof(B.dup, 0), dofs.global_dof(B.dup, 1),
                         dofs.global_dof(M.dup, 0), dofs.global_dof(M.dup, 1)};
        ie.props = ip;
        structures.interfaces.push_back(ie);
    }
    return {6, begin, structures.interfaces.size()};
}

// SEEDS the interface σ_n0 with the K0 geostatic stress (orientation-aware). Geostatic
// σ_xx=K0σ'_v, σ_yy=σ'_v, σ_xy=0 ⇒ on a plane with normal (nx,ny), σ_n = σ'_v·(K0·nx² + ny²)
// (compression<0). Vertical interface (n=±x) → K0σ'_v (identical to seed_interface_k0);
// horizontal (n=±y) → σ'_v. Carries the K0 discontinuity across the split line at Δu_n=0 →
// wished-in-place equilibrium. `(nx,ny)` = the interface unit normal (from the segment
// tangent).
inline void seed_soil_interface_k0(Structures& structures, const mesh::Mesh& mesh, const K0Options& opt,
                                   double nx, double ny, const InterfaceRange& rng) {
    const double fac = opt.k0 * nx * nx + ny * ny;  // σ_n / σ'_v
    if (rng.order == 15) {
        const auto ncp = iface::nc_points5();
        for (std::size_t i = rng.begin; i < rng.end; ++i)
            for (int q = 0; q < 5; ++q) {
                const double y = mesh.y[structures.interfaces5[i].soil_nodes[ncp[q].node]];
                const double sv = -opt.unit_weight * (opt.surface_elevation - y);
                structures.interfaces5[i].sigma_n0[q] = fac * sv;
            }
        return;
    }
    const auto ncp = iface::nc_points();
    for (std::size_t i = rng.begin; i < rng.end; ++i)
        for (int q = 0; q < 3; ++q) {
            const double y = mesh.y[structures.interfaces[i].soil_nodes[ncp[q].node]];
            const double sv = -opt.unit_weight * (opt.surface_elevation - y);
            structures.interfaces[i].sigma_n0[q] = fac * sv;
        }
}

}  // namespace katai::core
