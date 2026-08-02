#pragma once
// Parent structural-state carry (Track 1a; engine-owned since Stage B9).
// Soil carries its committed stresses through `initial_state`, but structural
// elements are total-displacement formulated (f = f(u_total, plastic state)) --
// a chained phase inherits the parent's structural state only if it is given
// (a) the converged displacement DATUM and (b) the committed plastic states.
// Without them every phase restarts its structures at u = 0 with zero plastic
// memory, and the SumMstage imbalance re-ramps the parent's structural
// tractions -- measured: an unchanged nil phase drifted the wall moment by 32%.
//
// Every size must match EXACTLY (same mesh, same DofMap layout, same structure
// set); anything else returns false and the consumer falls back -- honestly,
// named in its message -- to the old increment-from-zero path rather than seed
// a half-matching state.

#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>  // Structures, StructuralInit
#include <katai/analysis/results.hpp>           // StructCarryState
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/geogrid.hpp>       // geogrid::kGaussCount
#include <katai/fem/elements/interface.hpp>     // iface::kPointCount(5)
#include <katai/fem/elements/plate.hpp>         // plate::kPlasticStateSize(5)

namespace katai::core {

// Fill `out` from the parent's carried structural state. Returns true only when
// every member matches the current structure set exactly; on false, `out` is
// untouched and the caller must not consume the carry.
inline bool build_structural_init(const StructCarryState& ps, const Structures& structures,
                                  const DofMap& dofs, StructuralInit& out) {
    size_t skin_total = 0;
    for (const auto& eb : structures.embedded_beams) skin_total += eb.skin.size();
    if (ps.full_disp.size() != (Eigen::Index)dofs.total_dofs() ||
        ps.anchor_plastic.size() != structures.anchors.size() ||
        ps.geogrid_plastic.size() !=
            structures.geogrids.size() * (size_t)geogrid::kGaussCount ||
        ps.interface_slip.size() !=
            structures.interfaces.size() * (size_t)iface::kPointCount ||
        ps.interface5_slip.size() !=
            structures.interfaces5.size() * (size_t)iface::kPointCount5 ||
        ps.embedded_skin_slip.size() != skin_total ||
        ps.embedded_foot_slip.size() != structures.embedded_beams.size() ||
        ps.plate_plastic.size() !=
            structures.plates.size() * (size_t)plate::kPlasticStateSize ||
        ps.plate5_plastic.size() !=
            structures.plates5.size() * (size_t)plate::kPlasticStateSize5)
        return false;
    // Datum in EQUATION space: the solvers' structural loops read free DOFs
    // only (a fixed DOF contributes zero there, exactly as it does in every
    // static solve).
    out.u_datum = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int g = 0; g < dofs.total_dofs(); ++g) {
        const int eq = dofs.equation(g);
        if (eq >= 0) out.u_datum(eq) = ps.full_disp(g);
    }
    out.anchor_plastic = ps.anchor_plastic;
    out.geogrid_plastic = ps.geogrid_plastic;
    out.interface_slip = ps.interface_slip;
    out.interface5_slip = ps.interface5_slip;
    out.embedded_skin_slip = ps.embedded_skin_slip;
    out.embedded_foot_slip = ps.embedded_foot_slip;
    out.plate_plastic = ps.plate_plastic;
    out.plate5_plastic = ps.plate5_plastic;
    return true;
}

} // namespace katai::core
