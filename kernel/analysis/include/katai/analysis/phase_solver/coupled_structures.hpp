#pragma once
// WHAT A COUPLED PHASE OWES THE READER ABOUT ITS STRUCTURES -- written once, for both the
// consolidation and the fully-coupled strategies, because it is one statement and not two.
//
// Both phases solved their structural elements as one ELASTIC matrix until 2026-10, and this file
// held what that owed the reader: a plate, anchor or geogrid past its capacity (K2D-A014), a
// geogrid in compression pushing back on the soil (K2D-A015), a joint past its Coulomb capacity
// that could not slip (K2D-A016) -- each measured against the input and reported. Both phases now
// evaluate their structures by their own laws at every iterate (StructureLaws, nonlinear_solver.hpp),
// so there is no elastic branch left to declare, the three codes are retired, and the structures
// are reported as the static phase reports them.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>      // Structures
#include <katai/analysis/results.hpp>
#include <katai/analysis/structural_diagrams.hpp>   // DiagSpec / IfaceDiag / force_diagram
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// THE STRUCTURES OF A COUPLED PHASE THAT SOLVED THEM BY THEIR OWN LAWS (StructureLaws): reported
// as the static phase reports its own -- the diagrams from the total displacement and the committed
// plastic state, the joints with the slip that was actually checked. There is no elastic branch
// left to declare.
inline void report_structures_by_their_laws(const Structures& structures,
                                            const std::vector<DiagSpec>* diagrams,
                                            const std::vector<IfaceDiag>* iface_diagrams,
                                            const katai::mesh::Mesh& mesh, const DofMap& dofs,
                                            const Eigen::VectorXd& disp_total,
                                            const StructuralInit& st, SolveResult& R) {
    if (diagrams)
        for (const auto& sp : *diagrams) {
            StructForce d = force_diagram(sp, structures, mesh, dofs, disp_total, st.anchor_plastic,
                                          st.geogrid_plastic, /*elastic=*/false, st.plate_plastic,
                                          st.plate5_plastic);
            // a plate yielded: force_diagram reads it (plate_on_cap, static_phase.hpp)
            const auto env = force_envelope(d.stations);
            d.max_N = env.max_abs_N; d.max_Q = env.max_abs_Q; d.max_M = env.max_abs_M;
            R.struct_forces.push_back(std::move(d));
        }
    if (iface_diagrams)
        for (const auto& is : *iface_diagrams) {
            InterfaceResult ir = force_diagram(is, structures, mesh, dofs, disp_total,
                                               st.interface_slip, st.interface5_slip, false, true,
                                               false, nullptr, &st.interface_pore,
                                               &st.interface5_pore);
            ir.slip_checked = true;   // the Coulomb return ran on every iterate of the phase
            for (const auto& s : ir.stations) {
                ir.max_abs_tau = std::max(ir.max_abs_tau, std::fabs(s.tau));
                ir.max_abs_sigma_n = std::max(ir.max_abs_sigma_n, std::fabs(s.sigma_n));
                ir.max_abs_slip = std::max(ir.max_abs_slip, std::fabs(s.slip));
                if (s.slipping) ir.any_slip = true;
            }
            R.interface_forces.push_back(std::move(ir));
        }
}

} // namespace katai::core
