#pragma once
// The force of the excess pore pressure a committed Gauss state carries.
//
// A staged Plastic phase ramps f - B about a baseline B, and the ramp is the phase's change only if
// B is the internal force the phase STARTS with. The soil's part of that baseline was assembled
// from the effective stress alone (assemble_internal_force), which was right while no phase could
// start with an excess pore pressure the solver did not re-derive itself. Now the pressure is Gauss
// state (GaussState::eps_vol_und, pw_carried) and a phase can inherit one -- an undrained phase
// followed by one that ignores undrained behaviour, or a Plastic phase after a consolidation. Left
// out of B, that pressure is ramped as a load: in the undrained phase the water's stiffness keeps
// the detour small, but in a phase without it the column is first pushed up by the whole pressure
// and yields in tension on the way (measured: "equilibrated 0% of the applied load" in a phase that
// changes nothing). This adds the pressure's force, evaluated with the expression the internal
// force uses, so residual(0) = 0 holds by the parent's own equilibrium.
//
// The coupled phases (consolidation, fully coupled) keep the effective-stress baseline: their
// imbalance f - B is what re-generates the parent's pressure in their own pore field at t = 0+.

#include <vector>

#include <Eigen/Core>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Adds, in equation space, the integral of B^T (pw m) over every active element, where
// pw = pw_carried + (Kw/n) eps_vol_und at each Gauss point (the second term for undrained materials
// only, with Kw/n from the stress point's profiled material). Axisymmetric kinematics r-weight the
// integral and include the hoop component, as the internal force does.
void add_excess_pore_force(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<MaterialModel>& materials,
                           const std::vector<MaterialProfile>& profile,
                           const std::vector<GaussState>& states,
                           const std::vector<char>& active_element, bool axisymmetric,
                           Eigen::VectorXd& rhs);

}  // namespace katai::core
