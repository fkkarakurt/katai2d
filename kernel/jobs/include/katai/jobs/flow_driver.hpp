#pragma once
// The groundwater-flow driver (layer 2, katai/jobs; historical name build_flow.hpp).
// Turn a project's groundwater-flow boundary conditions into a steady-state seepage solve
// on the FE mesh (PLAXIS "Groundwater flow" calculation). Confined and unconfined regimes are
// both handled by the variable-k_rel free-surface solver with seepage-face active sets (a fully
// saturated domain simply keeps k_rel = 1 everywhere). Returns the nodal head / pore fields for
// display and for coupling into the mechanical solve (assemble_pore_load_from_head +
// assemble_gravity_from_head), plus the boundary discharge and a mass-balance check.
//
// Core formulation + validation: docs/references/seepage-formulation.md (Darcy/Laplace, Charny,
// seepage face); the GUI-path numbers are pinned in tests/test_seepage_gui.cpp.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

namespace katai::app {

struct FlowResult {
    Eigen::VectorXd head;        // nodal hydraulic head h [m] (size node_count)
    std::vector<double> pore;    // nodal pore pressure u = gamma_w * max(0, h - y) [kPa]
    double discharge = 0.0;      // TOTAL inflow into the domain [m3/day per m]: through
                                 // prescribed-head boundaries and prescribed-flux edges alike
    double balance_err = 0.0;    // |sum of all nodal fluxes| / inflow (mass conservation, ~0)
    int iterations = 0;
    bool ok = false;
    std::string message;
};

// Water unit weight [kN/m^3] (= build_problem.hpp kGammaWater; distinct name so the two headers
// can be included in either order / independently).
inline constexpr double kFlowGammaWater = 9.81;

// Definition in kernel/jobs/src/flow_driver.cpp (section 5.2).
FlowResult solve_groundwater_flow(const model::Project& pr, const katai::mesh::Mesh& mesh);

}  // namespace katai::app
