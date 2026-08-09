// The groundwater-flow driver body, compiled ONCE (section 5.2).
#include <katai/jobs/flow_driver.hpp>

#include <algorithm>
#include <cmath>

#include <katai/analysis/seepage.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/boundary_extraction.hpp>   // collect_chain (flux-edge node chains)
#include <katai/linsolve/direct_solver.hpp>

namespace katai::app {

FlowResult solve_groundwater_flow(const model::Project& pr, const katai::mesh::Mesh& mesh,
                                  const model::Phase* phase) {
    FlowResult R;
    if (mesh.element_count == 0) { R.message = "Empty mesh -- generate a mesh first (Mesh tab)."; return R; }
    if (pr.materials.empty())    { R.message = "No materials."; return R; }

    // Permeability per material [m/day]. Any material actually used by the mesh must have a
    // positive kx and ky (PLAXIS Groundwater tab); guard honestly instead of solving a singular
    // or silently-impermeable system.
    std::vector<katai::core::Permeability> perm(pr.materials.size(), {1.0, 1.0});
    std::vector<char> used(pr.materials.size(), 0);
    for (int e = 0; e < mesh.element_count; ++e) {
        const int m = mesh.element_material[e];
        if (m < 0 || m >= (int)pr.materials.size()) { R.message = "An element has no valid material."; return R; }
        used[m] = 1;
    }
    // Per-material van Genuchten/Mualem retention -- the SAME unsaturated k_rel(suction) the transient
    // + fully-coupled flow uses (build_problem.hpp), so the free-surface / phreatic location and
    // discharge are consistent between a steady "Groundwater flow" calculation and a coupled analysis
    // on the same model. In the saturated zone (psi >= 0) k_rel = 1, so confined flow is unchanged.
    std::vector<katai::core::WaterRetention> ret(pr.materials.size());
    for (size_t m = 0; m < pr.materials.size(); ++m) {
        perm[m] = {pr.materials[m].kx, pr.materials[m].ky};
        ret[m] = {pr.materials[m].gw_ga, pr.materials[m].gw_gn, pr.materials[m].gw_gl,
                  pr.materials[m].gw_Sres, 1.0};
        // NonPorous = IMPERMEABLE barrier (PLAXIS: no flow through a non-porous region).
        // Removing its elements from the assembly would leave isolated nodes and make H
        // singular; instead the conductivity is pinned to a negligible floor (a fixed
        // 1e-8 m/day, >= 1e6 times below neighbouring soil) -- the flow net routes AROUND
        // the barrier and the system stays well-posed. The user's kx/ky on this material
        // is NOT read in flow (so the k > 0 requirement is lifted too).
        // Undrained (C) joins it: a total-stress cluster takes no part in a flow calculation, and
        // PLAXIS does not offer its permeability for entry at all ("the input field for
        // permeabilities are greyed out when the drainage type is either Non-porous or
        // Undrained C"). Reading kx/ky there would be reading a number nobody was asked for.
        if (pr.materials[m].drainage == model::Drainage::NonPorous ||
            pr.materials[m].drainage == model::Drainage::UndrainedC) {
            perm[m] = {1e-8, 1e-8};
            continue;
        }
        if (used[m] && (perm[m].kx <= 0.0 || perm[m].ky <= 0.0)) {
            R.message = "Material '" + pr.materials[m].name + "' has no permeability -- set kx, ky "
                        "(> 0) on its Groundwater tab before a flow calculation.";
            return R;
        }
    }

    // Match each boundary node to the flow BCs of every polygon edge it lies on (the same
    // coincident-edge union as the deformation BCs). Priority at corners: a prescribed head wins
    // over a seepage face, which wins over closed (a reservoir corner is at reservoir head).
    std::vector<int> fixed_nodes, seepage_nodes;
    std::vector<double> fixed_values;
    for (int node : mesh.boundary_nodes) {
        const double xn = mesh.x[node], yn = mesh.y[node];
        double best_d = 1e300;
        struct Hit { int flow; double head, d; };
        std::vector<Hit> hits;
        for (const auto& P : pr.polygons) {
            const int n = (int)P.x.size();
            if (n < 3 || (int)P.edge_flow.size() != n) continue;   // no flow BCs on this polygon
            for (int i = 0; i < n; ++i) {
                const double ax = P.x[i], ay = P.y[i], bx = P.x[(i + 1) % n], by = P.y[(i + 1) % n];
                const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey; if (l2 < 1e-18) continue;
                double t = ((xn - ax) * ex + (yn - ay) * ey) / l2; t = std::clamp(t, 0.0, 1.0);
                const double d = std::hypot(xn - (ax + t * ex), yn - (ay + t * ey));
                best_d = std::fmin(best_d, d);
                const double h = (int)P.edge_head.size() == n ? P.edge_head[i] : 0.0;
                hits.push_back({P.edge_flow[i], h, d});
            }
        }
        const double tol = best_d + 1e-9 + 1e-7 * std::fmax(std::fabs(xn), std::fabs(yn));
        bool has_head = false, has_seep = false; double hval = 0.0, hd = 1e300;
        for (const auto& h : hits) {
            if (h.d > tol) continue;
            if (h.flow == (int)model::FlowBCType::Head && h.d < hd) { has_head = true; hval = h.head; hd = h.d; }
            if (h.flow == (int)model::FlowBCType::Seepage) has_seep = true;
        }
        if (has_head)      { fixed_nodes.push_back(node); fixed_values.push_back(hval); }
        else if (has_seep) { seepage_nodes.push_back(node); }
    }
    // Prescribed-flux (Neumann) edges: the manual's boundary term q (Scientific Manual Eqs.
    // 3-31/3-34). Each edge is integrated ONCE, node-indexed, and the solvers scatter it into
    // whatever equation numbering they currently have -- they rebuild it as the free surface and
    // the seepage-face active set move.
    std::vector<double> nodal_flux(mesh.node_count, 0.0);
    bool any_flux = false;
    for (const auto& P : pr.polygons) {
        const int n = (int)P.x.size();
        if (n < 3 || (int)P.edge_flow.size() != n) continue;
        for (int i = 0; i < n; ++i) {
            if (P.edge_flow[i] != (int)model::FlowBCType::Flux) continue;
            const double q = (int)P.edge_flux.size() == n ? P.edge_flux[i] : 0.0;
            if (q == 0.0) continue;   // a zero flux IS the natural (closed) condition
            const auto chain = katai::mesh::collect_chain(mesh, P.x[i], P.y[i],
                                                          P.x[(i + 1) % n], P.y[(i + 1) % n]);
            katai::core::accumulate_boundary_flux(mesh, chain, q, nodal_flux);
            any_flux = true;
        }
    }

    // Hydraulic conditions drawn INSIDE the model (PLAXIS Reference sec. 5.9). A well prescribes a
    // discharge along its line; a drain prescribes the head at its nodes. Both are one-sided in
    // general -- a well stops extracting once the head reaches h_min, a Normal drain does nothing
    // where the ground is already drier than the drain -- so both are resolved by an active set
    // around the solve rather than by a single pass.
    //
    //   level[i] = a node that MAY be held at a head, with which way it clamps:
    //     from_above (a drain): active while the free head would rise ABOVE the level
    //     from_below (a well's h_min): active while the extraction would pull the head BELOW it
    //   A Vacuum drain is `always` -- it holds its head in both directions, which is what a vacuum
    //   consolidation does to the ground.
    struct LevelNode { int node; double level; bool from_above; bool always; std::size_t owner; };
    std::vector<LevelNode> level;
    std::vector<double> well_flux(mesh.node_count, 0.0);   // the wells' own nodal discharge
    bool any_well = false, any_drain = false;
    std::string hydro_note;
    for (std::size_t hi = 0; hi < pr.hydros.size(); ++hi) {
        if (phase && !phase->active_hydro(hi)) continue;
        const model::HydroLine& H = pr.hydros[hi];
        const auto chain = katai::mesh::collect_chain(mesh, H.x1, H.y1, H.x2, H.y2);
        if (chain.size() < 2) {
            R.message = "The " + std::string(model::hydro_kind_name(H.kind)) + " \"" + H.name +
                        "\" does not lie on the mesh: its line resolved to " +
                        std::to_string(chain.size()) +
                        " node(s), so it would take water from nowhere. Draw it inside the soil.";
            return R;
        }
        if (H.kind == model::HydroKind::Drain) {
            any_drain = true;
            const bool vacuum = H.behaviour == (int)model::DrainBehaviour::Vacuum;
            for (int n : chain)
                level.push_back({n, H.head, true, vacuum, hi});
        } else {
            any_well = true;
            // The discharge is spread along the line. PLAXIS distributes it over the layers a well
            // crosses "as a function of the saturated permeability and the intersected depth"; this
            // build spreads it by length alone, which is the same statement for a well inside one
            // material and is said out loud when it is not.
            const double len = std::hypot(H.x2 - H.x1, H.y2 - H.y1);
            if (len < 1e-12) { R.message = "A well needs a line of non-zero length."; return R; }
            const double sign = H.behaviour == (int)model::WellBehaviour::Infiltration ? 1.0 : -1.0;
            katai::core::accumulate_boundary_flux(mesh, chain, sign * std::fabs(H.q) / len, well_flux);
            if (sign < 0.0)   // only an extraction well can pull the ground down to h_min
                for (int n : chain) level.push_back({n, H.h_min, false, false, hi});
            int first_mat = -1;
            bool mixed = false;
            for (int e = 0; e < mesh.element_count && !mixed; ++e) {
                bool touches = false;
                for (int k = 0; k < mesh.nodes_per_element && !touches; ++k)
                    for (int n : chain) if (mesh.node_of(e, k) == n) { touches = true; break; }
                if (!touches) continue;
                if (first_mat < 0) first_mat = mesh.element_material[e];
                else if (mesh.element_material[e] != first_mat) mixed = true;
            }
            if (mixed && hydro_note.empty())
                hydro_note = " Well \"" + H.name + "\" crosses more than one material: its "
                             "discharge is spread along its length, not weighted by each layer's "
                             "permeability.";
        }
    }
    for (int n = 0; n < mesh.node_count; ++n)
        if (well_flux[n] != 0.0) { nodal_flux[n] += well_flux[n]; any_flux = true; }

    if (fixed_nodes.empty() && !any_drain) {
        R.message = (any_flux || any_well)
                        ? "A prescribed discharge alone does not fix the head: a flow problem with "
                          "only fluxes and wells has no unique solution. Prescribe a head on at "
                          "least one boundary as well."
                        : "No flow boundary conditions: right-click a soil edge in Flow conditions "
                          "and prescribe a head on at least one boundary (everything else is closed "
                          "by default).";
        return R;
    }

    const katai::core::SeepageLinearSolve lin =
        [](const katai::math::CsrMatrix& A, const Eigen::VectorXd& b) {
            auto s = katai::linsolve::make_direct_solver(katai::linsolve::MatrixType::RealSymmetricPositiveDefinite);
            s->factorize(A);
            return s->solve(b);
        };

    // Free-surface transition width ~ half the mean element size (smaller = sharper phreatic
    // surface and more accurate discharge; too small stalls the Picard iteration). Computed
    // before the solve because the solve now happens more than once: a well or a drain is a
    // one-sided condition, and which of them are holding is not known until the head is.
    double area_sum = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        area_sum += 0.5 * std::fabs((mesh.x[b] - mesh.x[a]) * (mesh.y[c] - mesh.y[a]) -
                                    (mesh.x[c] - mesh.x[a]) * (mesh.y[b] - mesh.y[a]));
    }
    const double mean_size = std::sqrt(2.0 * area_sum / std::max(1, mesh.element_count));
    katai::core::UnconfinedOptions opt;
    opt.transition = 0.5 * mean_size;
    opt.retention = &ret;                              // van Genuchten k_rel (consistent with coupled)
    opt.relax = seepage_nodes.empty() ? 0.15 : 0.05;   // seepage face needs the stronger damping
    opt.max_iter = 1200;
    opt.tol = 1e-5;

    Eigen::VectorXd head_full;
    int iterations = 0;

    // ONE steady solve with a given Dirichlet set and a given nodal flux: the confined fast path
    // when it is self-consistent, the unconfined + seepage-face Picard iteration otherwise. With
    // no hydraulic conditions in the model this is called exactly once, with exactly the arguments
    // the driver passed before wells and drains existed.
    const auto solve_once = [&](const std::vector<int>& fx, const std::vector<double>& fv,
                                const std::vector<double>& flux, bool with_flux) -> bool {
        bool have_head = false;
        // CONFINED fast path: one linear solve with k_rel = 1 (saturated). If the resulting
        // pressure head psi = h - y is non-negative everywhere (and no seepage faces are
        // declared), the linear solution is already the exact steady state -- k_rel == 1 is
        // self-consistent, so the Picard iteration (and its accuracy floor) is unnecessary.
        if (seepage_nodes.empty()) {
            katai::core::DofMap fdofs(mesh.node_count, 1);
            std::vector<double> hp(mesh.node_count, 0.0);
            for (size_t i = 0; i < fx.size(); ++i) {
                fdofs.fix_node_component(fx[i], 0);
                hp[fx[i]] = fv[i];
            }
            fdofs.finalize();
            if (fdofs.equation_count() > 0) {
                katai::math::SparseMatrixBuilder builder(fdofs.equation_count());
                Eigen::VectorXd rhs = Eigen::VectorXd::Zero(fdofs.equation_count());
                katai::core::assemble_seepage(mesh, fdofs, perm, hp, builder, rhs);
                for (int n = 0; n < mesh.node_count; ++n) {
                    if (flux[n] == 0.0) continue;
                    const int eq = fdofs.equation(fdofs.global_dof(n, 0));
                    if (eq >= 0) rhs[eq] += flux[n];
                }
                const Eigen::VectorXd hf = lin(builder.build(), rhs);
                head_full.resize(mesh.node_count);
                double hscale = 1.0, psi_min = 1e300;
                for (int n = 0; n < mesh.node_count; ++n) {
                    const int eq = fdofs.equation(fdofs.global_dof(n, 0));
                    head_full[n] = eq >= 0 ? hf[eq] : hp[n];
                    hscale = std::fmax(hscale, std::fabs(head_full[n]));
                }
                for (int n = 0; n < mesh.node_count; ++n)
                    psi_min = std::fmin(psi_min, head_full[n] - mesh.y[n]);
                if (psi_min >= -1e-9 * hscale) { have_head = true; iterations = 1; }
            }
        }
        opt.nodal_flux.clear();
        if (with_flux) opt.nodal_flux = flux;   // the manual's boundary term q
        if (!have_head) {   // unconfined / seepage-face regime: variable-k_rel Picard + active set
            const auto res = katai::core::solve_unconfined_seepage_face(mesh, fx, fv, seepage_nodes,
                                                                        perm, lin, opt);
            iterations = res.iterations;
            if (!res.converged) {
                R.message = "Groundwater flow did not converge (residual " +
                            std::to_string(res.residual) +
                            "). Check that the prescribed heads are consistent with the geometry.";
                return false;
            }
            head_full = res.head;
        }
        return true;
    };

    // The hydraulic conditions' ACTIVE SET. A drain starts holding -- that is its ordinary state,
    // and it keeps the system non-singular when a drain is the only head in the model; a well
    // starts pumping its rated discharge and meets h_min only if the ground cannot supply it.
    std::vector<char> clamped;              // which nodes ended up held at a level
    std::vector<char> on(level.size(), 0);
    for (std::size_t i = 0; i < level.size(); ++i) on[i] = level[i].from_above ? 1 : 0;
    bool set_settled = level.empty();
    for (int pass = 0; pass < 30; ++pass) {
        std::vector<int> fx = fixed_nodes;
        std::vector<double> fv = fixed_values;
        std::vector<double> flux = nodal_flux;
        clamped.assign(mesh.node_count, 0);
        for (std::size_t i = 0; i < level.size(); ++i) {
            if (!on[i] || clamped[level[i].node]) continue;
            clamped[level[i].node] = 1;
            fx.push_back(level[i].node);
            fv.push_back(level[i].level);
        }
        // A well node held at h_min no longer receives its share of the discharge -- the manual's
        // "when the groundwater head reduces below the h_min level no further extraction will
        // occur". What the ground can still give at that head is then the clamp's own flux.
        for (int n = 0; n < mesh.node_count; ++n)
            if (clamped[n] && well_flux[n] != 0.0) flux[n] -= well_flux[n];
        if (!solve_once(fx, fv, flux, any_flux)) return R;
        if (level.empty()) break;

        const Eigen::VectorXd Qc = katai::core::compute_nodal_flux(mesh, perm, head_full, opt);
        double qscale = 1e-30;
        for (int n = 0; n < mesh.node_count; ++n) qscale = std::fmax(qscale, std::fabs(Qc[n]));
        const double qtol = 1e-9 * qscale;
        bool changed = false;
        for (std::size_t i = 0; i < level.size(); ++i) {
            if (level[i].always) continue;            // a vacuum drain holds in both directions
            const int n = level[i].node;
            const double htol = 1e-9 * std::fmax(1.0, std::fabs(level[i].level));
            if (!on[i]) {
                const bool binds = level[i].from_above ? head_full[n] > level[i].level + htol
                                                       : head_full[n] < level[i].level - htol;
                if (binds) { on[i] = 1; changed = true; }
            } else if (level[i].from_above) {
                // A Normal drain only takes water away. If holding its head would push water INTO
                // the ground, the ground is already drier than the drain and the drain does
                // nothing: "pore pressures lower than the equivalent to the given head are not
                // affected by the drain".
                if (Qc[n] > qtol) { on[i] = 0; changed = true; }
            } else {
                // The h_min clamp is a limit, not a target: once the ground can supply the well's
                // rated share at that head, the well is no longer limited and pumps it again.
                if (-Qc[n] >= std::fabs(well_flux[n]) - qtol) { on[i] = 0; changed = true; }
            }
        }
        if (!changed) { set_settled = true; break; }
    }
    if (!set_settled) {
        R.message = "The wells and drains did not settle into a steady set: their one-sided "
                    "conditions are cycling. Check that a well's h_min is below the head the "
                    "ground can supply it at, and that drains are not stacked on the same nodes.";
        return R;
    }

    R.iterations = iterations;
    R.head = head_full;

    // Discharge through the prescribed-head boundaries + global mass balance. At steady state the
    // inflow equals the outflow (incl. seepage faces); the sum over ALL nodes is zero up to solver
    // tolerance (conductivity rows sum to zero).
    const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, R.head, opt);
    double q_in = 0.0, q_all = 0.0;
    for (int n : fixed_nodes) if (Q[n] > 0.0) q_in += Q[n];
    for (int n = 0; n < mesh.node_count; ++n) q_all += Q[n];
    // Total inflow INTO the domain: through prescribed-head boundaries and through prescribed
    // flux edges alike. Reporting only the head part was complete while a flux boundary could
    // not be asked for; with one, a model whose water enters through the flux edge and leaves
    // through a head boundary would report "discharge = 0" and read as if nothing flowed. For a
    // model without flux edges the second term is exactly zero, so every existing number is
    // unchanged -- KV-FLW-001 (Charny) is the witness.
    double flux_in = 0.0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (nodal_flux[n] > 0.0) flux_in += nodal_flux[n];
    R.discharge = q_in + flux_in;
    R.balance_err = q_in > 1e-30 ? std::fabs(q_all) / q_in : 0.0;

    // What the wells and drains actually took out. This is not the same number as the discharge
    // the user asked a well for: a well limited by h_min extracts what the ground can give at
    // that head, and a Normal drain extracts whatever reaches it. Reporting only the rated
    // discharge would hide exactly the case the limit exists for.
    R.hydro_discharge = 0.0;
    R.hydro_limited = 0;
    for (int n = 0; n < mesh.node_count; ++n) {
        if (!clamped.empty() && clamped[n]) {
            if (Q[n] < 0.0) R.hydro_discharge += -Q[n];
        } else if (well_flux[n] < 0.0) {
            R.hydro_discharge += -well_flux[n];
        }
    }
    for (std::size_t i = 0; i < level.size(); ++i)
        if (!level[i].from_above && on[i]) ++R.hydro_limited;

    R.pore.assign(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        R.pore[n] = kFlowGammaWater * std::fmax(0.0, R.head[n] - mesh.y[n]);

    R.ok = true;
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "Flow solved: discharge Q = %.4g m3/day/m, %d iterations, mass balance %.2g.",
                  R.discharge, R.iterations, R.balance_err);
    R.message = buf;
    // Walls and interfaces are not flow barriers in this build (PLAXIS Scientific Manual sec. 3.4
    // gives interfaces their own flow setting: impermeable screen, semi-permeable, or drain).
    // The flow net here runs straight THROUGH a cut-off wall, which is a real difference between
    // the drawing and the calculation -- so the run says it rather than leaving it to be noticed.
    {
        bool barrier = false;
        for (const auto& st : pr.structs)
            barrier = barrier || st.kind == model::StructKind::Plate ||
                      st.kind == model::StructKind::Interface;
        if (barrier)
            R.message += " Note: walls and interfaces do not block flow in this build -- the "
                         "flow net is continuous through them. Model an impermeable barrier as a "
                         "thin Non-porous region.";
    }
    if (any_well || any_drain) {
        std::snprintf(buf, sizeof(buf), " Wells and drains removed %.4g m3/day/m.",
                      R.hydro_discharge);
        R.message += buf;
        if (R.hydro_limited > 0)
            R.message += " A well reached its minimum head h_min and is extracting less than its "
                         "rated discharge: the ground cannot supply it.";
        R.message += hydro_note;
    }
    return R;
}

}  // namespace katai::app
