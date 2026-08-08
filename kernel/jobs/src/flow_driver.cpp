// The groundwater-flow driver body, compiled ONCE (section 5.2).
#include <katai/jobs/flow_driver.hpp>

#include <algorithm>
#include <cmath>

#include <katai/analysis/seepage.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/boundary_extraction.hpp>   // collect_chain (flux-edge node chains)
#include <katai/linsolve/direct_solver.hpp>

namespace katai::app {

FlowResult solve_groundwater_flow(const model::Project& pr, const katai::mesh::Mesh& mesh) {
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
        if (pr.materials[m].drainage == model::Drainage::NonPorous) {
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

    if (fixed_nodes.empty()) {
        R.message = any_flux
                        ? "A prescribed flux alone does not fix the head: a flow problem with only "
                          "Neumann boundaries has no unique solution. Prescribe a head on at least "
                          "one boundary as well."
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

    Eigen::VectorXd head_full;
    int iterations = 0;
    bool have_head = false;

    // CONFINED fast path: one linear solve with k_rel = 1 (saturated). If the resulting pressure
    // head psi = h - y is non-negative everywhere (and no seepage faces are declared), the linear
    // solution is already the exact steady state -- k_rel == 1 is self-consistent, so the Picard
    // iteration (and its under-relaxation accuracy floor) is unnecessary.
    if (seepage_nodes.empty()) {
        katai::core::DofMap fdofs(mesh.node_count, 1);
        std::vector<double> hp(mesh.node_count, 0.0);
        for (size_t i = 0; i < fixed_nodes.size(); ++i) {
            fdofs.fix_node_component(fixed_nodes[i], 0);
            hp[fixed_nodes[i]] = fixed_values[i];
        }
        fdofs.finalize();
        if (fdofs.equation_count() > 0) {
            katai::math::SparseMatrixBuilder builder(fdofs.equation_count());
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(fdofs.equation_count());
            katai::core::assemble_seepage(mesh, fdofs, perm, hp, builder, rhs);
            for (int n = 0; n < mesh.node_count; ++n) {
                if (nodal_flux[n] == 0.0) continue;
                const int eq = fdofs.equation(fdofs.global_dof(n, 0));
                if (eq >= 0) rhs[eq] += nodal_flux[n];
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

    // Free-surface transition width ~ half the mean element size (smaller = sharper phreatic
    // surface and more accurate discharge; too small stalls the Picard iteration).
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
    if (any_flux) opt.nodal_flux = nodal_flux;   // the manual's boundary term q

    if (!have_head) {   // unconfined / seepage-face regime: variable-k_rel Picard + active set
        const auto res = katai::core::solve_unconfined_seepage_face(mesh, fixed_nodes, fixed_values,
                                                                    seepage_nodes, perm, lin, opt);
        iterations = res.iterations;
        if (!res.converged) {
            R.message = "Groundwater flow did not converge (residual " + std::to_string(res.residual) +
                        "). Check that the prescribed heads are consistent with the geometry.";
            return R;
        }
        head_full = res.head;
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

    R.pore.assign(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        R.pore[n] = kFlowGammaWater * std::fmax(0.0, R.head[n] - mesh.y[n]);

    R.ok = true;
    char buf[160];
    std::snprintf(buf, sizeof(buf),
                  "Flow solved: discharge Q = %.4g m3/day/m, %d iterations, mass balance %.2g.",
                  R.discharge, R.iterations, R.balance_err);
    R.message = buf;
    return R;
}

}  // namespace katai::app
