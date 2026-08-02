// The transient-flow bodies (W1 saturated backward-Euler + W2 Richards
// mass-conservative modified Picard), compiled ONCE (section 5.2 batch 3c): the
// element-templated detail implementations instantiate here for tri6 and tri15;
// consumers see declarations only.
#include <katai/analysis/transient_flow.hpp>

namespace katai::core {

namespace detail {
template <class E>
TransientFlowResult transient_flow_impl(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& specific_storage, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head,
    const std::vector<double>* flux, double dt, int nsteps,
    const TransientFlowSolveFactory& solve_factory) {
    constexpr int N = E::kNodeCount;
    const int nc = mesh.node_count;

    // Free-node numbering (Dirichlet nodes excluded).
    std::vector<int> feq(nc, -1);
    int nf = 0;
    for (int n = 0; n < nc; ++n)
        if (is_prescribed.empty() || !is_prescribed[n]) feq[n] = nf++;
    auto prescribed = [&](int n) { return !is_prescribed.empty() && is_prescribed[n]; };

    // Assembly of the full (node x node) S (storage) + A = S + dt.H (system) -- one Gauss loop.
    math::SparseMatrixBuilder Sbuild(nc), Abuild(nc);
    Sbuild.reserve((std::size_t)mesh.element_count * N * N);
    Abuild.reserve((std::size_t)mesh.element_count * N * N);
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        typename E::NodeCoords X;
        for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        const Permeability& pm = perm[mesh.element_material[e]];
        const double Ss = specific_storage[mesh.element_material[e]];
        Eigen::Matrix<double, N, N> He = Eigen::Matrix<double, N, N>::Zero();
        Eigen::Matrix<double, N, N> Se = Eigen::Matrix<double, N, N>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
            const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
            const double w = gp[g].weight * sd.det_jacobian;
            Eigen::Matrix<double, 1, N> gx, gy;
            for (int i = 0; i < N; ++i) { gx(0, i) = sd.B(0, 2 * i); gy(0, i) = sd.B(1, 2 * i + 1); }
            He.noalias() += w * (pm.kx * gx.transpose() * gx + pm.ky * gy.transpose() * gy);
            Se.noalias() += w * Ss * (Nsh * Nsh.transpose());
        }
        for (int a = 0; a < N; ++a) {
            const int na = mesh.node_of(e, a);
            for (int b = 0; b < N; ++b) {
                const int nb = mesh.node_of(e, b);
                Sbuild.add_entry(na, nb, Se(a, b));
                Abuild.add_entry(na, nb, Se(a, b) + dt * He(a, b));
            }
        }
    }
    const math::CsrMatrix Sfull = Sbuild.build();
    const math::CsrMatrix Afull = Abuild.build();

    // A_ff (free x free) -- the block to factorize. Afull's free rows/columns are re-indexed.
    math::SparseMatrixBuilder Affb(std::max(1, nf));
    for (int r = 0; r < nc; ++r) {
        if (feq[r] < 0) continue;
        for (int idx = Afull.row_ptr[r]; idx < Afull.row_ptr[r + 1]; ++idx) {
            const int c = Afull.col_indices[idx];
            if (feq[c] >= 0) Affb.add_entry(feq[r], feq[c], Afull.values[idx]);
        }
    }
    const math::CsrMatrix Aff = Affb.build();

    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> solve;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu;
    if (solve_factory) {
        solve = solve_factory(Aff);
    } else {
        Eigen::MatrixXd Ad = Eigen::MatrixXd::Zero(std::max(1, nf), std::max(1, nf));
        for (int r = 0; r < nf; ++r)
            for (int idx = Aff.row_ptr[r]; idx < Aff.row_ptr[r + 1]; ++idx)
                Ad(r, Aff.col_indices[idx]) = Aff.values[idx];
        lu.compute(Ad);
        solve = [&lu](const Eigen::VectorXd& b) { return lu.solve(b); };
    }

    Eigen::VectorXd q = Eigen::VectorXd::Zero(nc);
    if (flux) for (int n = 0; n < nc; ++n) q(n) = (*flux)[n];

    // Initial state h^0 = initial_head; t=0 is recorded.
    Eigen::VectorXd hfull(nc);
    for (int n = 0; n < nc; ++n) hfull(n) = initial_head[n];
    TransientFlowResult res;
    res.times.push_back(0.0);
    res.head.push_back(hfull);

    for (int step = 0; step < nsteps; ++step) {
        const double t = (step + 1) * dt;
        // Prescribed-value extension g (0 at free nodes, head_bc(n,t) at Dirichlet nodes).
        Eigen::VectorXd g = Eigen::VectorXd::Zero(nc);
        for (int n = 0; n < nc; ++n) if (prescribed(n)) g(n) = head_bc(n, t);
        // RHS_free = (S.h^n + dt.q)_f - (A.g)_f   [A.g on the free rows = A_fp.g_p, the lift].
        const Eigen::VectorXd bfull = Sfull * hfull + dt * q;
        const Eigen::VectorXd lift = Afull * g;
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(std::max(1, nf));
        for (int n = 0; n < nc; ++n) if (feq[n] >= 0) rhs(feq[n]) = bfull(n) - lift(n);
        const Eigen::VectorXd hf = (nf > 0) ? solve(rhs) : Eigen::VectorXd::Zero(0);
        for (int n = 0; n < nc; ++n) hfull(n) = prescribed(n) ? g(n) : hf(feq[n]);
        res.times.push_back(t);
        res.head.push_back(hfull);
    }
    return res;
}
}  // namespace detail

TransientFlowResult solve_transient_flow(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& specific_storage, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head, double dt, int nsteps,
    const std::vector<double>* flux,
    const TransientFlowSolveFactory& solve_factory) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::transient_flow_impl<Tri15Element>(mesh, perm, specific_storage, is_prescribed,
                                                         head_bc, initial_head, flux, dt, nsteps, solve_factory);
    return detail::transient_flow_impl<Tri6Element>(mesh, perm, specific_storage, is_prescribed,
                                                    head_bc, initial_head, flux, dt, nsteps, solve_factory);
}

namespace detail {
template <class E>
UnsaturatedFlowResult unsaturated_flow_impl(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& porosity, const std::vector<double>& ss_sat,
    const std::vector<WaterRetention>& retention, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head,
    const std::vector<double>* flux, double dt, int nsteps,
    const TransientFlowSolveFactory& solve_factory, int max_picard, double picard_tol) {
    constexpr int N = E::kNodeCount;
    const int nc = mesh.node_count;

    std::vector<int> feq(nc, -1);
    int nf = 0;
    for (int n = 0; n < nc; ++n)
        if (is_prescribed.empty() || !is_prescribed[n]) feq[n] = nf++;
    auto prescribed = [&](int n) { return !is_prescribed.empty() && is_prescribed[n]; };
    const auto gp = E::gauss_points();

    // Node->material (for the output saturation only; last element wins -- visualisation only).
    std::vector<int> node_mat(nc, 0);
    for (int e = 0; e < mesh.element_count; ++e)
        for (int k = 0; k < N; ++k) node_mat[mesh.node_of(e, k)] = mesh.element_material[e];

    Eigen::VectorXd qv = Eigen::VectorXd::Zero(nc);
    if (flux) for (int n = 0; n < nc; ++n) qv(n) = (*flux)[n];

    // CONSISTENT assembly (theta, k_rel, capacitance at the Gauss point -- no tri6
    // corner-lumping degeneracy). A = H + M_C/dt (node^2), H (node^2, residual/flux),
    // storage_vec (node), water volume = int theta. Suction psi = y-h.
    //   storage_vec_i = (1/dt)[ int N_i(theta-theta_prev) + (int N' S_s.S_e N.(h-h_prev))_i ],
    //   theta = n.S.
    auto assemble = [&](const Eigen::VectorXd& h, const Eigen::VectorXd& hprev,
                        math::SparseMatrixBuilder& Ab, math::SparseMatrixBuilder& Hb,
                        Eigen::VectorXd& storage_vec, double& water_vol) {
        storage_vec.setZero(nc);
        water_vol = 0.0;
        for (int e = 0; e < mesh.element_count; ++e) {
            typename E::NodeCoords X;
            for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
            const Permeability& pm = perm[mesh.element_material[e]];
            const WaterRetention& w = retention[mesh.element_material[e]];
            const double ne = porosity[mesh.element_material[e]];
            const double ss = ss_sat[mesh.element_material[e]];
            Eigen::Matrix<double, N, 1> he, hpe;
            for (int i = 0; i < N; ++i) { he(i) = h[mesh.node_of(e, i)]; hpe(i) = hprev[mesh.node_of(e, i)]; }
            Eigen::Matrix<double, N, N> He = Eigen::Matrix<double, N, N>::Zero();
            Eigen::Matrix<double, N, N> Ce = Eigen::Matrix<double, N, N>::Zero();   // kapasitans C
            Eigen::Matrix<double, N, N> Sse = Eigen::Matrix<double, N, N>::Zero();  // S_s·S_e mass
            Eigen::Matrix<double, N, 1> bth = Eigen::Matrix<double, N, 1>::Zero();  // ∫N(θ−θ_prev)
            for (int g = 0; g < E::kGaussCount; ++g) {
                const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
                const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
                const double wdet = gp[g].weight * sd.det_jacobian;
                double h_gp = 0.0, hp_gp = 0.0, y_gp = 0.0;
                for (int i = 0; i < N; ++i) { h_gp += Nsh(i) * he(i); hp_gp += Nsh(i) * hpe(i); y_gp += Nsh(i) * X(i, 1); }
                const double psi = y_gp - h_gp, psip = y_gp - hp_gp;       // emme (current / prev)
                const double Se = effective_saturation(w, psi);
                const double kr = relative_permeability(w, Se);
                const double C = -ne * moisture_capacity(w, psi) + ss * Se;
                const double theta = ne * saturation(w, psi), theta_p = ne * saturation(w, psip);
                Eigen::Matrix<double, 1, N> gx, gy;
                for (int i = 0; i < N; ++i) { gx(0, i) = sd.B(0, 2 * i); gy(0, i) = sd.B(1, 2 * i + 1); }
                He.noalias()  += wdet * kr * (pm.kx * gx.transpose() * gx + pm.ky * gy.transpose() * gy);
                Ce.noalias()  += (wdet * C)       * (Nsh * Nsh.transpose());
                Sse.noalias() += (wdet * ss * Se) * (Nsh * Nsh.transpose());
                bth.noalias() += (wdet * (theta - theta_p)) * Nsh;
                water_vol += wdet * theta;
            }
            const Eigen::Matrix<double, N, 1> bss = Sse * (he - hpe);   // ∫Nᵀ S_s S_e N (h−h_prev)
            for (int a = 0; a < N; ++a) {
                const int na = mesh.node_of(e, a);
                storage_vec(na) += bth(a) + bss(a);
                for (int b = 0; b < N; ++b) {
                    const int nb = mesh.node_of(e, b);
                    Hb.add_entry(na, nb, He(a, b));
                    Ab.add_entry(na, nb, He(a, b) + Ce(a, b) / dt);
                }
            }
        }
        storage_vec /= dt;
    };

    Eigen::VectorXd h(nc);
    for (int n = 0; n < nc; ++n) h(n) = initial_head[n];
    auto nodal_sat = [&](const Eigen::VectorXd& hh) {
        Eigen::VectorXd S(nc);
        for (int n = 0; n < nc; ++n) S(n) = saturation(retention[node_mat[n]], mesh.y[n] - hh(n));
        return S;
    };

    UnsaturatedFlowResult R;
    {  // t=0 record
        math::SparseMatrixBuilder Ab(nc), Hb(nc);
        Ab.reserve((std::size_t)mesh.element_count * N * N); Hb.reserve((std::size_t)mesh.element_count * N * N);
        Eigen::VectorXd sv; double vol = 0.0;
        assemble(h, h, Ab, Hb, sv, vol);
        R.times.push_back(0.0); R.head.push_back(h); R.saturation.push_back(nodal_sat(h));
        R.water_volume.push_back(vol); R.darcy_influx.push_back(0.0);
    }

    for (int step = 0; step < nsteps; ++step) {
        const double t = (step + 1) * dt;
        for (int n = 0; n < nc; ++n) if (prescribed(n)) h(n) = head_bc(n, t);   // prescribed = g(t)
        const Eigen::VectorXd hn = h;   // h^n (prescribed values current)
        double water_vol = 0.0;
        Eigen::VectorXd Hh = Eigen::VectorXd::Zero(nc), storage_vec = Eigen::VectorXd::Zero(nc);
        bool ok = false; int it = 0;
        for (; it < max_picard; ++it) {
            math::SparseMatrixBuilder Ab(nc), Hb(nc);
            Ab.reserve((std::size_t)mesh.element_count * N * N); Hb.reserve((std::size_t)mesh.element_count * N * N);
            assemble(h, hn, Ab, Hb, storage_vec, water_vol);
            const math::CsrMatrix A = Ab.build();
            Hh = Hb.build() * h;
            math::SparseMatrixBuilder Affb(std::max(1, nf));   // extract A_ff (free^2)
            for (int r = 0; r < nc; ++r) {
                if (feq[r] < 0) continue;
                for (int idx = A.row_ptr[r]; idx < A.row_ptr[r + 1]; ++idx) {
                    const int c = A.col_indices[idx];
                    if (feq[c] >= 0) Affb.add_entry(feq[r], feq[c], A.values[idx]);
                }
            }
            const math::CsrMatrix Aff = Affb.build();
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(std::max(1, nf));
            for (int n = 0; n < nc; ++n)
                if (feq[n] >= 0) rhs(feq[n]) = qv(n) - Hh(n) - storage_vec(n);
            Eigen::VectorXd delta = Eigen::VectorXd::Zero(std::max(1, nf));
            if (nf > 0) {
                if (solve_factory) { delta = solve_factory(Aff)(rhs); }
                else {
                    Eigen::MatrixXd Ad = Eigen::MatrixXd::Zero(nf, nf);
                    for (int r = 0; r < nf; ++r)
                        for (int idx = Aff.row_ptr[r]; idx < Aff.row_ptr[r + 1]; ++idx)
                            Ad(r, Aff.col_indices[idx]) = Aff.values[idx];
                    delta = Eigen::PartialPivLU<Eigen::MatrixXd>(Ad).solve(rhs);
                }
            }
            double dmax = 0.0, scale = 1e-30;
            for (int n = 0; n < nc; ++n)
                if (feq[n] >= 0) { h(n) += delta(feq[n]); dmax = std::fmax(dmax, std::fabs(delta(feq[n]))); scale = std::fmax(scale, std::fabs(h(n))); }
            if (dmax < picard_tol * scale + 1e-12) { ok = true; ++it; break; }
        }
        R.max_picard_iters = std::max(R.max_picard_iters, it);
        if (!ok) R.converged = false;
        // Net boundary water discharge = reaction sum R = sum_{prescribed}(storage_i + (H.h)_i)
        // (pos = inflow). Global balance: sum R = sum_all storage = d(total store)/dt -> matches
        // int dtheta in the time integral (mass conservation).
        double infl = 0.0;
        for (int n = 0; n < nc; ++n) if (prescribed(n)) infl += storage_vec(n) + Hh(n);
        R.darcy_influx.push_back(infl);
        R.times.push_back(t); R.head.push_back(h); R.saturation.push_back(nodal_sat(h));
        R.water_volume.push_back(water_vol);
        if (!R.converged) break;
    }
    return R;
}
}  // namespace detail

UnsaturatedFlowResult solve_transient_unsaturated_flow(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& porosity, const std::vector<double>& ss_sat,
    const std::vector<WaterRetention>& retention, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head, double dt, int nsteps,
    const std::vector<double>* flux, const TransientFlowSolveFactory& solve_factory,
    int max_picard, double picard_tol) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::unsaturated_flow_impl<Tri15Element>(mesh, perm, porosity, ss_sat, retention,
            is_prescribed, head_bc, initial_head, flux, dt, nsteps, solve_factory, max_picard, picard_tol);
    return detail::unsaturated_flow_impl<Tri6Element>(mesh, perm, porosity, ss_sat, retention,
        is_prescribed, head_bc, initial_head, flux, dt, nsteps, solve_factory, max_picard, picard_tol);
}

}  // namespace katai::core
