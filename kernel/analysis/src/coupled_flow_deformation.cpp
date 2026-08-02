// The fully-coupled flow-deformation bodies (LE Picard + elastoplastic monolithic
// Newton-Picard), compiled ONCE (section 5.2 batch 3b): the element-templated detail
// implementations instantiate here for tri6 and tri15; consumers see declarations only.
#include <katai/analysis/coupled_flow_deformation.hpp>

namespace katai::core {

namespace detail {
template <class E>
CoupledFlowResult coupled_flow_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory, int max_picard, double picard_tol) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();

    std::vector<int> pore_eq(mesh.node_count, -1);
    int npore = 0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (!drained_node[n]) pore_eq[n] = npore++;
    const int NT = ndisp + npore;
    const auto gp = E::gauss_points();
    const Eigen::Vector3d mvec(1.0, 1.0, 0.0);

    // Assembly of the saddle A=[K L_chi; L_chi' -(dt.H_kr+S_st)] + H_kr (npore^2, for the
    // continuity RHS dt.H_kr.p_n), coefficients evaluated at the suction psi=-p/gamma_w
    // interpolated at the Gauss point from the FULL nodal pore field `pfull`.
    auto assemble = [&](const Eigen::VectorXd& pfull, math::SparseMatrixBuilder& Ab,
                        math::SparseMatrixBuilder& Hb) {
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            typename E::NodeCoords X;
            for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
            const int e_mat = mesh.element_material[e];
            const MaterialModel& mat = materials[e_mat];
            const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
            const Eigen::Matrix3d M = mat.elastic_plane_strain();   // uniform case: hoisted, as before
            const Permeability& pm = perm[e_mat];
            const WaterRetention& w = retention[e_mat];
            const double ne = porosity[e_mat];
            Eigen::Matrix<double, ND, ND> Ke = Eigen::Matrix<double, ND, ND>::Zero();
            Eigen::Matrix<double, ND, N> Le = Eigen::Matrix<double, ND, N>::Zero();
            Eigen::Matrix<double, N, N> He = Eigen::Matrix<double, N, N>::Zero();
            Eigen::Matrix<double, N, N> Se = Eigen::Matrix<double, N, N>::Zero();
            for (int g = 0; g < E::kGaussCount; ++g) {
                const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
                const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
                const double wdet = gp[g].weight * sd.det_jacobian;
                double p_gp = 0.0;
                for (int i = 0; i < N; ++i) p_gp += Nsh(i) * pfull[mesh.node_of(e, i)];
                // Suction: p is TENSION-POSITIVE => p>0 (water in tension) = suction =
                // unsaturated; p<=0 (compressive pore pressure, under load/water table) =
                // saturated (S_e=1). psi = +p/gamma_w.
                const double psi = p_gp / gamma_w;
                const double Se_eff = effective_saturation(w, psi); // Bishop chi = S_eff
                const double kr = relative_permeability(w, Se_eff);
                const double Sdeg = saturation(w, psi);             // degree of saturation
                const double dS_dpw = -ne * moisture_capacity(w, psi) / gamma_w;  // n.dS/dp_w >= 0
                const double Cp = Sdeg / kw_over_n + dS_dpw;        // storage coefficient (pos)
                Eigen::Matrix<double, 2, N> G;
                for (int i = 0; i < N; ++i) { G(0, i) = sd.B(0, 2 * i); G(1, i) = sd.B(1, 2 * i + 1); }
                // Depth-varying E' at the stress point (uniform() keeps the hoisted M -> bit-for-bit).
                Eigen::Matrix3d Mg = M;
                if (!prof.uniform()) {
                    double yg = 0.0;
                    for (int i = 0; i < N; ++i) yg += Nsh(i) * X(i, 1);
                    MaterialModel mgg = mat;
                    mgg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, yg);
                    Mg = mgg.elastic_plane_strain();
                }
                Ke.noalias() += wdet * sd.B.transpose() * Mg * sd.B;
                Le.noalias() += (wdet * Se_eff) * sd.B.transpose() * (mvec * Nsh.transpose());
                He.noalias() += wdet * kr * (G.row(0).transpose() * (pm.kx / gamma_w) * G.row(0)
                                          + G.row(1).transpose() * (pm.ky / gamma_w) * G.row(1));
                Se.noalias() += (wdet * Cp) * (Nsh * Nsh.transpose());
            }
            for (int a = 0; a < N; ++a) {
                const int na = mesh.node_of(e, a);
                for (int ca = 0; ca < 2; ++ca) {
                    const int ra = dofs.equation(dofs.global_dof(na, ca));
                    if (ra < 0) continue;
                    for (int b = 0; b < N; ++b) {
                        const int nb = mesh.node_of(e, b);
                        for (int cb = 0; cb < 2; ++cb) {
                            const int rb = dofs.equation(dofs.global_dof(nb, cb));
                            if (rb >= 0) Ab.add_entry(ra, rb, Ke(2 * a + ca, 2 * b + cb));
                        }
                        const int pb = pore_eq[nb];
                        if (pb >= 0) {
                            Ab.add_entry(ra, ndisp + pb, Le(2 * a + ca, b));
                            Ab.add_entry(ndisp + pb, ra, Le(2 * a + ca, b));
                        }
                    }
                }
            }
            for (int a = 0; a < N; ++a) {
                const int pa = pore_eq[mesh.node_of(e, a)];
                if (pa < 0) continue;
                for (int b = 0; b < N; ++b) {
                    const int pb = pore_eq[mesh.node_of(e, b)];
                    if (pb < 0) continue;
                    Ab.add_entry(ndisp + pa, ndisp + pb, -(dt * He(a, b) + Se(a, b)));
                    Hb.add_entry(pa, pb, He(a, b));
                }
            }
        }
    };

    auto to_full_pore = [&](const Eigen::VectorXd& p) {
        Eigen::VectorXd pf = Eigen::VectorXd::Zero(mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n)
            if (pore_eq[n] >= 0) pf(n) = p(pore_eq[n]);
            else pf(n) = 0.0;   // drained: p=0
        return pf;
    };

    Eigen::VectorXd v = Eigen::VectorXd::Zero(ndisp);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(std::max(1, npore));
    for (int n = 0; n < mesh.node_count; ++n)
        if (pore_eq[n] >= 0) p(pore_eq[n]) = initial_pore[n];

    CoupledFlowResult R;
    auto nodal_sat = [&](const Eigen::VectorXd& pf) {
        Eigen::VectorXd S(mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n) {
            const int m0 = mesh.element_material[0];  // node->material approximation (output only)
            // Degree of saturation S (matches the field label): S_res-based, not S_e (Bishop chi
            // comes separately from Se_eff).
            S(n) = saturation(retention[m0], pf(n) / gamma_w);
        }
        return S;
    };
    auto record = [&](double t) {
        R.series.times.push_back(t);
        R.series.displacement.push_back(expand_to_full(dofs, v));
        R.series.pore.push_back(to_full_pore(p));
        R.saturation.push_back(nodal_sat(to_full_pore(p)));
    };
    record(0.0);

    for (int step = 0; step < nsteps && R.converged; ++step) {
        Eigen::VectorXd v_new = v, p_new = p;
        bool ok = false; int it = 0;
        for (; it < max_picard; ++it) {
            // Picard: evaluate the coefficients at the latest estimate p_new; solve the
            // increment from (v^n, p^n).
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            assemble(to_full_pore(p_new), Ab, Hb);
            const math::CsrMatrix A = Ab.build();
            const math::CsrMatrix Hc = Hb.build();
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(NT);
            if (step == 0 && load_increment) rhs.head(ndisp) = *load_increment;
            if (npore > 0) rhs.tail(npore) += dt * (Hc * p);   // dt.H_kr.p_n (start-of-step pore)
            Eigen::VectorXd d;
            if (solve_factory) { d = solve_factory(A)(rhs); }
            else {
                Eigen::MatrixXd Ad = Eigen::MatrixXd::Zero(NT, NT);
                for (int r = 0; r < NT; ++r)
                    for (int idx = A.row_ptr[r]; idx < A.row_ptr[r + 1]; ++idx)
                        Ad(r, A.col_indices[idx]) = A.values[idx];
                d = Eigen::PartialPivLU<Eigen::MatrixXd>(Ad).solve(rhs);
            }
            const Eigen::VectorXd v_prev = v_new, p_prev = p_new;
            v_new = v + d.head(ndisp);
            if (npore > 0) p_new = p + d.tail(npore);
            // Convergence: the change between successive Picard estimates.
            double dchg = (v_new - v_prev).norm() + (npore > 0 ? (p_new - p_prev).norm() : 0.0);
            double scale = v_new.norm() + (npore > 0 ? p_new.norm() : 0.0) + 1e-30;
            if (dchg <= picard_tol * scale) { ok = true; ++it; break; }
        }
        R.max_picard_iters = std::max(R.max_picard_iters, it);
        if (!ok) { R.converged = false; break; }
        v = v_new; if (npore > 0) p = p_new;
        record((step + 1) * dt);
    }
    return R;
}
}  // namespace detail

CoupledFlowResult solve_coupled_flow_deformation(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory,
    int max_picard, double picard_tol,
    const std::vector<MaterialProfile>& profile) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::coupled_flow_impl<Tri15Element>(mesh, dofs, materials, profile, perm, retention,
            porosity, gamma_w, kw_over_n, drained_node, initial_pore, dt, nsteps, active, load_increment,
            solve_factory, max_picard, picard_tol);
    return detail::coupled_flow_impl<Tri6Element>(mesh, dofs, materials, profile, perm, retention,
        porosity, gamma_w, kw_over_n, drained_node, initial_pore, dt, nsteps, active, load_increment,
        solve_factory, max_picard, picard_tol);
}

namespace detail {
template <class E>
CoupledFlowPlasticResult coupled_flow_plastic_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();
    const int ngp = E::kGaussCount;

    std::vector<int> pore_eq(mesh.node_count, -1);
    int npore = 0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (!drained_node[n]) pore_eq[n] = npore++;
    const int NT = ndisp + npore;

    std::vector<GaussState> committed = initial_state;
    if ((int)committed.size() != mesh.element_count * ngp)
        committed.assign((size_t)mesh.element_count * ngp, GaussState{});
    std::vector<GaussState> trial = committed;

    const auto gp = E::gauss_points();
    const Eigen::Vector3d mvec(1.0, 1.0, 0.0);

    // Assemble A = [K_T L_χ; L_χᵀ −(ΔtH_kr+S_st)] + separate L_χ / H_kr (for residual/RHS) + f_int(trial),
    // given mechanical increment `dv` (drives plasticity) and the FULL nodal pore field `pfull` (drives
    // unsaturated coefficients, Picard-lagged). trial Gauss states are written. Mirrors consolidation_plastic's
    // assemble with the W3 unsaturated coefficients (Bishop χ=S_eff on L, Mualem k_rel on H, n·dS/dp_w in S).
    auto assemble = [&](const Eigen::VectorXd& dv, const Eigen::VectorXd& pfull,
                        math::SparseMatrixBuilder& Ab, math::SparseMatrixBuilder& Hb,
                        math::SparseMatrixBuilder& Lb, Eigen::VectorXd& f_int, TangentMode mode) {
        f_int.setZero(ndisp);
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            typename E::NodeCoords X;
            for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
            const int e_mat = mesh.element_material[e];
            const MaterialModel& mat = materials[e_mat];
            const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
            const Permeability& pm = perm[e_mat];
            const WaterRetention& w = retention[e_mat];
            const double ne = porosity[mesh.element_material[e]];
            Eigen::Matrix<double, ND, 1> du_e = Eigen::Matrix<double, ND, 1>::Zero();
            int eqd[ND];
            for (int a = 0; a < N; ++a)
                for (int c = 0; c < 2; ++c) {
                    const int eq = dofs.equation(dofs.global_dof(mesh.node_of(e, a), c));
                    eqd[2 * a + c] = eq;
                    if (eq >= 0) du_e(2 * a + c) = dv(eq);
                }
            Eigen::Matrix<double, ND, ND> Ke = Eigen::Matrix<double, ND, ND>::Zero();
            Eigen::Matrix<double, ND, 1> fe = Eigen::Matrix<double, ND, 1>::Zero();
            Eigen::Matrix<double, ND, N> Le = Eigen::Matrix<double, ND, N>::Zero();
            Eigen::Matrix<double, N, N> He = Eigen::Matrix<double, N, N>::Zero();
            Eigen::Matrix<double, N, N> Se = Eigen::Matrix<double, N, N>::Zero();
            for (int g = 0; g < ngp; ++g) {
                const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
                const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
                const double wdet = gp[g].weight * sd.det_jacobian;
                const Eigen::Vector3d deps = sd.B * du_e;
                GaussState& tr = trial[(size_t)e * ngp + g];
                Eigen::Matrix3d Dt;
                // E'(y) / c'(y) at the stress point (uniform() -> the element material, bit-for-bit).
                const MaterialModel* mpp = &mat;
                MaterialModel mgg;
                if (!prof.uniform()) {
                    double yg = 0.0;
                    for (int i = 0; i < N; ++i) yg += Nsh(i) * X(i, 1);
                    mgg = mat;
                    mgg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, yg);
                    mgg.cohesion = profile_at(mat.cohesion, prof.c_inc, prof.y_ref, yg);
                    mpp = &mgg;
                }
                // dt: this coupled flow-deformation time step [day] -- SoftSoilCreep sees the
                // creep + flow interaction through it; the other models do not read the parameter.
                integrate_point(*mpp, committed[(size_t)e * ngp + g], deps, tr, Dt, mode, dt);
                // Suction psi=+p/gamma_w (TENSION-POSITIVE: p>0 suction=unsaturated; p<=0
                // saturated S_e=1). Bishop chi=S_eff.
                double p_gp = 0.0;
                for (int i = 0; i < N; ++i) p_gp += Nsh(i) * pfull[mesh.node_of(e, i)];
                const double psi = p_gp / gamma_w;
                const double Se_eff = effective_saturation(w, psi);
                const double kr = relative_permeability(w, Se_eff);
                const double Sdeg = saturation(w, psi);
                const double dS_dpw = -ne * moisture_capacity(w, psi) / gamma_w;  // n.dS/dp_w >= 0
                const double Cp = Sdeg / kw_over_n + dS_dpw;                       // storage coefficient
                fe.noalias() += wdet * sd.B.transpose() * tr.stress;
                Ke.noalias() += wdet * sd.B.transpose() * Dt * sd.B;
                Eigen::Matrix<double, 2, N> G;
                for (int i = 0; i < N; ++i) { G(0, i) = sd.B(0, 2 * i); G(1, i) = sd.B(1, 2 * i + 1); }
                Le.noalias() += (wdet * Se_eff) * sd.B.transpose() * (mvec * Nsh.transpose());
                He.noalias() += wdet * kr * (G.row(0).transpose() * (pm.kx / gamma_w) * G.row(0)
                                          + G.row(1).transpose() * (pm.ky / gamma_w) * G.row(1));
                Se.noalias() += (wdet * Cp) * (Nsh * Nsh.transpose());
            }
            for (int a = 0; a < N; ++a)
                for (int ca = 0; ca < 2; ++ca) {
                    const int ra = eqd[2 * a + ca];
                    if (ra < 0) continue;
                    f_int(ra) += fe(2 * a + ca);
                    for (int b = 0; b < N; ++b) {
                        for (int cb = 0; cb < 2; ++cb) {
                            const int rb = eqd[2 * b + cb];
                            if (rb >= 0) Ab.add_entry(ra, rb, Ke(2 * a + ca, 2 * b + cb));
                        }
                        const int pb = pore_eq[mesh.node_of(e, b)];
                        if (pb >= 0) {
                            Ab.add_entry(ra, ndisp + pb, Le(2 * a + ca, b));
                            Ab.add_entry(ndisp + pb, ra, Le(2 * a + ca, b));
                            Lb.add_entry(ra, pb, Le(2 * a + ca, b));   // separate L_χ (residual L_χ·Δp)
                        }
                    }
                }
            for (int a = 0; a < N; ++a) {
                const int pa = pore_eq[mesh.node_of(e, a)];
                if (pa < 0) continue;
                for (int b = 0; b < N; ++b) {
                    const int pb = pore_eq[mesh.node_of(e, b)];
                    if (pb < 0) continue;
                    Ab.add_entry(ndisp + pa, ndisp + pb, -(dt * He(a, b) + Se(a, b)));
                    Hb.add_entry(pa, pb, He(a, b));
                }
            }
        }
    };

    auto to_full_pore = [&](const Eigen::VectorXd& pp) {
        Eigen::VectorXd pf = Eigen::VectorXd::Zero(mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n)
            pf(n) = (pore_eq[n] >= 0) ? pp(pore_eq[n]) : 0.0;   // drained: p=0
        return pf;
    };

    Eigen::VectorXd v = Eigen::VectorXd::Zero(ndisp);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(std::max(1, npore));
    if (!initial_pore.empty())
        for (int n = 0; n < mesh.node_count; ++n)
            if (pore_eq[n] >= 0) p(pore_eq[n]) = initial_pore[n];

    CoupledFlowPlasticResult R;
    auto nodal_sat = [&](const Eigen::VectorXd& pf) {
        Eigen::VectorXd S(mesh.node_count);
        const int m0 = mesh.element_material[0];   // node->material approximation (output only)
        for (int n = 0; n < mesh.node_count; ++n)  // degree of saturation S (matches the label), not S_e
            S(n) = saturation(retention[m0], pf(n) / gamma_w);
        return S;
    };
    auto record = [&](double t) {
        R.series.times.push_back(t);
        R.series.displacement.push_back(expand_to_full(dofs, v));
        const Eigen::VectorXd pf = to_full_pore(p);
        R.series.pore.push_back(pf);
        R.saturation.push_back(nodal_sat(pf));
    };
    record(0.0);

    // Tangent mode (same as consolidation_plastic): continuum when HS is present (robust+cheap),
    // otherwise consistent.
    bool any_hs = false;
    for (const auto& mm : materials) if (mm.type == MaterialType::HardeningSoil) any_hs = true;
    const TangentMode tmode = any_hs ? TangentMode::kContinuum : TangentMode::kConsistent;

    for (int step = 0; step < nsteps && R.converged; ++step) {
        // Base internal force of the committed (start-of-step) effective state (coefficients
        // evaluated at the committed pore).
        Eigen::VectorXd Bbase;
        { math::SparseMatrixBuilder a0(NT), h0(std::max(1, npore)), l0(ndisp, std::max(1, npore));
          assemble(Eigen::VectorXd::Zero(ndisp), to_full_pore(p), a0, h0, l0, Bbase, tmode); }
        Eigen::VectorXd dv = Eigen::VectorXd::Zero(ndisp), dpv = Eigen::VectorXd::Zero(std::max(1, npore));
        Eigen::VectorXd df = Eigen::VectorXd::Zero(ndisp);
        if (step == 0 && load_increment) df = *load_increment;
        bool step_ok = false; int it = 0;
        for (; it < max_newton; ++it) {
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore)), Lb(ndisp, std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            Eigen::VectorXd f_int;
            assemble(dv, to_full_pore(p + dpv), Ab, Hb, Lb, f_int, tmode);   // coefficients at p+dp (Picard lag)
            const math::CsrMatrix A = Ab.build();
            Eigen::VectorXd r = Eigen::VectorXd::Zero(NT);
            r.head(ndisp) = df - (f_int - Bbase);
            if (npore > 0) {
                const math::CsrMatrix Hc = Hb.build();
                const math::CsrMatrix Lc = Lb.build();
                r.head(ndisp).noalias() -= Lc * dpv;                            // − L_χ·Δp
                Eigen::VectorXd x(NT); x.head(ndisp) = dv; x.tail(npore) = dpv;
                const Eigen::VectorXd Ax = A * x;
                r.tail(npore) = dt * (Hc * p) - Ax.tail(npore);                 // ΔtH_kr·p_n − (L_χᵀΔv − S*Δp)
                const double ref = std::max({df.norm(), (dt * (Hc * p)).norm(), 1e-30});
                if (r.norm() <= newton_tol * ref + 1e-9 * (Bbase.norm() + 1.0)) { step_ok = true; break; }
            } else {
                const double ref = std::max(df.norm(), 1e-30);
                if (r.norm() <= newton_tol * ref + 1e-9 * (Bbase.norm() + 1.0)) { step_ok = true; break; }
            }
            if (!solve_factory) { step_ok = false; break; }
            const auto solve = solve_factory(A);
            const Eigen::VectorXd d = solve(r);
            dv += d.head(ndisp);
            if (npore > 0) dpv += d.tail(npore);
        }
        R.max_newton_iters = std::max(R.max_newton_iters, it);
        if (!step_ok) { R.converged = false; break; }
        v += dv; if (npore > 0) p += dpv; committed = trial;
        record((step + 1) * dt);
    }
    R.committed = committed;
    return R;
}
}  // namespace detail

CoupledFlowPlasticResult solve_coupled_flow_deformation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps,
    const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const std::vector<MaterialProfile>& profile) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::coupled_flow_plastic_impl<Tri15Element>(mesh, dofs, materials, profile, perm,
            retention, porosity, gamma_w, kw_over_n, drained_node, initial_state, initial_pore, dt,
            nsteps, active, load_increment, solve_factory, max_newton, newton_tol);
    return detail::coupled_flow_plastic_impl<Tri6Element>(mesh, dofs, materials, profile, perm,
        retention, porosity, gamma_w, kw_over_n, drained_node, initial_state, initial_pore, dt,
        nsteps, active, load_increment, solve_factory, max_newton, newton_tol);
}

}  // namespace katai::core
