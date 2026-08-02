// The Biot consolidation bodies (LE + elastoplastic), compiled ONCE (section 5.2
// batch 3a): the element-templated detail implementations instantiate here for tri6
// and tri15, and every consumer sees declarations only.
#include <katai/analysis/consolidation.hpp>

namespace katai::core {

namespace detail {
template <class E>
ConsolidationResult consolidation_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const std::vector<MaterialModel>& materials,
                                       const std::vector<MaterialProfile>& profile,
                                       const std::vector<Permeability>& perm,
                                       double gamma_w, double kw_over_n,
                                       const std::vector<char>& drained_node,
                                       const std::vector<double>& initial_pore,
                                       double dt, int nsteps,
                                       const std::vector<char>& active,
                                       const Eigen::VectorXd* load_increment,
                                       const ConsolidationSolveFactory& solve_factory) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();

    // Pore DOF numbering (drained node = fixed p=0 -> -1).
    std::vector<int> pore_eq(mesh.node_count, -1);
    int npore = 0;
    for (int n = 0; n < mesh.node_count; ++n)
        if (!drained_node[n]) pore_eq[n] = npore++;

    // The combined (saddle-point) system A = [K  L; L'  -(dt.H+S)] -- FULL, symmetric. Pore
    // equations are offset by ndisp. dt is fixed -> pattern+values built once, A factorized once.
    // H is also stored on its own (npore x npore): the continuity right-hand side dt.H.p_n needs
    // the product H.p at every step.
    const int NT = ndisp + npore;
    math::SparseMatrixBuilder Abuild(NT);
    math::SparseMatrixBuilder Hbuild(std::max(1, npore));
    Abuild.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
    Hbuild.reserve((std::size_t)mesh.element_count * N * N);

    const auto gp = E::gauss_points();
    const Eigen::Vector3d mvec(1.0, 1.0, 0.0);
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active.empty() && !active[e]) continue;   // staged: excavated / not-yet-placed soil
        typename E::NodeCoords X;
        for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        const int e_mat = mesh.element_material[e];
        const MaterialModel& mat = materials[e_mat];
        const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
        const Eigen::Matrix3d M = mat.elastic_plane_strain();   // uniform case: hoisted, as before
        const Permeability& pm = perm[e_mat];

        Eigen::Matrix<double, ND, ND> Ke = Eigen::Matrix<double, ND, ND>::Zero();
        Eigen::Matrix<double, ND, N> Le = Eigen::Matrix<double, ND, N>::Zero();
        Eigen::Matrix<double, N, N> He = Eigen::Matrix<double, N, N>::Zero();
        Eigen::Matrix<double, N, N> Se = Eigen::Matrix<double, N, N>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
            const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
            const double w = gp[g].weight * sd.det_jacobian;
            Eigen::Matrix<double, 2, N> G;
            for (int i = 0; i < N; ++i) { G(0, i) = sd.B(0, 2 * i); G(1, i) = sd.B(1, 2 * i + 1); }
            // Depth-varying E' at the stress point (uniform() keeps the hoisted M -> bit-for-bit).
            Eigen::Matrix3d Mg = M;
            if (!prof.uniform()) {
                double y = 0.0;
                for (int i = 0; i < N; ++i) y += Nsh(i) * X(i, 1);
                MaterialModel mg = mat;
                mg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, y);
                Mg = mg.elastic_plane_strain();
            }
            Ke.noalias() += w * sd.B.transpose() * Mg * sd.B;
            Le.noalias() += w * sd.B.transpose() * (mvec * Nsh.transpose());
            He.noalias() += w * (G.row(0).transpose() * (pm.kx / gamma_w) * G.row(0)
                               + G.row(1).transpose() * (pm.ky / gamma_w) * G.row(1));
            Se.noalias() += w * (1.0 / kw_over_n) * (Nsh * Nsh.transpose());
        }
        // Scatter into the combined system (full symmetric) + H (for the RHS product).
        for (int a = 0; a < N; ++a) {
            const int na = mesh.node_of(e, a);
            for (int ca = 0; ca < 2; ++ca) {
                const int ra = dofs.equation(dofs.global_dof(na, ca));
                if (ra < 0) continue;
                for (int b = 0; b < N; ++b) {
                    const int nb = mesh.node_of(e, b);
                    for (int cb = 0; cb < 2; ++cb) {
                        const int rb = dofs.equation(dofs.global_dof(nb, cb));
                        if (rb >= 0) Abuild.add_entry(ra, rb, Ke(2 * a + ca, 2 * b + cb));   // K
                    }
                    const int pb = pore_eq[nb];
                    if (pb >= 0) {                                                            // L + Lᵀ
                        Abuild.add_entry(ra, ndisp + pb, Le(2 * a + ca, b));
                        Abuild.add_entry(ndisp + pb, ra, Le(2 * a + ca, b));
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
                Abuild.add_entry(ndisp + pa, ndisp + pb, -(dt * He(a, b) + Se(a, b)));        // −(ΔtH+S)
                Hbuild.add_entry(pa, pb, He(a, b));
            }
        }
    }

    const math::CsrMatrix A = Abuild.build();
    const math::CsrMatrix Hcsr = Hbuild.build();

    // Solve backend: a caller-provided factory (factor once -> repeated back-solve, e.g. PARDISO
    // mtype=-2) when given; otherwise a dense Eigen LU (MKL-free reference, small meshes only).
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> solve_step;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu;
    if (solve_factory) {
        solve_step = solve_factory(A);
    } else {
        Eigen::MatrixXd Adense = Eigen::MatrixXd::Zero(NT, NT);
        for (int r = 0; r < NT; ++r)
            for (int idx = A.row_ptr[r]; idx < A.row_ptr[r + 1]; ++idx)
                Adense(r, A.col_indices[idx]) = A.values[idx];
        lu.compute(Adense);
        solve_step = [&lu](const Eigen::VectorXd& b) { return lu.solve(b); };
    }

    // Initial state: v=0, p=p0 (initial_pore at the free pore nodes).
    Eigen::VectorXd v = Eigen::VectorXd::Zero(ndisp);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(std::max(1, npore));
    p.setZero();
    for (int n = 0; n < mesh.node_count; ++n)
        if (pore_eq[n] >= 0) p(pore_eq[n]) = initial_pore[n];

    ConsolidationResult res;
    auto record = [&](double t) {
        res.times.push_back(t);
        res.displacement.push_back(expand_to_full(dofs, v));
        Eigen::VectorXd pfull = Eigen::VectorXd::Zero(mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n)
            if (pore_eq[n] >= 0) pfull(n) = p(pore_eq[n]);
        res.pore.push_back(pfull);
    };
    record(0.0);
    for (int step = 0; step < nsteps; ++step) {
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(NT);
        // Applied load increment Δf (staged-construction surcharge / fill): applied in full at the
        // FIRST step, i.e. at t=0+ of the consolidation phase (PLAXIS practice). With the storage
        // matrix S≈0 (near-incompressible pore fluid) the t=0+ response is the UNDRAINED one — the
        // load generates an excess pore pressure that the subsequent (Δf=0) steps then dissipate.
        if (step == 0 && load_increment) rhs.head(ndisp) = *load_increment;
        if (npore > 0) rhs.tail(npore) += dt * (Hcsr * p);   // continuity right-hand side dt.H.p_n
        const Eigen::VectorXd d = solve_step(rhs);
        v += d.head(ndisp);
        if (npore > 0) p += d.tail(npore);
        record((step + 1) * dt);
    }
    return res;
}
}  // namespace detail

namespace detail {
template <class E>
ConsolidationPlasticResult consolidation_plastic_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory, int max_newton, double newton_tol) {
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

    // Assemble the coupled system A = [K_T L; Lᵀ −(ΔtH+S)] (full symmetric for the disp/pore couple;
    // K_T itself is NONsymmetric for non-associated plasticity) + internal force f_int(trial), given the
    // current increment dv. trial Gauss states are written. Returns the sparse A + Hcsr (for the RHS).
    auto assemble = [&](const Eigen::VectorXd& dv, math::SparseMatrixBuilder& Ab,
                        math::SparseMatrixBuilder& Hb, math::SparseMatrixBuilder& Lb,
                        Eigen::VectorXd& f_int, TangentMode mode) {
        f_int.setZero(ndisp);
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            typename E::NodeCoords X;
            for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
            const int e_mat = mesh.element_material[e];
            const MaterialModel& mat = materials[e_mat];
            const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
            const Permeability& pm = perm[e_mat];
            // element translational increment du_e (fixed DOFs contribute 0)
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
                const double w = gp[g].weight * sd.det_jacobian;
                const Eigen::Vector3d deps = sd.B * du_e;
                GaussState& tr = trial[(size_t)e * ngp + g];
                Eigen::Matrix3d Dt;
                // E'(y) / c'(y) at the stress point (uniform() -> the element material, bit-for-bit).
                const MaterialModel* mp = &mat;
                MaterialModel mg;
                if (!prof.uniform()) {
                    double y = 0.0;
                    for (int i = 0; i < N; ++i) y += Nsh(i) * X(i, 1);
                    mg = mat;
                    mg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, y);
                    mg.cohesion = profile_at(mat.cohesion, prof.c_inc, prof.y_ref, y);
                    mp = &mg;
                }
                // dt: this consolidation time step [day] -- SoftSoilCreep sees the creep +
                // dissipation interaction through it; the other models do not read the parameter.
                integrate_point(*mp, committed[(size_t)e * ngp + g], deps, tr, Dt, mode, dt);
                fe.noalias() += w * sd.B.transpose() * tr.stress;
                Ke.noalias() += w * sd.B.transpose() * Dt * sd.B;
                Eigen::Matrix<double, 2, N> G;
                for (int i = 0; i < N; ++i) { G(0, i) = sd.B(0, 2 * i); G(1, i) = sd.B(1, 2 * i + 1); }
                Le.noalias() += w * sd.B.transpose() * (mvec * Nsh.transpose());
                He.noalias() += w * (G.row(0).transpose() * (pm.kx / gamma_w) * G.row(0)
                                   + G.row(1).transpose() * (pm.ky / gamma_w) * G.row(1));
                Se.noalias() += w * (1.0 / kw_over_n) * (Nsh * Nsh.transpose());
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
                            Lb.add_entry(ra, pb, Le(2 * a + ca, b));   // separate L (residual L·Δp)
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

    Eigen::VectorXd v = Eigen::VectorXd::Zero(ndisp);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(std::max(1, npore));
    if (!initial_pore.empty())   // optional initial excess pore (classic Terzaghi: dissipate u0, no load)
        for (int n = 0; n < mesh.node_count; ++n)
            if (pore_eq[n] >= 0) p(pore_eq[n]) = initial_pore[n];
    ConsolidationPlasticResult R;
    auto record = [&](double t) {
        R.series.times.push_back(t);
        R.series.displacement.push_back(expand_to_full(dofs, v));
        Eigen::VectorXd pf = Eigen::VectorXd::Zero(mesh.node_count);
        for (int n = 0; n < mesh.node_count; ++n) if (pore_eq[n] >= 0) pf(n) = p(pore_eq[n]);
        R.series.pore.push_back(pf);
    };
    record(0.0);

    // Per material, the integration tangent mode (MC: closed-form consistent; HS: continuum -> linear
    // but robust + far cheaper than the per-iterate numerical-consistent FD; LE: exact).
    bool any_hs = false;
    for (const auto& mm : materials) if (mm.type == MaterialType::HardeningSoil) any_hs = true;
    const TangentMode tmode = any_hs ? TangentMode::kContinuum : TangentMode::kConsistent;

    for (int step = 0; step < nsteps && R.converged; ++step) {
        // baseline internal force of the committed (start-of-step) effective state
        Eigen::VectorXd Bbase;
        { math::SparseMatrixBuilder a0(NT), h0(std::max(1, npore)), l0(ndisp, std::max(1, npore));
          Eigen::VectorXd dvz = Eigen::VectorXd::Zero(ndisp);
          assemble(dvz, a0, h0, l0, Bbase, tmode); }   // trial == committed here -> Bbase = f_int(committed)
        Eigen::VectorXd dv = Eigen::VectorXd::Zero(ndisp), dpv = Eigen::VectorXd::Zero(std::max(1, npore));
        Eigen::VectorXd df = Eigen::VectorXd::Zero(ndisp);
        if (step == 0 && load_increment) df = *load_increment;
        bool step_ok = false;
        for (int it = 0; it < max_newton; ++it) {
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore)), Lb(ndisp, std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            Eigen::VectorXd f_int;
            assemble(dv, Ab, Hb, Lb, f_int, tmode);
            const math::CsrMatrix A = Ab.build();
            // Newton RHS = −R: equilibrium R_u = (f_int − Bbase) + L·Δp − Δf; continuity R_p computed
            // from the assembled bottom block (A·x).tail = Lᵀ Δv − S* Δp (S* = ΔtH + S).
            Eigen::VectorXd r = Eigen::VectorXd::Zero(NT);
            r.head(ndisp) = df - (f_int - Bbase);
            if (npore > 0) {
                const math::CsrMatrix Hc = Hb.build();
                const math::CsrMatrix Lc = Lb.build();
                r.head(ndisp).noalias() -= Lc * dpv;                        // − L·Δp
                Eigen::VectorXd x(NT); x.head(ndisp) = dv; x.tail(npore) = dpv;
                const Eigen::VectorXd Ax = A * x;
                r.tail(npore) = dt * (Hc * p) - Ax.tail(npore);             // ΔtH·p_n − (Lᵀ Δv − S*Δp)
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
        if (!step_ok) { R.converged = false; break; }
        v += dv; if (npore > 0) p += dpv; committed = trial;
        record((step + 1) * dt);
    }
    R.committed = committed;
    return R;
}
}  // namespace detail

ConsolidationResult solve_consolidation(const mesh::Mesh& mesh, const DofMap& dofs,
                                               const std::vector<MaterialModel>& materials,
                                               const std::vector<Permeability>& perm,
                                               double gamma_w, double kw_over_n,
                                               const std::vector<char>& drained_node,
                                               const std::vector<double>& initial_pore,
                                               double dt, int nsteps,
                                               const std::vector<char>& active,
                                               const Eigen::VectorXd* load_increment,
                                               const ConsolidationSolveFactory& solve_factory,
                                               const std::vector<MaterialProfile>& profile) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_impl<Tri15Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                        kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                        active, load_increment, solve_factory);
    return detail::consolidation_impl<Tri6Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                   kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                   active, load_increment, solve_factory);
}

ConsolidationPlasticResult solve_consolidation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const std::vector<MaterialProfile>& profile) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_plastic_impl<Tri15Element>(
            mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
            initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol);
    return detail::consolidation_plastic_impl<Tri6Element>(
        mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
        initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol);
}

}  // namespace katai::core
