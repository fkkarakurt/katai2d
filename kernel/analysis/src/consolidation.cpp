// The Biot consolidation bodies (LE + elastoplastic), compiled ONCE (section 5.2
// batch 3a): the element-templated detail implementations instantiate here for tri6
// and tri15, and every consumer sees declarations only.
#include <katai/math/solve_error.hpp>   // SingularSystem: a refused solve is an answer
#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/nonlinear_solver.hpp>   // kStagnationAccept

#include <functional>
#include <limits>

namespace katai::core {

namespace detail {


template <class E>
ConsolidationResult consolidation_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const std::vector<MaterialModel>& materials,
                                       const std::vector<MaterialProfile>& profile,
                                       const std::vector<Permeability>& perm,
                                       double gamma_w, const PoreFluidStiffness& kw_over_n,
                                       const std::vector<char>& drained_node,
                                       const std::vector<double>& initial_pore,
                                       double dt, int nsteps,
                                       const std::vector<char>& active,
                                       const Eigen::VectorXd* load_increment,
                                       const ConsolidationSolveFactory& solve_factory,
                                       const math::CsrMatrix* struct_k,
                                       const std::vector<int>* pore_tie) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();

    // Pore DOF numbering (drained node = fixed p=0 -> -1; tied seam nodes share one equation).
    std::vector<int> pore_eq;
    const int npore = number_pore_equations(mesh.node_count, drained_node, pore_tie, pore_eq);

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
            Se.noalias() += w * (1.0 / kw_over_n.at(e_mat)) * (Nsh * Nsh.transpose());
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

    // Structural elements' elastic stiffness into the SAME mechanical block (see the header).
    add_structural_block(Abuild, struct_k);

    const math::CsrMatrix A = Abuild.build();
    const math::CsrMatrix Hcsr = Hbuild.build();

    // Solve backend: a caller-provided factory (factor once -> repeated back-solve, e.g. PARDISO
    // mtype=-2) when given; otherwise a dense Eigen LU (MKL-free reference, small meshes only).
    // A refused linear solve is a property of THIS system, not a fatal error: the coupled tangent
    // is singular (commonly an insufficiently restrained model), and a verifying backend refuses to
    // return a vector that does not satisfy it rather than returning a meaningless one. The static
    // solver has caught exactly this since it began checking its answers (nonlinear_solver.cpp:
    // "Only SingularSystem is caught"); here it used to escape the call stack and ABORT the
    // process, which is the one answer a solver may never give. An empty result is this core's
    // "no result", and the phase strategy already reports it.
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> solve_step;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu;
    if (solve_factory) {
        try { solve_step = solve_factory(A); }
        catch (const math::SingularSystem&) { return {}; }
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
        // FIRST step, i.e. at t=0+ of the consolidation phase. With the storage
        // matrix S≈0 (near-incompressible pore fluid) the t=0+ response is the UNDRAINED one — the
        // load generates an excess pore pressure that the subsequent (Δf=0) steps then dissipate.
        if (step == 0 && load_increment) rhs.head(ndisp) = *load_increment;
        if (npore > 0) rhs.tail(npore) += dt * (Hcsr * p);   // continuity right-hand side dt.H.p_n
        Eigen::VectorXd d;
        try { d = solve_step(rhs); }
        catch (const math::SingularSystem&) { return {}; }
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
    const std::vector<Permeability>& perm, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory, int max_newton, double newton_tol,
    const math::CsrMatrix* struct_k, const std::vector<int>* pore_tie) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();
    const int ngp = E::kGaussCount;

    std::vector<int> pore_eq;
    const int npore = number_pore_equations(mesh.node_count, drained_node, pore_tie, pore_eq);
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
                        Eigen::VectorXd& f_int, TangentMode mode, double dt) {
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
                Se.noalias() += w * (1.0 / kw_over_n.at(e_mat)) * (Nsh * Nsh.transpose());
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
        // The structural elements are linear here, so their internal force is K_s.dv exactly --
        // no state, no return mapping. Adding it to BOTH the tangent and the internal force is
        // what makes residual(dv) = 0 hold at the solution rather than approximately.
        add_structural_block(Ab, struct_k);
        if (struct_k) f_int.noalias() += (*struct_k) * dv;
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

    // ONE TIME STEP of length h carrying the load increment df_step, from the committed state.
    // Monolithic coupled Newton with the robustness the static solver has:
    //
    //  * TWO CRITERIA, EACH ON ITS OWN SCALE. Equilibrium is a force balance and is measured against
    //    the forces in play -- the committed internal force, the load increment and the force of the
    //    pore pressure (L p) -- and continuity is a volume balance, measured against its own terms.
    //    They used to share one reference, the larger of |df| and |dt H p|, which in a dissipation
    //    step with no load is the continuity term alone: measured on an embankment over Hardening
    //    Soil sand, a step whose out-of-balance force was 1.7e-3 kN against 237 kN in play (7e-6)
    //    was held to 3e-8 kN, oscillated there for the whole iteration budget, and stopped the
    //    phase as "did not converge".
    //  * A BACKTRACKING LINE SEARCH on the equilibrium residual (the continuity block is linear,
    //    so any full step satisfies it).
    //  * STAGNATION ACCEPTANCE at kStagnationAccept of the force scale (K2D-A019), counted.
    // Writes the converged increment to (dv, dpv) and leaves the trial states at it.
    const auto one_step = [&](double h, const Eigen::VectorXd& df_step, Eigen::VectorXd& dv,
                              Eigen::VectorXd& dpv, bool& stagnant) -> bool {
        stagnant = false;
        Eigen::VectorXd Bbase;
        { math::SparseMatrixBuilder a0(NT), h0(std::max(1, npore)), l0(ndisp, std::max(1, npore));
          Eigen::VectorXd dvz = Eigen::VectorXd::Zero(ndisp);
          assemble(dvz, a0, h0, l0, Bbase, tmode, h); }   // trial == committed -> Bbase = f_int(committed)
        dv = Eigen::VectorXd::Zero(ndisp);
        dpv = Eigen::VectorXd::Zero(std::max(1, npore));
        struct Eval {
            math::CsrMatrix A;
            Eigen::VectorXd r;
            double ru = 0.0, rp = 0.0, ref_u = 0.0, ref_p = 0.0;
        };
        const auto evaluate = [&](const Eigen::VectorXd& dv_, const Eigen::VectorXd& dpv_) {
            Eval ev;
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore)), Lb(ndisp, std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            Eigen::VectorXd f_int;
            assemble(dv_, Ab, Hb, Lb, f_int, tmode, h);
            ev.A = Ab.build();
            // Newton RHS = -R: equilibrium R_u = (f_int - Bbase) + L dp - df; continuity from the
            // assembled bottom block (A x).tail = Lt dv - S* dp (S* = h H + S).
            ev.r = Eigen::VectorXd::Zero(NT);
            ev.r.head(ndisp) = df_step - (f_int - Bbase);
            double pore_force = 0.0;
            if (npore > 0) {
                const math::CsrMatrix Hc = Hb.build();
                const math::CsrMatrix Lc = Lb.build();
                ev.r.head(ndisp).noalias() -= Lc * dpv_;
                Eigen::VectorXd x(NT);
                x.head(ndisp) = dv_;
                x.tail(npore) = dpv_;
                const Eigen::VectorXd Ax = ev.A * x;
                const Eigen::VectorXd hHp = h * (Hc * p);
                ev.r.tail(npore) = hHp - Ax.tail(npore);
                ev.rp = ev.r.tail(npore).norm();
                // The volume balance's own terms: the flow in the step and the volume change the
                // displacement makes (L^T dv). The assembled bottom row alone is no scale -- it IS
                // the balance, zero to round-off once a step has converged.
                Eigen::VectorXd xv = Eigen::VectorXd::Zero(NT);
                xv.head(ndisp) = dv_;
                ev.ref_p = std::max(hHp.norm(), (ev.A * xv).tail(npore).norm());   // |L^T dv|
                pore_force = (Lc * (p + dpv_)).norm();
            }
            ev.ru = ev.r.head(ndisp).norm();
            ev.ref_u = df_step.norm() + Bbase.norm() + pore_force;
            return ev;
        };
        const auto converged = [&](const Eval& ev) {
            return ev.ru <= newton_tol * ev.ref_u + 1e-12 &&
                   (npore == 0 || ev.rp <= newton_tol * ev.ref_p + 1e-15);
        };
        // The line search compares the two balances together, each relative to its own scale:
        // a step with no load starts in force balance (ru = 0) and out of volume balance, and a
        // search on the force balance alone would reject every step that restores the volume.
        const auto worse = [&](const Eval& a, const Eval& b) {   // is a no better than b?
            const double su = std::max({a.ref_u, b.ref_u, 1e-300});
            const double sp = std::max({a.ref_p, b.ref_p, 1e-300});
            const double ma = a.ru / su + (npore > 0 ? a.rp / sp : 0.0);
            const double mb = b.ru / su + (npore > 0 ? b.rp / sp : 0.0);
            return ma >= mb;
        };
        const auto rp_ok = [&](const Eval& ev) {
            return npore == 0 || ev.rp <= newton_tol * ev.ref_p + 1e-15;
        };
        Eval cur = evaluate(dv, dpv);
        Eigen::VectorXd best_dv = dv, best_dpv = dpv;
        bool best_rp_ok = rp_ok(cur);
        double best_ru = best_rp_ok ? cur.ru : std::numeric_limits<double>::infinity();
        for (int it = 0; it < max_newton && !converged(cur); ++it) {
            if (!solve_factory) return false;
            // A refused linear solve ends THIS time step, exactly as a non-converged Newton
            // does -- the tangent is singular at this iterate, which on this path means the soil
            // has reached its capacity under the load being consolidated. Only SingularSystem is
            // caught: a malformed request or a broken backend raises SolveError and still
            // propagates, because turning one of those into "did not converge" would publish a
            // modelling answer for a bug.
            Eigen::VectorXd d;
            try {
                const auto solve = solve_factory(cur.A);
                d = solve(cur.r);
            } catch (const math::SingularSystem&) { break; }
            double alpha = 1.0;
            Eigen::VectorXd dv_n, dpv_n;
            Eval next;
            for (int ls = 0; ls < 6; ++ls) {
                dv_n = dv + alpha * d.head(ndisp);
                dpv_n = dpv;
                if (npore > 0) dpv_n += alpha * d.tail(npore);
                next = evaluate(dv_n, dpv_n);
                if (!worse(next, cur) || converged(next)) break;
                alpha *= 0.5;
            }
            dv = std::move(dv_n);
            dpv = std::move(dpv_n);
            cur = std::move(next);
            if (rp_ok(cur) && cur.ru < best_ru) {
                best_ru = cur.ru;
                best_dv = dv;
                best_dpv = dpv;
                best_rp_ok = true;
            }
        }
        if (converged(cur)) return true;
        if (!best_rp_ok || !(best_ru <= kStagnationAccept * cur.ref_u)) return false;
        dv = best_dv;
        dpv = best_dpv;
        (void)evaluate(dv, dpv);   // the trial states of the iterate being kept
        stagnant = true;
        return true;
    };

    // A step that fails is cut in two, its load increment split with it, down to 1/64 of the
    // step, as a static phase cuts its load step back. Measured on the same embankment: a Mohr-
    // Coulomb fill activated in the first step of a consolidation phase could not be carried in
    // one Newton solve of 40 iterations at any number of time steps, while the same fill in a
    // Plastic phase -- which cuts its load -- was. Every sub-step that converges is committed and
    // recorded, so a cut step adds its intermediate times to the series.
    int depth_used = 0;
    const std::function<bool(double, double, const Eigen::VectorXd&, int)> advance =
        [&](double t0, double h, const Eigen::VectorXd& df_step, int depth) -> bool {
        Eigen::VectorXd dv, dpv;
        bool stagnant = false;
        if (one_step(h, df_step, dv, dpv, stagnant)) {
            v += dv;
            if (npore > 0) p += dpv;
            committed = trial;
            if (stagnant) {
                ++R.stagnation_accepted;
            }
            record(t0 + h);
            return true;
        }
        trial = committed;
        if (depth >= 6) return false;
        depth_used = std::max(depth_used, depth + 1);
        const Eigen::VectorXd half = 0.5 * df_step;
        return advance(t0, 0.5 * h, half, depth + 1) && advance(t0 + 0.5 * h, 0.5 * h, half, depth + 1);
    };
    const Eigen::VectorXd no_load = Eigen::VectorXd::Zero(ndisp);
    for (int step = 0; step < nsteps && R.converged; ++step) {
        const Eigen::VectorXd& df = (step == 0 && load_increment) ? *load_increment : no_load;
        if (!advance(step * dt, dt, df, 0)) R.converged = false;
    }
    R.steps_cut = depth_used;
    R.committed = committed;
    return R;
}
}  // namespace detail

ConsolidationResult solve_consolidation(const mesh::Mesh& mesh, const DofMap& dofs,
                                               const std::vector<MaterialModel>& materials,
                                               const std::vector<Permeability>& perm,
                                               double gamma_w, const PoreFluidStiffness& kw_over_n,
                                               const std::vector<char>& drained_node,
                                               const std::vector<double>& initial_pore,
                                               double dt, int nsteps,
                                               const std::vector<char>& active,
                                               const Eigen::VectorXd* load_increment,
                                               const ConsolidationSolveFactory& solve_factory,
                                               const std::vector<MaterialProfile>& profile,
                                               const math::CsrMatrix* struct_k,
                                               const std::vector<int>* pore_tie) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_impl<Tri15Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                        kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                        active, load_increment, solve_factory, struct_k,
                                                        pore_tie);
    return detail::consolidation_impl<Tri6Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                   kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                   active, load_increment, solve_factory, struct_k,
                                                   pore_tie);
}

ConsolidationPlasticResult solve_consolidation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const std::vector<MaterialProfile>& profile, const math::CsrMatrix* struct_k,
    const std::vector<int>* pore_tie) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_plastic_impl<Tri15Element>(
            mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
            initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
            struct_k, pore_tie);
    return detail::consolidation_plastic_impl<Tri6Element>(
        mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
        initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
        struct_k, pore_tie);
}

}  // namespace katai::core
