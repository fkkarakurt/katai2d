// The fully-coupled flow-deformation bodies (LE Picard + elastoplastic monolithic
// Newton-Picard), compiled ONCE (section 5.2 batch 3b): the element-templated detail
// implementations instantiate here for tri6 and tri15; consumers see declarations only.
#include <katai/analysis/cancel.hpp>
#include <katai/math/solve_error.hpp>   // SingularSystem: a refused solve is an answer
#include <katai/analysis/coupled_flow_deformation.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <katai/analysis/nonlinear_solver.hpp>   // StructureLaws
#include <katai/analysis/internal_forces.hpp>    // the kinematics policies (PlaneStrainKin, AxisymKin)

namespace katai::core {

namespace detail {
template <class E>
CoupledFlowResult coupled_flow_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory, int max_picard, double picard_tol,
    const math::CsrMatrix* struct_k, const std::vector<int>* pore_tie) {
    constexpr int N = E::kNodeCount;
    constexpr int ND = 2 * N;
    const int ndisp = dofs.equation_count();

    std::vector<int> pore_eq;
    const int npore = number_pore_equations(mesh.node_count, drained_node, pore_tie, pore_eq);
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
                const double Cp = Sdeg / kw_over_n.at(e_mat) + dS_dpw;   // storage coefficient (pos)
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
        katai::analysis::poll_cancel();
        Eigen::VectorXd v_new = v, p_new = p;
        bool ok = false; int it = 0;
        for (; it < max_picard; ++it) {
            katai::analysis::poll_cancel();
            // Picard: evaluate the coefficients at the latest estimate p_new; solve the
            // increment from (v^n, p^n).
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            assemble(to_full_pore(p_new), Ab, Hb);
            // The structural elements' elastic stiffness into the SAME mechanical block as the
            // soil's (consolidation.hpp: it arrives as a matrix, not as a system).
            add_structural_block(Ab, struct_k);
            const math::CsrMatrix A = Ab.build();
            const math::CsrMatrix Hc = Hb.build();
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(NT);
            if (step == 0 && load_increment) rhs.head(ndisp) = *load_increment;
            if (npore > 0) rhs.tail(npore) += dt * (Hc * p);   // dt.H_kr.p_n (start-of-step pore)
            Eigen::VectorXd d;
            // A refused solve ends the run rather than the process -- see consolidation.cpp, the
            // same discipline and the same reason (an empty result is this core's "no result").
            if (solve_factory) {
                try { d = solve_factory(A)(rhs); }
                catch (const math::SingularSystem&) { return {}; }
            }
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
    const std::vector<double>& porosity, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory,
    int max_picard, double picard_tol,
    const std::vector<MaterialProfile>& profile, const math::CsrMatrix* struct_k,
    const std::vector<int>* pore_tie) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::coupled_flow_impl<Tri15Element>(mesh, dofs, materials, profile, perm, retention,
            porosity, gamma_w, kw_over_n, drained_node, initial_pore, dt, nsteps, active, load_increment,
            solve_factory, max_picard, picard_tol, struct_k, pore_tie);
    return detail::coupled_flow_impl<Tri6Element>(mesh, dofs, materials, profile, perm, retention,
        porosity, gamma_w, kw_over_n, drained_node, initial_pore, dt, nsteps, active, load_increment,
        solve_factory, max_picard, picard_tol, struct_k, pore_tie);
}

namespace detail {
// Kin: the kinematics (PlaneStrainKin / AxisymKin, internal_forces.hpp), as in the consolidation
// twin -- the strain operator and its weight (per radian in axisymmetry), the material
// integration, the stress the element balances and the direction of the pore pressure (with the
// hoop in axisymmetry). The flow rows of B and the storage are the same under the same weight.
template <class E, class Kin = PlaneStrainKin>
CoupledFlowPlasticResult coupled_flow_plastic_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const math::CsrMatrix* struct_k, const std::vector<int>* pore_tie,
    StructureLaws* structure_laws) {
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
    const typename Kin::Strain mvec = Kin::pore_vector();

    // The length of the time step being solved: the phase's dt, or a part of it where a step that
    // did not close is cut (below).
    double h_step = dt;
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
                const auto sd = Kin::template gradients<E>(X, gp[g].xi, gp[g].eta);
                const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
                const double wdet = gp[g].weight * sd.weight;
                const typename Kin::Strain deps = sd.B * du_e;
                GaussState& tr = trial[(size_t)e * ngp + g];
                typename Kin::Tangent Dt;
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
                Kin::integrate(*mpp, committed[(size_t)e * ngp + g], deps, tr, Dt, mode, h_step);
                // Suction psi=+p/gamma_w (TENSION-POSITIVE: p>0 suction=unsaturated; p<=0
                // saturated S_e=1). Bishop chi=S_eff.
                double p_gp = 0.0;
                for (int i = 0; i < N; ++i) p_gp += Nsh(i) * pfull[mesh.node_of(e, i)];
                const double psi = p_gp / gamma_w;
                const double Se_eff = effective_saturation(w, psi);
                // The volumetric strain, of which the void ratio is read, and the CHANGE OF PERMEABILITY,
                // log10(k / k0) = de / ck, from the state at the start of the time step, as in the consolidation twin.
                tr.eps_vol = committed[(size_t)e * ngp + g].eps_vol + mvec.dot(deps);
                double kfac = 1.0;
                if (pm.ck > 0.0)
                    kfac = std::pow(10.0, (1.0 + mat.e_init) * std::expm1(committed[(size_t)e * ngp + g].eps_vol) / pm.ck);
                const double kr = kfac * relative_permeability(w, Se_eff);
                const double Sdeg = saturation(w, psi);
                const double dS_dpw = -ne * moisture_capacity(w, psi) / gamma_w;  // n.dS/dp_w >= 0
                const double Cp = Sdeg / kw_over_n.at(e_mat) + dS_dpw;             // storage coefficient
                fe.noalias() += wdet * sd.B.transpose() * Kin::stress(tr);
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
                    Ab.add_entry(ndisp + pa, ndisp + pb, -(h_step * He(a, b) + Se(a, b)));
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

    // The force a step leaves out of balance is carried into the next one, as in the consolidation
    // twin (see consolidation.cpp): dropped, it is lost once per time step.
    Eigen::VectorXd oob = Eigen::VectorXd::Zero(ndisp), oob_next = Eigen::VectorXd::Zero(ndisp);
    // ONE TIME STEP of length h_step carrying the load df, from the committed state; leaves the
    // trial states at the converged iterate.
    const auto one_step = [&](const Eigen::VectorXd& df, Eigen::VectorXd& dv, Eigen::VectorXd& dpv) -> bool {
        // Base internal force of the committed (start-of-step) effective state (coefficients
        // evaluated at the committed pore). No time passes in it: see the consolidation twin.
        Eigen::VectorXd Bbase;
        { math::SparseMatrixBuilder a0(NT), h0(std::max(1, npore)), l0(ndisp, std::max(1, npore));
          const double h_of_step = h_step;
          h_step = 0.0;
          assemble(Eigen::VectorXd::Zero(ndisp), to_full_pore(p), a0, h0, l0, Bbase, tmode);
          h_step = h_of_step; }
        // The structures by their own laws: the pore pressure a joint reads beside it is the
        // step's start, in the share the skeleton feels (Bishop chi); their force is in the base.
        if (structure_laws) {
            const Eigen::VectorXd pf = to_full_pore(p);
            structure_laws->set_soil_pore(std::vector<double>(pf.data(), pf.data() + pf.size()),
                                          [&](int mat, double p_gp) {
                                              return effective_saturation(retention[(size_t)mat], p_gp / gamma_w);
                                          });
            structure_laws->add(v, Bbase, nullptr);
        }
        dv = Eigen::VectorXd::Zero(ndisp);
        dpv = Eigen::VectorXd::Zero(std::max(1, npore));
        // The balance of an iterate, the matrix at it, and whether it has closed. The trial states
        // are left at the iterate evaluated last.
        // Converged: the force balance and the volume balance, each against its own scale, as in the consolidation
        // twin. One norm over the two, with a floor tied to the force the soil already carries, let a step close at
        // its first iterate once the flow of the step fell under that floor, with nothing moved: the pressure and
        // the settlement froze and the phase reported itself solved. Measured on a column of Mohr-Coulomb sand seated
        // at 1000 kPa under 10 kPa: 18% short of the end settlement, and of a jointed rock 75%, with the excess pressure
        // stopped at 2.9 and 9.6 kPa where the consolidation solver reaches 0.1; more of it on more steps. `tight` is
        // the first correction's test: it closes a step only well inside the tolerance (see that twin).
        struct Ev { Eigen::VectorXd r; math::CsrMatrix A; bool conv = false, tight = false; };
        const auto evaluate = [&](const Eigen::VectorXd& dv_, const Eigen::VectorXd& dpv_) {
            Ev ev;
            math::SparseMatrixBuilder Ab(NT), Hb(std::max(1, npore)), Lb(ndisp, std::max(1, npore));
            Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
            Eigen::VectorXd f_int;
            assemble(dv_, to_full_pore(p + dpv_), Ab, Hb, Lb, f_int, tmode);   // coefficients at p+dp (Picard lag)
            // Linear structures: their internal force is K_s.dv exactly, so it goes into BOTH the
            // tangent and f_int -- the same argument as the consolidation twin.
            add_structural_block(Ab, struct_k);
            if (struct_k) f_int.noalias() += (*struct_k) * dv_;
            if (structure_laws) structure_laws->add(v + dv_, f_int, &Ab);   // force and tangent at v + dv
            ev.A = Ab.build();
            ev.r = Eigen::VectorXd::Zero(NT);
            ev.r.head(ndisp) = df + oob - (f_int - Bbase);
            if (npore > 0) {
                const math::CsrMatrix Hc = Hb.build();
                const math::CsrMatrix Lc = Lb.build();
                ev.r.head(ndisp).noalias() -= Lc * dpv_;                        // − L_χ·Δp
                Eigen::VectorXd x(NT); x.head(ndisp) = dv_; x.tail(npore) = dpv_;
                const Eigen::VectorXd Ax = ev.A * x;
                const Eigen::VectorXd hHp = h_step * (Hc * p);
                ev.r.tail(npore) = hHp - Ax.tail(npore);                       // ΔtH_kr·p_n − (L_χᵀΔv − S*Δp)
                Eigen::VectorXd xv = Eigen::VectorXd::Zero(NT);
                xv.head(ndisp) = dv_;
                const Eigen::VectorXd Ltdv = (ev.A * xv).tail(npore);          // L_χᵀΔv
                const double ru = ev.r.head(ndisp).norm(), rp = ev.r.tail(npore).norm();
                const double ref_u = (df + oob).norm() + Bbase.norm() + (Lc * (p + dpv_)).norm();
                const double ref_p = std::max(hHp.norm(), Ltdv.norm());
                ev.conv = ru <= newton_tol * ref_u + 1e-12 && rp <= newton_tol * ref_p + 1e-15;
                ev.tight = ru <= 0.1 * newton_tol * ref_u + 1e-12 && rp <= 0.1 * newton_tol * ref_p + 1e-15;
            } else {
                const double ru = ev.r.norm(), ref_u = (df + oob).norm() + Bbase.norm();
                ev.conv = ru <= newton_tol * ref_u + 1e-12;
                ev.tight = ru <= 0.1 * newton_tol * ref_u + 1e-12;
            }
            return ev;
        };
        bool step_ok = false; int it = 0;
        Ev cur = evaluate(dv, dpv);
        std::vector<double> recent{cur.r.norm()};   // the balance at the last few iterates
        for (; it < max_newton; ++it) {
            katai::analysis::poll_cancel();
            if (cur.conv && (it >= 2 || cur.tight)) { step_ok = true; oob_next = cur.r.head(ndisp); break; }
            if (std::getenv("KATAI_NL_DEBUG"))
                std::fprintf(stderr, "    coupled h=%.4g it %d: |r| %.3e (disp %.3e) |Bbase| %.3e |df| %.3e\n", h_step, it,
                             cur.r.norm(), cur.r.head(ndisp).norm(), Bbase.norm(), df.norm());
            if (!solve_factory) { step_ok = false; break; }
            // A refused solve ends THIS time step, not the process (see consolidation.cpp).
            Eigen::VectorXd d;
            try {
                const auto solve = solve_factory(cur.A);
                d = solve(cur.r);
            } catch (const math::SingularSystem&) { step_ok = false; break; }
            // A NON-MONOTONE LINE SEARCH (Grippo, Lampariello & Lucidi 1986, memory four). The full
            // Newton step is taken where it leaves the balance below the worst of the last four
            // iterates -- every step of a smooth problem, so such a phase iterates exactly as it did;
            // otherwise the step is halved until it does. Where no halving does, the full step is
            // taken after all: a joint slipping in a dissipation step is closed by full steps, and
            // taking the shortest one instead crawled away from the answer, the balance growing at
            // every iterate. Measured on a raft under twenty times its plastic moment: with full steps
            // only, the iterates cycled between two states with the balance at 25% of the load, on
            // every cut of the time step.
            const double bar = *std::max_element(recent.begin(), recent.end());
            double alpha = 1.0;
            Ev nxt;
            bool found = false;
            for (int ls = 0; ls < 8; ++ls) {
                katai::analysis::poll_cancel();
                Eigen::VectorXd dv2 = dv + alpha * d.head(ndisp);
                Eigen::VectorXd dpv2 = dpv;
                if (npore > 0) dpv2 += alpha * d.tail(npore);
                nxt = evaluate(dv2, dpv2);
                if (nxt.conv || nxt.r.norm() < bar) { found = true; break; }
                alpha *= 0.5;
            }
            if (!found) {   // the full step, evaluated again so the trial states are at it
                alpha = 1.0;
                Eigen::VectorXd dv2 = dv + d.head(ndisp);
                Eigen::VectorXd dpv2 = dpv;
                if (npore > 0) dpv2 += d.tail(npore);
                nxt = evaluate(dv2, dpv2);
            }
            dv += alpha * d.head(ndisp);
            if (npore > 0) dpv += alpha * d.tail(npore);
            cur = std::move(nxt);
            recent.push_back(cur.r.norm());
            if (recent.size() > 4) recent.erase(recent.begin());
        }
        if (!step_ok && cur.conv && it == max_newton) { step_ok = true; oob_next = cur.r.head(ndisp); }
        R.max_newton_iters = std::max(R.max_newton_iters, it);
        return step_ok;
    };
    // A STEP THAT DOES NOT CLOSE IS CUT IN TWO, its load split with it, down to 1/64 of the step,
    // as the consolidation twin does (consolidation.cpp, `advance`). A plate reaching its plastic
    // moment, an anchor its capacity or a joint its Coulomb limit makes the tangent jump, and a
    // whole time step taken across the jump did not converge -- the phase used to stop there.
    // Every sub-step that converges is committed and recorded; a step that closes at once is
    // solved and recorded exactly as before.
    const std::function<bool(double, double, double, const Eigen::VectorXd&, int)> advance =
        [&](double t0, double h, double t1, const Eigen::VectorXd& df, int depth) -> bool {
        h_step = h;
        Eigen::VectorXd dv, dpv;
        if (one_step(df, dv, dpv)) {
            v += dv; if (npore > 0) p += dpv; committed = trial;
            if (structure_laws) structure_laws->commit();   // evaluated last at the converged dv
            oob = oob_next;
            record(t1);
            return true;
        }
        trial = committed;
        if (depth >= 6) return false;
        const Eigen::VectorXd half = 0.5 * df;
        return advance(t0, 0.5 * h, t0 + 0.5 * h, half, depth + 1) &&
               advance(t0 + 0.5 * h, 0.5 * h, t1, half, depth + 1);
    };
    for (int step = 0; step < nsteps && R.converged; ++step) {
        katai::analysis::poll_cancel();
        Eigen::VectorXd df = Eigen::VectorXd::Zero(ndisp);
        if (step == 0 && load_increment) df = *load_increment;
        if (!advance(step * dt, dt, (step + 1) * dt, df, 0)) R.converged = false;
    }
    R.committed = committed;
    return R;
}
}  // namespace detail

CoupledFlowPlasticResult solve_coupled_flow_deformation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps,
    const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const std::vector<MaterialProfile>& profile, const math::CsrMatrix* struct_k,
    const std::vector<int>* pore_tie, StructureLaws* structure_laws, bool axisymmetric) {
    if (axisymmetric) {
        if (mesh.nodes_per_element == Tri15Element::kNodeCount)
            return detail::coupled_flow_plastic_impl<Tri15Element, detail::AxisymKin>(mesh, dofs, materials,
                profile, perm, retention, porosity, gamma_w, kw_over_n, drained_node, initial_state,
                initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton,
                newton_tol, struct_k, pore_tie, structure_laws);
        return detail::coupled_flow_plastic_impl<Tri6Element, detail::AxisymKin>(mesh, dofs, materials,
            profile, perm, retention, porosity, gamma_w, kw_over_n, drained_node, initial_state,
            initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
            struct_k, pore_tie, structure_laws);
    }
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::coupled_flow_plastic_impl<Tri15Element>(mesh, dofs, materials, profile, perm,
            retention, porosity, gamma_w, kw_over_n, drained_node, initial_state, initial_pore, dt,
            nsteps, active, load_increment, solve_factory, max_newton, newton_tol, struct_k, pore_tie,
            structure_laws);
    return detail::coupled_flow_plastic_impl<Tri6Element>(mesh, dofs, materials, profile, perm,
        retention, porosity, gamma_w, kw_over_n, drained_node, initial_state, initial_pore, dt,
        nsteps, active, load_increment, solve_factory, max_newton, newton_tol, struct_k, pore_tie,
        structure_laws);
}

}  // namespace katai::core
