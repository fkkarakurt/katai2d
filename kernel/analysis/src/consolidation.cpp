// The Biot consolidation bodies (LE + elastoplastic), compiled ONCE (section 5.2
// batch 3a): the element-templated detail implementations instantiate here for tri6
// and tri15, and every consumer sees declarations only.
#include <katai/math/solve_error.hpp>   // SingularSystem: a refused solve is an answer
#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/nonlinear_solver.hpp>   // kStagnationAccept
#include <katai/analysis/internal_forces.hpp>    // the kinematics policies (PlaneStrainKin, AxisymKin)
#include <katai/math/thread_pool.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
                                       const std::vector<int>* pore_tie, bool ramp_load,
    const Eigen::VectorXd* instant_load) {
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
        // ... or, `ramp_load`, in equal parts over the steps: the load of a phase that has a
        // duration is put on DURING that time, linearly, and what the phase ends with is then
        // the pore pressure the loading generated less what drained while it was being placed.
        if (load_increment && ramp_load) rhs.head(ndisp) = *load_increment / (double)nsteps;
        else if (step == 0 && load_increment) rhs.head(ndisp) = *load_increment;
        // What the phase INHERITS is not placed over time: the pore pressure the phase before it
        // left is a state, re-generated in this phase's own pore field at its first instant.
        if (step == 0 && instant_load) rhs.head(ndisp) += *instant_load;
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
// Kin: the kinematics (PlaneStrainKin / AxisymKin, internal_forces.hpp) -- the strain operator and
// its integration weight (r-weighted, per radian, in axisymmetry), the material integration, the
// stress the element balances and the direction of the pore pressure (with the hoop in
// axisymmetry). The flow rows of B (dN/dx, dN/dy) and the storage are the same in both, under the
// same weight.
template <class E, class Kin = PlaneStrainKin>
ConsolidationPlasticResult consolidation_plastic_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<MaterialProfile>& profile,
    const std::vector<Permeability>& perm, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active, const Eigen::VectorXd* load_increment,
    const ConsolidationSolveFactory& solve_factory, int max_newton, double newton_tol,
    const math::CsrMatrix* struct_k, const std::vector<int>* pore_tie, bool ramp_load,
    const Eigen::VectorXd* instant_load, StructureLaws* structure_laws) {
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

    // Assemble the coupled system A = [K_T L; Lᵀ −(ΔtH+S)] (full symmetric for the disp/pore couple;
    // K_T itself is NONsymmetric for non-associated plasticity) + internal force f_int(trial), given the
    // current increment dv. trial Gauss states are written. Returns the sparse A + Hcsr (for the RHS).
    // A FRACTION OF THE ELASTIC OPERATOR ADDED TO THE TANGENT OF EVERY PLASTIC POINT, off until
    // a linear solve is refused or returns something that is not a number. A point returned to
    // a vertex hands back a tangent of zero and its neighbours a rank-deficient one; a node
    // surrounded by them makes the matrix singular with no mechanism formed. The static solver
    // has retried such an increment with this floor since it met the case (kVertexFloor); a
    // consolidation step did not, and on a finer mesh of an embankment over soft clay the
    // second iterate of the first time step came back as NaN at every step length. Only the
    // tangent changes: the internal force, and so the balance a step is driven to, is the same.
    double tangent_floor = 0.0;
    // THE ELASTIC STIFFNESS, for a step whose tangent has stopped giving a direction of descent.
    // Every plastic point hands the matrix its elastic operator instead of its plastic one: the
    // initial-stiffness iteration, slow and sure, that a Plastic phase falls back on in the same
    // situation. Measured on an embankment over soft clay on a 1.2 m mesh: the balance of a time
    // step sat at 3.4e-3 of its scale while the line search ended at 1e-4 of the step, iteration
    // after iteration, on each half of the step down to a sixty-fourth, and the phase stopped
    // there. As with the floor, only the matrix changes.
    bool elastic_matrix = false;
    // THE SMALL-STRAIN STIFFNESS OF A TIME STEP, READ FROM A FIXED INCREMENT: the same device,
    // for the same reason, as in the static solver (see small_freeze in nonlinear_solver.cpp).
    // The stiffness of a step of Hardening Soil with small-strain stiffness follows the direction
    // of the step's strain increment; a time step that stalls on it is continued with the
    // stiffness read from the displacement increment it had reached, and that increment is
    // replaced by the answer and the step closed again, twice.
    bool has_small = false;
    for (const auto& mm : materials)
        if (mm.type == MaterialType::HardeningSoil && mm.hs.G0_ref > 0.0) has_small = true;
    bool small_frozen = false;
    Eigen::VectorXd dv_small;
    // TWO-PHASE ASSEMBLY, as in the static solver: the element matrices -- which is where the
    // stress integration is, and with it nearly all of the time of a consolidation phase -- are
    // computed in parallel, each element writing only its own block and its own stress points;
    // they are then added to the system one element after the other, in element order, so the
    // sums are the same whatever the number of threads. Measured on fill placed on soft clay
    // with 15-noded elements: 449 of the 473 s of the phase were in this loop, on one core.
    // Developer timers of the phase (KATAI_NL_PROFILE): where its wall-clock time goes.
    double prof_par = 0.0;
    struct ElemOut {
        Eigen::Matrix<double, ND, ND> Ke;
        Eigen::Matrix<double, ND, 1> fe;
        Eigen::Matrix<double, ND, N> Le;
        Eigen::Matrix<double, N, N> He, Se;
        int eqd[ND];
    };
    std::vector<ElemOut> elem_out;
    // The structures by their own laws read the TOTAL displacement of the phase, v + dv: v is the
    // committed one (declared below, with the time march), and the increment of the evaluation
    // the matrix is built from is kept for the matrix.
    const Eigen::VectorXd* v_phase = nullptr;
    Eigen::VectorXd struct_dv;
    auto assemble = [&](const Eigen::VectorXd& dv, Eigen::VectorXd& f_int, TangentMode mode,
                        double dt) {
        f_int.setZero(ndisp);
        elem_out.resize((size_t)mesh.element_count);
        const auto pp0 = std::chrono::steady_clock::now();
        math::parallel_for(mesh.element_count, [&](int e_begin, int e_end) {
        for (int e = e_begin; e < e_end; ++e) {
            if (!active.empty() && !active[e]) continue;
            typename E::NodeCoords X;
            for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
            const int e_mat = mesh.element_material[e];
            const MaterialModel& mat = materials[e_mat];
            const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
            const Permeability& pm = perm[e_mat];
            // element translational increment du_e (fixed DOFs contribute 0)
            Eigen::Matrix<double, ND, 1> du_e = Eigen::Matrix<double, ND, 1>::Zero();
            int* eqd = elem_out[(size_t)e].eqd;
            for (int a = 0; a < N; ++a)
                for (int c = 0; c < 2; ++c) {
                    const int eq = dofs.equation(dofs.global_dof(mesh.node_of(e, a), c));
                    eqd[2 * a + c] = eq;
                    if (eq >= 0) du_e(2 * a + c) = dv(eq);
                }
            const bool read_small = small_frozen && dv_small.size() == dv.size();
            Eigen::Matrix<double, ND, 1> du_e_small = Eigen::Matrix<double, ND, 1>::Zero();
            if (read_small)
                for (int a = 0; a < ND; ++a)
                    if (eqd[a] >= 0) du_e_small(a) = dv_small(eqd[a]);
            ElemOut& eo_w = elem_out[(size_t)e];
            auto& Ke = eo_w.Ke; auto& fe = eo_w.fe; auto& Le = eo_w.Le; auto& He = eo_w.He; auto& Se = eo_w.Se;
            Ke.setZero(); fe.setZero(); Le.setZero(); He.setZero(); Se.setZero();
            for (int g = 0; g < ngp; ++g) {
                const auto sd = Kin::template gradients<E>(X, gp[g].xi, gp[g].eta);
                const auto Nsh = E::shape_functions(gp[g].xi, gp[g].eta);
                const double w = gp[g].weight * sd.weight;
                const typename Kin::Strain deps = sd.B * du_e;
                GaussState& tr = trial[(size_t)e * ngp + g];
                typename Kin::Tangent Dt;
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
                typename Kin::Report rep;
                typename Kin::Strain deps_small;
                if (read_small) deps_small = sd.B * du_e_small;
                Kin::integrate(*mp, committed[(size_t)e * ngp + g], deps, tr, Dt, mode, dt,
                               (tangent_floor > 0.0 || elastic_matrix) ? &rep : nullptr, 0.0,
                               read_small ? &deps_small : nullptr);
                if (elastic_matrix && rep.plastic) Dt = rep.elastic;
                else if (tangent_floor > 0.0 && rep.plastic) Dt += tangent_floor * rep.elastic;
                // The volumetric strain is accumulated here as in every other phase: the void
                // ratio is read from it (the dilatancy cut-off, the change of permeability), and
                // a consolidation phase used to leave it where the phase before had left it.
                tr.eps_vol = committed[(size_t)e * ngp + g].eps_vol + mvec.dot(deps);   // m.deps: with the hoop in axisymmetry
                // CHANGE OF PERMEABILITY: log10(k / k0) = de / ck, de the change in void ratio
                // since the initial state. It is read from the COMMITTED state -- the void ratio
                // at the start of the time step -- so the flow block stays linear within the
                // step and the iteration is what it was.
                double kfac = 1.0;
                if (pm.ck > 0.0) {
                    const double de =
                        (1.0 + mat.e_init) * std::expm1(committed[(size_t)e * ngp + g].eps_vol);
                    kfac = std::pow(10.0, de / pm.ck);
                }
                if (std::getenv("KATAI_NL_DEBUG") && (!Dt.allFinite() || !tr.stress.allFinite())) {
                    const GaussState& c0 = committed[(size_t)e * ngp + g];
                    std::fprintf(stderr, "    NONFINITE el %d gp %d mat %d type %d  Dt finite %d  stress finite %d  "
                                 "committed s=(%.4g %.4g %.4g) zz=%.4g pp=%.4g  deps=(%.3e %.3e %.3e)\n",
                                 e, g, e_mat, (int)mat.type, (int)Dt.allFinite(), (int)tr.stress.allFinite(),
                                 c0.stress(0), c0.stress(1), c0.stress(2), c0.stress_zz, c0.pp,
                                 deps(0), deps(1), deps(2));
                }
                fe.noalias() += w * sd.B.transpose() * Kin::stress(tr);
                Ke.noalias() += w * sd.B.transpose() * Dt * sd.B;
                Eigen::Matrix<double, 2, N> G;
                for (int i = 0; i < N; ++i) { G(0, i) = sd.B(0, 2 * i); G(1, i) = sd.B(1, 2 * i + 1); }
                Le.noalias() += w * sd.B.transpose() * (mvec * Nsh.transpose());
                He.noalias() += w * kfac * (G.row(0).transpose() * (pm.kx / gamma_w) * G.row(0)
                                   + G.row(1).transpose() * (pm.ky / gamma_w) * G.row(1));
                Se.noalias() += w * (1.0 / kw_over_n.at(e_mat)) * (Nsh * Nsh.transpose());
            }
        }
        });   // parallel_for (rethrows any exception to the caller)
        prof_par += std::chrono::duration<double>(std::chrono::steady_clock::now() - pp0).count();
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            const ElemOut& eo_r = elem_out[(size_t)e];
            const int* eqd = eo_r.eqd;
            const auto& Ke = eo_r.Ke; const auto& fe = eo_r.fe; const auto& Le = eo_r.Le;
            const auto& He = eo_r.He; const auto& Se = eo_r.Se;
            for (int a = 0; a < N; ++a)
                for (int ca = 0; ca < 2; ++ca) {
                    const int ra = eqd[2 * a + ca];
                    if (ra < 0) continue;
                    f_int(ra) += fe(2 * a + ca);
                }
        }
        // The structural elements are linear here, so their internal force is K_s.dv exactly --
        // no state, no return mapping. Adding it to BOTH the tangent and the internal force is
        // what makes residual(dv) = 0 hold at the solution rather than approximately.
        if (struct_k) f_int.noalias() += (*struct_k) * dv;
        // By their own laws: their internal force at the total displacement, from their committed
        // state -- an anchor or a geogrid at its capacity, a joint or a pile skin slipping.
        if (structure_laws) {
            structure_laws->add(*v_phase + dv, f_int, nullptr);
            struct_dv = dv;
        }
    };
    // THE MATRIX IS BUILT ONLY WHERE IT IS SOLVED WITH. The balance of an iterate needs the
    // internal force and the coupling blocks, and those blocks -- L, H and the storage -- do
    // not change inside a time step (the permeability is read at its start). Every trial of a
    // line search used to add the whole system to a sparse builder, sort it and throw it away:
    // 2489 systems for the 553 that were solved with, on the measurement above. The element
    // matrices of the LAST assembly are kept; the two functions below read them.
    const auto system_matrix = [&](double dt) {
        math::SparseMatrixBuilder Ab(NT);
        Ab.reserve((std::size_t)mesh.element_count * (3 * N) * (3 * N));
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            const ElemOut& eo_r = elem_out[(size_t)e];
            const int* eqd = eo_r.eqd;
            for (int a = 0; a < N; ++a)
                for (int ca = 0; ca < 2; ++ca) {
                    const int ra = eqd[2 * a + ca];
                    if (ra < 0) continue;
                    for (int b = 0; b < N; ++b) {
                        for (int cb = 0; cb < 2; ++cb) {
                            const int rb = eqd[2 * b + cb];
                            if (rb >= 0) Ab.add_entry(ra, rb, eo_r.Ke(2 * a + ca, 2 * b + cb));
                        }
                        const int pb = pore_eq[mesh.node_of(e, b)];
                        if (pb >= 0) {
                            Ab.add_entry(ra, ndisp + pb, eo_r.Le(2 * a + ca, b));
                            Ab.add_entry(ndisp + pb, ra, eo_r.Le(2 * a + ca, b));
                        }
                    }
                }
            for (int a = 0; a < N; ++a) {
                const int pa = pore_eq[mesh.node_of(e, a)];
                if (pa < 0) continue;
                for (int b = 0; b < N; ++b) {
                    const int pb = pore_eq[mesh.node_of(e, b)];
                    if (pb < 0) continue;
                    Ab.add_entry(ndisp + pa, ndisp + pb, -(dt * eo_r.He(a, b) + eo_r.Se(a, b)));
                }
            }
        }
        add_structural_block(Ab, struct_k);
        if (structure_laws) {   // their tangent where the internal force was last evaluated
            Eigen::VectorXd scratch = Eigen::VectorXd::Zero(ndisp);
            structure_laws->add(*v_phase + struct_dv, scratch, &Ab);
        }
        return Ab.build();
    };
    // The coupling blocks of a time step: L (force of the pore pressure), its transpose (volume
    // change of the displacement), H (flow) and dt H + S (flow and storage of the step).
    struct Coupling { math::CsrMatrix L, Lt, H, HS; };
    const auto coupling_blocks = [&](double dt) {
        const int np1 = std::max(1, npore);
        math::SparseMatrixBuilder Lb(ndisp, np1), Ltb(np1, ndisp), Hb(np1), HSb(np1);
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!active.empty() && !active[e]) continue;
            const ElemOut& eo_r = elem_out[(size_t)e];
            const int* eqd = eo_r.eqd;
            for (int a = 0; a < N; ++a)
                for (int ca = 0; ca < 2; ++ca) {
                    const int ra = eqd[2 * a + ca];
                    if (ra < 0) continue;
                    for (int b = 0; b < N; ++b) {
                        const int pb = pore_eq[mesh.node_of(e, b)];
                        if (pb < 0) continue;
                        Lb.add_entry(ra, pb, eo_r.Le(2 * a + ca, b));
                        Ltb.add_entry(pb, ra, eo_r.Le(2 * a + ca, b));
                    }
                }
            for (int a = 0; a < N; ++a) {
                const int pa = pore_eq[mesh.node_of(e, a)];
                if (pa < 0) continue;
                for (int b = 0; b < N; ++b) {
                    const int pb = pore_eq[mesh.node_of(e, b)];
                    if (pb < 0) continue;
                    Hb.add_entry(pa, pb, eo_r.He(a, b));
                    HSb.add_entry(pa, pb, dt * eo_r.He(a, b) + eo_r.Se(a, b));
                }
            }
        }
        return Coupling{Lb.build(), Ltb.build(), Hb.build(), HSb.build()};
    };

    Eigen::VectorXd v = Eigen::VectorXd::Zero(ndisp);
    v_phase = &v;
    struct_dv = Eigen::VectorXd::Zero(ndisp);
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
    const bool debug = std::getenv("KATAI_NL_DEBUG") != nullptr;
    // THE FORCE A STEP LEAVES OUT OF BALANCE IS CARRIED INTO THE NEXT ONE. A time step is solved
    // as an increment from the state it starts in, and it is accepted with up to the tolerated
    // error still out of balance. Dropped there, that remainder is lost for good -- and a phase
    // of some hundreds of steps loses it some hundreds of times. With the tolerance at 1e-6 the
    // sum was nothing; at the 1e-2 of a hardening soil it was measured as half the settlement of
    // a column consolidated to rest (337 mm against the 683 mm of the same column loaded
    // drained). Carried, each step's remainder is the next step's load and the error of the
    // phase is that of its last step.
    Eigen::VectorXd oob = Eigen::VectorXd::Zero(ndisp), oob_next = Eigen::VectorXd::Zero(ndisp);
    double prof_asm = 0, prof_build = 0, prof_factor = 0, prof_solve = 0; int prof_neval = 0, prof_nfac = 0;
    int step_depth = 0;          // how many times the step being solved has been cut in two
    double kept_share = 0.0;     // the balance a stalled step was kept at, as a share of its scale
    const auto one_step = [&](double h, const Eigen::VectorXd& df_step, Eigen::VectorXd& dv,
                              Eigen::VectorXd& dpv, bool& stagnant) -> bool {
        stagnant = false;
        small_frozen = false;   // each time step starts with the stiffness read from itself
        if (structure_laws) {   // the excess pore pressure a joint reads beside it: the step's start
            std::vector<double> pn(mesh.node_count, 0.0);
            for (int n = 0; n < mesh.node_count; ++n) if (pore_eq[n] >= 0) pn[n] = p(pore_eq[n]);
            structure_laws->set_soil_pore(pn);
        }
        Eigen::VectorXd Bbase;
        { Eigen::VectorXd dvz = Eigen::VectorXd::Zero(ndisp);
          assemble(dvz, Bbase, tmode, h); }   // trial == committed -> Bbase = f_int(committed)
        const Coupling cpl = npore > 0 ? coupling_blocks(h) : Coupling{};
        int eval_count = 0, elem_of_eval = -1;   // which evaluation the kept element matrices are of
        dv = Eigen::VectorXd::Zero(ndisp);
        dpv = Eigen::VectorXd::Zero(std::max(1, npore));
        struct Eval {
            int stamp = -1;   // the evaluation this is; the matrix is built from it on demand
            Eigen::VectorXd r;
            double ru = 0.0, rp = 0.0, ref_u = 0.0, ref_p = 0.0;
        };
        const auto evaluate = [&](const Eigen::VectorXd& dv_, const Eigen::VectorXd& dpv_) {
            Eval ev;
            Eigen::VectorXd f_int;
            const auto pt0 = std::chrono::steady_clock::now();
            assemble(dv_, f_int, tmode, h);
            prof_asm += std::chrono::duration<double>(std::chrono::steady_clock::now() - pt0).count();
            ++prof_neval;
            ev.stamp = ++eval_count;
            elem_of_eval = ev.stamp;
            // Newton RHS = -R: equilibrium R_u = (f_int - Bbase) + L dp - df; continuity from the
            // assembled bottom block (A x).tail = Lt dv - S* dp (S* = h H + S).
            ev.r = Eigen::VectorXd::Zero(NT);
            ev.r.head(ndisp) = df_step + oob - (f_int - Bbase);
            double pore_force = 0.0;
            if (npore > 0) {
                ev.r.head(ndisp).noalias() -= cpl.L * dpv_;
                // the bottom row of the system times the iterate: L^T dv - (h H + S) dp
                const Eigen::VectorXd Ltdv = cpl.Lt * dv_;
                const Eigen::VectorXd hHp = h * (cpl.H * p);
                ev.r.tail(npore) = hHp - (Ltdv - cpl.HS * dpv_);
                ev.rp = ev.r.tail(npore).norm();
                // The volume balance's own terms: the flow in the step and the volume change the
                // displacement makes (L^T dv). The assembled bottom row alone is no scale -- it IS
                // the balance, zero to round-off once a step has converged.
                ev.ref_p = std::max(hHp.norm(), Ltdv.norm());   // |L^T dv|
                pore_force = (cpl.L * (p + dpv_)).norm();
            }
            ev.ru = ev.r.head(ndisp).norm();
            ev.ref_u = (df_step + oob).norm() + Bbase.norm() + pore_force;
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
            // A balance that is not a number is worse than any that is. Written the other way
            // round -- "not worse unless it compares greater" -- a trial step that overflowed
            // the stress integration passed as an improvement, because nothing compares
            // greater than a NaN: measured on an embankment over soft clay on a 1.2 m mesh,
            // the second iterate of the first time step was NaN and every one after it.
            if (!std::isfinite(ma)) return true;
            return ma >= mb;
        };
        const auto rp_ok = [&](const Eval& ev) {
            return npore == 0 || ev.rp <= newton_tol * ev.ref_p + 1e-15;
        };
        Eval cur = evaluate(dv, dpv);
        Eigen::VectorXd best_dv = dv, best_dpv = dpv;
        bool best_rp_ok = rp_ok(cur);
        double best_ru = best_rp_ok ? cur.ru : std::numeric_limits<double>::infinity();
        // AN ITERATION THAT HAS STOPPED MOVING IS STOPPED. Where the tangent no longer gives a
        // direction of descent the line search ends at its shortest length and the balance
        // stays where it is; the rest of the budget -- forty iterations of fourteen assemblies,
        // then the same again on each half of the step down to a sixty-fourth -- bought
        // nothing (measured: 2 000 iterations on one time step, the balance at 1.035e-3 of its
        // scale throughout). Six such iterations in a row end it; what is then done with the
        // iterate is decided below, as for one that ran out of budget.
        int stalled = 0;
        std::vector<double> recent_ru;   // the balance at each iterate since the last change of method
        // ANDERSON ACCELERATION OF THE ELASTIC-STIFFNESS ITERATION. With the elastic operator
        // each step is K_e^-1 r, a fixed-point iteration that contracts by the ratio of the
        // ground's stiffness to the elastic one -- a percent or less per step where the ground
        // is yielding (measured on fill placed on soft clay: 8.8e-3 of the scale to 8.6e-3 in
        // seven iterations). The steps taken so far say how the step itself changes with the
        // iterate; the combination of the last few that makes the present step smallest in the
        // least-squares sense is taken out of it (Anderson 1965; Walker & Ni 2011, SIAM J.
        // Numer. Anal. 49). The memory is six steps and is dropped whenever the method changes
        // or a step has to be searched.
        constexpr int kAndersonDepth = 6;
        std::vector<Eigen::VectorXd> aa_dx, aa_dd;
        Eigen::VectorXd aa_x_prev, aa_d_prev;
        bool aa_has = false;
        tangent_floor = 0.0;
        elastic_matrix = false;
        int budget = max_newton;
        small_frozen = false;
        int small_refresh = 0;
        int it = 0;
        for (;;) {
        for (; it < budget && !converged(cur); ++it) {
            if (stalled >= 6 && has_small && !small_frozen) {
                // Stalled on a stiffness that follows the iterate: continued with the
                // small-strain stiffness read from the increment reached so far.
                small_frozen = true;
                dv_small = dv;
                stalled = 0;
                recent_ru.clear();
                budget = it + max_newton;
                if (debug)
                    std::fprintf(stderr, "    cons h=%.4g it %d: stalled; continuing with the "
                                         "small-strain stiffness read from this iterate\n", h, it);
                cur = evaluate(dv, dpv);
                if (converged(cur)) break;
            }
            if (stalled >= 6) {
                // Stalled on the tangent: once more with the elastic stiffness, which is given
                // the iterations a linearly converging method needs. Stalled on that too, the
                // iterate is judged as it stands.
                if (elastic_matrix) break;
                elastic_matrix = true;
                stalled = 0;
                recent_ru.clear();
                aa_dx.clear(); aa_dd.clear(); aa_has = false;
                budget = it + 4 * max_newton;
                if (debug)
                    std::fprintf(stderr, "    cons h=%.4g it %d: stalled; continuing with the "
                                         "elastic stiffness\n", h, it);
                cur = evaluate(dv, dpv);
            }
            if (!solve_factory) return false;
            // A refused linear solve ends THIS time step, exactly as a non-converged Newton
            // does -- the tangent is singular at this iterate, which on this path means the soil
            // has reached its capacity under the load being consolidated. Only SingularSystem is
            // caught: a malformed request or a broken backend raises SolveError and still
            // propagates, because turning one of those into "did not converge" would publish a
            // modelling answer for a bug.
            Eigen::VectorXd d;
            bool refused = false;
            try {
                // The element matrices kept are those of the last evaluation; if that is not
                // this iterate (it always has been, the search leaves the trial states at the
                // iterate it accepts) the iterate is evaluated again.
                if (cur.stamp != elem_of_eval) cur = evaluate(dv, dpv);
                const auto pb0 = std::chrono::steady_clock::now();
                const math::CsrMatrix A_now = system_matrix(h);
                const auto pt0 = std::chrono::steady_clock::now();
                prof_build += std::chrono::duration<double>(pt0 - pb0).count();
                const auto solve = solve_factory(A_now);
                const auto pt1 = std::chrono::steady_clock::now();
                d = solve(cur.r);
                prof_factor += std::chrono::duration<double>(pt1 - pt0).count();
                prof_solve += std::chrono::duration<double>(std::chrono::steady_clock::now() - pt1).count();
                ++prof_nfac;
            } catch (const math::SingularSystem&) { refused = true; }
            if (refused || !d.allFinite()) {
                // Once, with the floor on the tangent of the plastic points; then it is the
                // step that cannot be taken.
                if (tangent_floor > 0.0) break;
                tangent_floor = kVertexFloor;
                if (debug)
                    std::fprintf(stderr, "    cons h=%.4g it %d: the solve was refused; retried "
                                         "with the tangent floor\n", h, it);
                cur = evaluate(dv, dpv);
                continue;
            }
            if (elastic_matrix) {
                Eigen::VectorXd x(NT);
                x.head(ndisp) = dv;
                if (npore > 0) x.tail(npore) = dpv;
                if (aa_has) {
                    aa_dx.push_back(x - aa_x_prev);
                    aa_dd.push_back(d - aa_d_prev);
                    if ((int)aa_dx.size() > kAndersonDepth) {
                        aa_dx.erase(aa_dx.begin());
                        aa_dd.erase(aa_dd.begin());
                    }
                }
                aa_x_prev = x;
                aa_d_prev = d;
                aa_has = true;
                if (!aa_dd.empty()) {
                    // Least squares on the displacement part: the pore part has its own units.
                    const int k = (int)aa_dd.size();
                    Eigen::MatrixXd F(ndisp, k);
                    for (int j = 0; j < k; ++j) F.col(j) = aa_dd[(size_t)j].head(ndisp);
                    const Eigen::VectorXd gam =
                        F.colPivHouseholderQr().solve(Eigen::VectorXd(d.head(ndisp)));
                    if (gam.allFinite()) {
                        Eigen::VectorXd acc = d;
                        for (int j = 0; j < k; ++j)
                            acc -= gam(j) * (aa_dx[(size_t)j] + aa_dd[(size_t)j]);
                        if (acc.allFinite()) d = acc;
                    }
                }
            }
            // Developer diagnostic (KATAI_CONS_TANGENT_CHECK): is the matrix the derivative of
            // the balance along the direction just solved for? Along a Newton direction the
            // balance at a fraction a of the step is (1 - a) r if it is; what is printed is the
            // departure from that as a share of a |r|, at four sizes of a, and the nodes that
            // carry it. A share that does not fall with a is a balance that is not smooth there.
            if (debug && std::getenv("KATAI_CONS_TANGENT_CHECK") != nullptr) {
                const Eigen::VectorXd r0 = cur.r.head(ndisp);
                for (double a : {1e-2, 1e-4, 1e-6, 1e-8}) {
                    Eigen::VectorXd dp_a = dpv;
                    if (npore > 0) dp_a += a * d.tail(npore);
                    const Eval ea = evaluate(dv + a * d.head(ndisp), dp_a);
                    const Eigen::VectorXd miss = ea.r.head(ndisp) - (1.0 - a) * r0;
                    std::fprintf(stderr, "      tangent check a=%.0e: departure %.3e of a|r|\n", a,
                                 miss.norm() / (a * std::max(r0.norm(), 1e-300)));
                    if (a == 1e-6) {
                        std::vector<std::pair<double, int>> big;
                        for (int q = 0; q < ndisp; ++q) big.push_back({std::fabs(miss(q)), q});
                        const size_t nb = std::min<size_t>(6, big.size());
                        std::partial_sort(big.begin(), big.begin() + nb, big.end(),
                                          [](const auto& x, const auto& y) { return x.first > y.first; });
                        for (size_t b = 0; b < nb; ++b)
                            for (int nd = 0; nd < mesh.node_count; ++nd)
                                for (int c = 0; c < 2; ++c)
                                    if (dofs.equation(dofs.global_dof(nd, c)) == big[b].second)
                                        std::fprintf(stderr, "        %.3e at node %d (%.3f, %.3f) %c\n",
                                                     miss(big[b].second), nd, mesh.x[nd], mesh.y[nd], c ? 'y' : 'x');
                    }
                }
                (void)evaluate(dv, dpv);   // the trial states back at the iterate
            }
            double alpha = 1.0;
            Eigen::VectorXd dv_n, dpv_n;
            Eval next;
            // The step lengths tried are remembered, and when none of them improves the balance
            // the LEAST BAD is taken, not the last. The last one is the shortest, and on a
            // hardening soil -- whose tangent is not consistent with its integration -- the
            // shortest step is not the mildest: measured on an embankment over soft clay, a
            // search that ran out at 1/64 took the balance from 9e-2 of its scale to 3e+3, and
            // the rest of the iteration budget went on coming back from there.
            double best_alpha = 1.0, best_measure = std::numeric_limits<double>::infinity();
            bool improved = false;
            // The search goes down to 2^-13 of the step. A soil whose stiffness is exponential in
            // the volumetric strain (the compression-index models) turns a trial strain a few
            // times its compression index into a stress many orders too large, so a step from a
            // soft tangent can be worse than the iterate at 1/64 and still an improvement at
            // 1/1000; stopping at 1/64 took such steps as they were (measured: the balance of a
            // dissipation step went from 0.16 of its scale to 1e+3, then to 1e+17).
            constexpr int kLineSearchHalvings = 14;
            for (int ls = 0; ls < kLineSearchHalvings; ++ls) {
                dv_n = dv + alpha * d.head(ndisp);
                dpv_n = dpv;
                if (npore > 0) dpv_n += alpha * d.tail(npore);
                next = evaluate(dv_n, dpv_n);
                if (!worse(next, cur) || converged(next)) { improved = true; break; }
                // THE ELASTIC-STIFFNESS ITERATION TAKES ITS FULL STEP, as it does in a Plastic
                // phase: it is a fixed-point iteration, its balance does not fall monotonically
                // on a response that is not smooth, and a search that insists it must takes
                // the least bad of fourteen shrinking steps instead and stays where it is
                // (measured on fill placed on soft clay: the balance held at 2.26e-2 of its
                // scale with steps of a five-thousandth). Only a step that throws the balance
                // far out -- tenfold -- is treated as a failure and searched.
                if (elastic_matrix && ls == 0) {
                    const double mc = cur.ru / std::max(cur.ref_u, 1e-300) +
                                      (npore > 0 ? cur.rp / std::max(cur.ref_p, 1e-300) : 0.0);
                    const double mn = next.ru / std::max(next.ref_u, 1e-300) +
                                      (npore > 0 ? next.rp / std::max(next.ref_p, 1e-300) : 0.0);
                    if (std::isfinite(mn) && mn < 10.0 * mc) { improved = true; break; }
                }
                const double m = next.ru / std::max(next.ref_u, 1e-300) +
                                 (npore > 0 ? next.rp / std::max(next.ref_p, 1e-300) : 0.0);
                if (m < best_measure) { best_measure = m; best_alpha = alpha; }
                alpha *= 0.5;
            }
            // WITH THE ELASTIC STIFFNESS A FULL STEP IS A SHORT ONE. The elastic operator is
            // stiffer than the ground it stands in for, so its step falls short of the balance by
            // the ratio of the two, every time: measured, the balance fell by 0.3% an iteration.
            // The step is therefore lengthened -- doubled while the balance keeps improving, up
            // to sixteen times -- which is the line search run the other way.
            if (improved && elastic_matrix && alpha == 1.0 && !converged(next)) {
                Eval longer;
                Eigen::VectorXd dv_l, dpv_l;
                double a = 2.0;
                for (; a <= 16.0; a *= 2.0) {
                    dv_l = dv + a * d.head(ndisp);
                    dpv_l = dpv;
                    if (npore > 0) dpv_l += a * d.tail(npore);
                    longer = evaluate(dv_l, dpv_l);
                    if (worse(longer, next)) break;
                    next = longer;
                    dv_n = dv_l;
                    dpv_n = dpv_l;
                    alpha = a;
                    if (converged(next)) break;
                }
                // The last evaluation set the trial states; if it was a rejected length, put
                // them back at the one kept.
                if (a <= 16.0 && alpha != a && !converged(next)) next = evaluate(dv_n, dpv_n);
            }
            // A STEP THAT ONLY MAKES THE BALANCE WORSE IS NOT TAKEN. The least bad of the lengths
            // tried used to be taken whatever it did. With stress-free points in the mesh --
            // freshly placed fill held on the tip of its cap answers a pull with no stiffness --
            // the matrix can be all but singular without being refused, the direction it gives
            // is then of no use at any length, and the least bad of fourteen lengths was a step
            // of a ten-thousandth that took the balance from 0.04 of its scale to 3.9e7; the
            // iterate never came back (measured on half a metre of fill placed on soft clay,
            // 15-noded). Where the best length still doubles the balance the iterate stays
            // where it is and the iteration is treated as stalled: the next method takes over.
            if (!improved) {
                const double m_cur = cur.ru / std::max(cur.ref_u, 1e-300) +
                                     (npore > 0 ? cur.rp / std::max(cur.ref_p, 1e-300) : 0.0);
                if (!(best_measure <= 2.0 * m_cur)) {
                    cur = evaluate(dv, dpv);   // the trial states back at the iterate
                    if (debug)
                        std::fprintf(stderr, "    cons h=%.4g it %d: no step length within twice "
                                             "the balance (best %.3e against %.3e); not taken\n",
                                     h, it, best_measure, m_cur);
                    if (elastic_matrix) { aa_dx.clear(); aa_dd.clear(); aa_has = false; }
                    stalled = 6;
                    if (elastic_matrix) break;
                    continue;
                }
            }
            if (!improved && best_alpha != alpha * 2.0) {   // the last one tried was not the best
                alpha = best_alpha;
                dv_n = dv + alpha * d.head(ndisp);
                dpv_n = dpv;
                if (npore > 0) dpv_n += alpha * d.tail(npore);
                next = evaluate(dv_n, dpv_n);
            } else if (!improved) {
                alpha *= 2.0;   // the length actually taken, for the trace
            }
            // No step length gave a finite balance: this time step cannot be taken from here,
            // and is cut in two like any other that does not close.
            if (!std::isfinite(next.ru) || !std::isfinite(next.rp)) {
                (void)evaluate(dv, dpv);   // leave the trial states at the last finite iterate
                return false;
            }
            if (elastic_matrix && alpha < 1.0) { aa_dx.clear(); aa_dd.clear(); aa_has = false; }
            stalled = (next.ru > 0.98 * cur.ru && alpha < 1e-2) ? stalled + 1 : 0;
            // ... AND ONE THAT MOVES WITHOUT GETTING ANYWHERE IS STOPPED TOO. The rule above
            // asks for six shortest steps in a row. An iteration can go on taking steps of a
            // tenth and a hundredth of the Newton direction for its whole budget with the
            // balance where it was -- measured on the first time step of an embankment built on
            // 15-noded elements: 9.07e-3 of its scale at iteration 20 and 9.01e-3 at 59, seven
            // seconds an iteration, the step then halved and the same again. If the best of the
            // last eight iterations is not a tenth below the best of the eight before them, the
            // iteration has stalled, whatever its step lengths were.
            recent_ru.push_back(next.ru);
            if (recent_ru.size() >= 16) {
                const auto b = recent_ru.end();
                const double last = *std::min_element(b - 8, b);
                const double before = *std::min_element(b - 16, b - 8);
                if (last > 0.9 * before) stalled = 6;
            }
            dv = std::move(dv_n);
            dpv = std::move(dpv_n);
            cur = std::move(next);
            if (debug)
                std::fprintf(stderr, "    cons h=%.4g it %d  ru/ref=%.3e  rp/ref=%.3e  alpha=%.4g\n", h,
                             it, cur.ru / std::max(cur.ref_u, 1e-300),
                             npore > 0 ? cur.rp / std::max(cur.ref_p, 1e-300) : 0.0, alpha);
            if (rp_ok(cur) && cur.ru < best_ru) {
                best_ru = cur.ru;
                best_dv = dv;
                best_dpv = dpv;
                best_rp_ok = true;
            }
        }
        // The answer of a fixed-stiffness pass becomes the increment the stiffness is read from,
        // and the step is closed again from it.
        if (!(small_frozen && small_refresh < 2 && converged(cur))) break;
        if ((dv - dv_small).norm() <= 1.0e-6 * std::max(dv.norm(), 1.0e-300)) break;
        dv_small = dv;
        ++small_refresh;
        stalled = 0;
        budget = it + max_newton;
        cur = evaluate(dv, dpv);
        }
        elastic_matrix = false;
        tangent_floor = 0.0;
        if (converged(cur)) { oob_next = cur.r.head(ndisp); return true; }
        // A stalled iterate is kept when it is within kStagnationAccept of the force scale, or --
        // for the soils solved to 1e-3 -- within three times their tolerance. Measured on a Soft
        // Soil column consolidated to rest, 3e-3 is 0.4% from the answer at 1e-6 and 1e-2 is 7%.
        // Counted and reported either way (K2D-A019), and what is left out of balance is
        // carried into the next step.
        // A STEP THAT HAS BEEN CUT AND STALLS AGAIN IS KEPT WITHIN ONE PERCENT. The tolerated
        // error geotechnical finite element analyses are run with by default is one percent of
        // the forces in play, and there are states this solver cannot take further than that
        // either: ground at its strength with no dilatancy has no equilibrium to a tighter
        // tolerance on some meshes (the non-associated flow rule; Tschuchnigg, Schweiger and
        // Sloan 2015), and freshly placed fill beside a spreading base is such ground. Measured
        // on fill placed on soft clay on 15-noded elements: the balance of the first time
        // steps comes down to 1.07e-2 of its scale on every method and stays there, at every
        // size of step down to a sixty-fourth, where the same model on 6-noded elements
        // closes. The first attempt at a step is held to the tolerance of the phase; once it
        // has been cut in two and stalls again, the best iterate is kept if it is within
        // kLastResortAccept, as an increment of a Plastic phase that cannot be cut further is.
        // It is counted and reported with what it leaves out of balance (K2D-A019), and that
        // remainder is carried into the next step. Measured on the 6-noded model, which does
        // close: the settlement moves by 0.7% at ten time steps and not at all at forty.
        const double stall_accept = step_depth >= 1
                                        ? std::max({kStagnationAccept, 3.0 * newton_tol, kLastResortAccept})
                                        : std::max(kStagnationAccept, 3.0 * newton_tol);
        if (!best_rp_ok || !(best_ru <= stall_accept * cur.ref_u)) return false;
        kept_share = best_ru / std::max(cur.ref_u, 1e-300);
        dv = best_dv;
        dpv = best_dpv;
        oob_next = evaluate(dv, dpv).r.head(ndisp);   // the trial states of the iterate being kept
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
        step_depth = depth;
        if (one_step(h, df_step, dv, dpv, stagnant)) {
            v += dv;
            if (npore > 0) p += dpv;
            committed = trial;
            if (structure_laws) structure_laws->commit();   // their trial state is at the step kept
            oob = oob_next;
            if (stagnant) {
                ++R.stagnation_accepted;
                R.worst_accepted = std::max(R.worst_accepted, kept_share);
            }
            record(t0 + h);
            return true;
        }
        trial = committed;
        if (debug)
            std::fprintf(stderr, "  cons t0=%.5g h=%.4g depth %d: step not closed%s\n", t0, h, depth,
                         depth >= 6 ? " -- GIVEN UP" : ", cut in two");
        if (depth >= 6) return false;
        depth_used = std::max(depth_used, depth + 1);
        const Eigen::VectorXd half = 0.5 * df_step;
        return advance(t0, 0.5 * h, half, depth + 1) && advance(t0 + 0.5 * h, 0.5 * h, half, depth + 1);
    };
    const Eigen::VectorXd no_load = Eigen::VectorXd::Zero(ndisp);
    Eigen::VectorXd part;   // the load of one step when it is placed over the whole duration
    if (load_increment && ramp_load) part = *load_increment / (double)nsteps;
    Eigen::VectorXd first;   // the first step also carries what the phase inherits, at once
    if (instant_load) {
        first = *instant_load;
        if (load_increment) first += ramp_load ? part : *load_increment;
    }
    for (int step = 0; step < nsteps && R.converged; ++step) {
        const Eigen::VectorXd& df = (step == 0 && instant_load) ? first
                                    : (load_increment && ramp_load) ? part
                                    : (step == 0 && load_increment) ? *load_increment : no_load;
        if (!advance(step * dt, dt, df, 0)) R.converged = false;
    }
    R.steps_cut = depth_used;
    R.committed = committed;
    if (std::getenv("KATAI_NL_PROFILE"))
        std::fprintf(stderr, "CONSPROF threads %d par %.1fs evals %d asm %.1fs build %.1fs | factors %d factor %.1fs "
                             "solve %.1fs%c",
                     math::parallel_threads(), prof_par, prof_neval, prof_asm, prof_build, prof_nfac, prof_factor, prof_solve, 10);
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
                                               const std::vector<int>* pore_tie, bool ramp_load,
    const Eigen::VectorXd* instant_load) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_impl<Tri15Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                        kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                        active, load_increment, solve_factory, struct_k,
                                                        pore_tie, ramp_load, instant_load);
    return detail::consolidation_impl<Tri6Element>(mesh, dofs, materials, profile, perm, gamma_w,
                                                   kw_over_n, drained_node, initial_pore, dt, nsteps,
                                                   active, load_increment, solve_factory, struct_k,
                                                   pore_tie, ramp_load, instant_load);
}

ConsolidationPlasticResult solve_consolidation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, double gamma_w, const PoreFluidStiffness& kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton, double newton_tol,
    const std::vector<MaterialProfile>& profile, const math::CsrMatrix* struct_k,
    const std::vector<int>* pore_tie, bool ramp_load,
    const Eigen::VectorXd* instant_load, bool axisymmetric, StructureLaws* structure_laws) {
    if (axisymmetric) {
        if (mesh.nodes_per_element == Tri15Element::kNodeCount)
            return detail::consolidation_plastic_impl<Tri15Element, detail::AxisymKin>(
                mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
                initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
                struct_k, pore_tie, ramp_load, instant_load, structure_laws);
        return detail::consolidation_plastic_impl<Tri6Element, detail::AxisymKin>(
            mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
            initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
            struct_k, pore_tie, ramp_load, instant_load, structure_laws);
    }
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::consolidation_plastic_impl<Tri15Element>(
            mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
            initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
            struct_k, pore_tie, ramp_load, instant_load, structure_laws);
    return detail::consolidation_plastic_impl<Tri6Element>(
        mesh, dofs, materials, profile, perm, gamma_w, kw_over_n, drained_node, initial_state,
        initial_pore, dt, nsteps, active, load_increment, solve_factory, max_newton, newton_tol,
        struct_k, pore_tie, ramp_load, instant_load, structure_laws);
}

}  // namespace katai::core
