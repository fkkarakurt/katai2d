#pragma once
// SHARED internal-force + consistent-tangent assembly — soil (tri6/tri15 × plane-strain/
// axisym) constitutive return mapping + ALL structural elements embedded in the soil (plate
// 3/5, anchor, geogrid, interface 3/5, embedded beam). Produces f_int and (on request) K_T as
// a PURE function of the total displacement (u_free + du_free) given the committed Gauss/
// structural state.
//
// WHY SHARED: this machinery is used by two solvers under WORD-FOR-WORD the same rules:
//   (1) static Newton  — analysis/nonlinear_solver.cpp  (solve_nonlinear)
//   (2) nonlinear dynamic Newmark+Newton — analysis/dynamics_nonlinear.cpp (solve_newmark_nonlinear)
// Two copies = silent-drift risk (the same decision was made for make_diagram/
// make_iface_diagram: the static tail and the dynamic envelope call the SAME lambda). This
// header is that single-source assembly. Moving it here does NOT change the behaviour of
// solve_nonlinear (bit-for-bit; all existing V&V is the oracle).
//
// A heavy template (pulls in all element headers + Eigen) → included ONLY by the two .cpp
// translation units above (NOT header-only), so 148-TU shared headers like build_problem are
// unaffected by it.

#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/nonlinear_solver.hpp>   // Structures + element structs + NewtonResult::Timings
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/fem/elements/geogrid.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/math/thread_pool.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {
namespace detail {

// Node coordinates of element e + global indices of its local DOFs.
template <class E>
inline void gather_element(const mesh::Mesh& mesh, const DofMap& dofs, int e,
                           typename E::NodeCoords& coords,
                           std::array<int, E::kDofCount>& element_dofs) {
    for (int k = 0; k < E::kNodeCount; ++k) {
        const int n = mesh.node_of(e, k);
        coords(k, 0) = mesh.x[n];
        coords(k, 1) = mesh.y[n];
        element_dofs[2 * k] = dofs.global_dof(n, 0);
        element_dofs[2 * k + 1] = dofs.global_dof(n, 1);
    }
}

// Kinematics policies — gather the kinematics-dependent part of the assembly (B matrix,
// strain/stress size, integration weight, material integration) into a single type. The outer
// Newton/adaptive-step loop is kinematics-independent.
struct PlaneStrainKin {
    static constexpr int kStrain = 3;
    using Strain = Eigen::Matrix<double, 3, 1>;
    using Tangent = Eigen::Matrix3d;
    template <class E>
    struct Grad { Eigen::Matrix<double, 3, E::kDofCount> B; double weight; };
    template <class E>
    static Grad<E> gradients(const typename E::NodeCoords& c, double xi, double eta) {
        const auto g = E::strain_displacement(c, xi, eta);
        return {g.B, g.det_jacobian};
    }
    static void integrate(const MaterialModel& m, const GaussState& comm,
                          const Strain& de, GaussState& tr, Tangent& t, TangentMode mode,
                          double creep_dt) {
        integrate_point(m, comm, de, tr, t, mode, creep_dt);
    }
    static Eigen::Matrix<double, 3, 1> stress(const GaussState& s) { return s.stress; }
    // Pore-pressure direction m: normal components only (shear unaffected).
    static Strain pore_vector() { return Strain(1.0, 1.0, 0.0); }
};

struct AxisymKin {
    static constexpr int kStrain = 4;
    using Strain = Eigen::Matrix<double, 4, 1>;
    using Tangent = Eigen::Matrix4d;
    template <class E>
    struct Grad { Eigen::Matrix<double, 4, E::kDofCount> B; double weight; };
    template <class E>
    static Grad<E> gradients(const typename E::NodeCoords& c, double xi, double eta) {
        const auto g = axisym::strain_displacement<E>(c, xi, eta);
        return {g.B, g.det_jacobian * g.radius};  // r-weighted
    }
    static void integrate(const MaterialModel& m, const GaussState& comm,
                          const Strain& de, GaussState& tr, Tangent& t, TangentMode mode,
                          double creep_dt) {
        integrate_point_axisym(m, comm, de, tr, t, mode, creep_dt);
    }
    static Eigen::Matrix<double, 4, 1> stress(const GaussState& s) {
        Eigen::Matrix<double, 4, 1> v;
        v << s.stress, s.stress_zz;
        return v;
    }
    // Pore-pressure direction m = [r, z, rz, theta]: the hoop is a real normal stress.
    static Strain pore_vector() {
        Strain m;
        m << 1.0, 1.0, 0.0, 1.0;
        return m;
    }
};

// Internal-force/tangent assembler templated over element (tri6/tri15) × kinematics
// (plane-strain/axisym). Element topology + fe/ke buffers LIVE for the whole solve (gathered
// once in the constructor); committed/trial state pointers come in via State on every call.
// The `assemble` call is the verbatim body of the ORIGINAL solve_nonlinear_impl assemble
// lambda (only captures → members/arguments).
template <class E, class Kin>
class InternalForceAssembler {
public:
    using ElementVector = Eigen::Matrix<double, E::kDofCount, 1>;
    static constexpr int n_gp = E::kGaussCount;

    // Path-dependent state pointers (committed = start of step; trial = WRITTEN here). Both
    // point into the caller's vectors; the caller owns them.
    struct State {
        const std::vector<GaussState>* committed = nullptr;
        std::vector<GaussState>* trial = nullptr;
        const std::vector<double>* anchor_c = nullptr;  std::vector<double>* anchor_t = nullptr;
        const std::vector<double>* geogrid_c = nullptr; std::vector<double>* geogrid_t = nullptr;
        const std::vector<double>* iface_c = nullptr;   std::vector<double>* iface_t = nullptr;
        const std::vector<double>* iface5_c = nullptr;  std::vector<double>* iface5_t = nullptr;
        const std::vector<double>* eskin_c = nullptr;   std::vector<double>* eskin_t = nullptr;
        const std::vector<double>* efoot_c = nullptr;   std::vector<double>* efoot_t = nullptr;
        // Plate M-N hinge state [ε_p,κ_p]×Gauss (plate::kPlasticStateSize(5) per element).
        const std::vector<double>* plate_c = nullptr;   std::vector<double>* plate_t = nullptr;
        const std::vector<double>* plate5_c = nullptr;  std::vector<double>* plate5_t = nullptr;
    };

    // Nonzero-Dirichlet (prescribed displacement ū) ramp — used ONLY by the static path.
    // presc==null → no fixed-DOF displacement at all (the dynamic seismic path, in the
    // relative frame, does not use this). factor = (cur_target − cur_lambda): this step's
    // increment fraction.
    // KNOWN HIDDEN LIMIT (declared): the ramp enters only the SOIL element loop — the
    // structural element loops (plate/anchor/geogrid/interface/embedded) never read fixed
    // DOFs (geq<0 → contribution 0), so a structural element on a boundary driven by a
    // prescribed displacement does NOT see that motion. Prescribed ū is today a kernel/test
    // path (not available as a deformation BC in the GUI); if added to the GUI, the TOTAL
    // λ·ū must be carried into the structural elements (total-displacement formulation; the
    // increment is not enough). Measured: in the first design of test_plate_plastic the
    // plate did not see the driven middle node (M ~ 0).
    struct Ramp {
        const Eigen::VectorXd* presc = nullptr;   // total_dofs
        double factor = 0.0;
    };

    InternalForceAssembler(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<MaterialModel>& materials,
                           const std::vector<char>& active_element,
                           const Structures& structures,
                           const std::vector<MaterialProfile>& profile)
        : mesh_(mesh), dofs_(dofs), materials_(materials), active_(active_element),
          structures_(structures), profile_(profile), neq_(dofs.equation_count()) {
        // Element topology is CONSTANT for the whole solve → coordinates + DOF mapping are
        // gathered once (drops the gather cost per iteration; read-only in the parallel phase).
        all_coords_.resize(mesh.element_count);
        all_edofs_.resize(mesh.element_count);
        for (int e = 0; e < mesh.element_count; ++e)
            gather_element<E>(mesh, dofs, e, all_coords_[e], all_edofs_[e]);
        fe_buf_.resize(mesh.element_count);
    }

    int equation_count() const { return neq_; }

    // TIME share of this increment [days] — read only by time-dependent constitutive models
    // (SoftSoilCreep); 0 = no creep (all old paths bit-for-bit). The static solver writes
    // time_interval·Δλ per increment; consolidation/dynamics write their own step duration
    // (converted to days).
    double dt_day = 0.0;

    // Internal-force (and, if build_tangent, consistent-tangent) assembly. Since Δε = B·du_e,
    // f_int is a pure function of du_free (of u_free+du_free for structural elements) given
    // the committed state (total-increment formulation) — this makes it cheap for the line
    // search to recompute the residual without building K_T. state.trial[] is written on
    // every call. tmode: (hs_consistent_mode ? kConsistent : kContinuum) while build_tangent,
    // otherwise kNone is used.
    Eigen::VectorXd assemble(const Eigen::VectorXd& u_free, const Eigen::VectorXd& du_free,
                             bool build_tangent, TangentMode tmode, const State& st,
                             const Ramp& ramp, math::SparseMatrixBuilder* builder,
                             NewtonResult::Timings* timings = nullptr) {
        using Clock = std::chrono::steady_clock;
        const auto t_asm = Clock::now();
        auto elapsed = [](Clock::time_point t0) {
            return std::chrono::duration<double>(Clock::now() - t0).count();
        };
        const mesh::Mesh& mesh = mesh_;
        const DofMap& dofs = dofs_;
        const auto& materials = materials_;
        const auto& active_element = active_;
        const auto& structures = structures_;
        const auto& profile = profile_;
        const int neq = neq_;
        const auto gauss = E::gauss_points();
        const std::vector<GaussState>& committed = *st.committed;
        std::vector<GaussState>& trial = *st.trial;
        const bool has_presc = ramp.presc != nullptr;

        Eigen::VectorXd f_int = Eigen::VectorXd::Zero(neq);

        // --- Soil elements: TWO-PHASE assembly. (1) Per-element computation in PARALLEL —
        // each element writes only to its OWN fe/ke buffer and its OWN Gauss trial states (no
        // data race; no shared accumulator → floating-point order unchanged). (2) Scatter
        // SEQUENTIAL — f_int accumulation and COO entry order are IDENTICAL to the serial
        // path = the result is deterministic, independent of the thread count.
        if (build_tangent && ke_buf_.size() != static_cast<size_t>(mesh.element_count))
            ke_buf_.resize(mesh.element_count);
        math::parallel_for(mesh.element_count, [&](int e_begin, int e_end) {
        for (int e = e_begin; e < e_end; ++e) {
            if (!active_element.empty() && !active_element[e]) continue;  // passive (excavated)
            const typename E::NodeCoords& coords = all_coords_[e];
            const std::array<int, E::kDofCount>& edofs = all_edofs_[e];
            const int e_mat = mesh.element_material[e];
            const MaterialModel& mat = materials[e_mat];
            const MaterialProfile prof =
                e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};

            ElementVector du_e = ElementVector::Zero();  // fixed DOF → 0
            for (int a = 0; a < E::kDofCount; ++a) {
                const int eq = dofs.equation(edofs[a]);
                if (eq >= 0)
                    du_e(a) = du_free(eq);
                else if (has_presc)  // prescribed displacement: this step's increment of the fixed DOF
                    du_e(a) = ramp.factor * (*ramp.presc)(edofs[a]);
            }

            ElementVector fe = ElementVector::Zero();
            typename E::ElementMatrix ke = E::ElementMatrix::Zero();
            for (int g = 0; g < n_gp; ++g) {
                const auto grad =
                    Kin::template gradients<E>(coords, gauss[g].xi, gauss[g].eta);
                const typename Kin::Strain dstrain = grad.B * du_e;
                const int gi = e * n_gp + g;
                typename Kin::Tangent dt;
                const TangentMode tm = !build_tangent ? TangentMode::kNone : tmode;
                // Depth-varying stiffness / cohesion are evaluated HERE, at the stress point (PLAXIS
                // does the same): the whole point of a gradient is that it varies WITHIN an element,
                // so an element-average would quietly flatten it on a coarse mesh. uniform() keeps the
                // per-element material by reference -> the constant-E path stays bit-for-bit.
                const MaterialModel* mp = &mat;
                MaterialModel mg;
                if (!prof.uniform()) {
                    const typename E::ShapeValues sh = E::shape_functions(gauss[g].xi, gauss[g].eta);
                    double y = 0.0;
                    for (int i = 0; i < E::kNodeCount; ++i) y += sh(i) * coords(i, 1);
                    mg = mat;
                    mg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, y);
                    mg.cohesion = profile_at(mat.cohesion, prof.c_inc, prof.y_ref, y);
                    mp = &mg;
                }
                const MaterialModel& matg = *mp;
                Kin::integrate(matg, committed[gi], dstrain, trial[gi], dt, tm, dt_day);
                typename Kin::Strain sigma = Kin::stress(trial[gi]);
                // Undrained (A): the constitutive model returns EFFECTIVE stress and its effective
                // tangent; add the pore fluid's volumetric (bulk) contribution. Total = sigma' +
                // (Kw/n) eps_v m and tangent += (Kw/n) m m^T = D_u. The excess pore pressure
                // u = -(Kw/n) eps_v is carried in the Gauss state (eps_vol) for post-processing /
                // staging. (See effective-stress-formulation.md.)
                // matg, not mat: kw_over_n is derived from E' (K' = E'/(3(1-2nu'))), so under a depth
                // gradient the pore-fluid stiffness varies with depth too. Using `mat` here would
                // freeze it at the reference value while the skeleton stiffened -- a silent mismatch.
                // The accumulated volumetric strain is state every material may need, not only
                // an undrained one: the dilatancy cut-off reads it to know the void ratio
                // (material_model.hpp, void_ratio_of). It used to be tracked only where it was
                // consumed, which is why a drained soil could dilate for ever -- nothing was
                // counting. Accumulating it always changes no existing result, because nothing
                // else reads it for a drained material; the suite is what says so.
                const typename Kin::Strain mvec = Kin::pore_vector();
                trial[gi].eps_vol = committed[gi].eps_vol + mvec.dot(dstrain);
                if (matg.undrained) {
                    // Undrained (A): the constitutive model returned EFFECTIVE stress; add the
                    // pore fluid's volumetric contribution. u = -(Kw/n) eps_v.
                    const double kwn = matg.kw_over_n(matg.undrained_poisson);
                    sigma += kwn * trial[gi].eps_vol * mvec;
                    if (build_tangent) dt += kwn * mvec * mvec.transpose();
                }
                const double w = gauss[g].weight * grad.weight;
                fe.noalias() += w * grad.B.transpose() * sigma;
                if (build_tangent)
                    ke.noalias() += w * grad.B.transpose() * dt * grad.B;
            }
            fe_buf_[e] = fe;
            if (build_tangent) ke_buf_[e] = ke;
        }
        });  // parallel_for (rethrows any exception to the caller)

        for (int e = 0; e < mesh.element_count; ++e) {  // sequential scatter (deterministic)
            if (!active_element.empty() && !active_element[e]) continue;
            const std::array<int, E::kDofCount>& edofs = all_edofs_[e];
            const ElementVector& fe = fe_buf_[e];
            for (int a = 0; a < E::kDofCount; ++a) {
                const int eq_a = dofs.equation(edofs[a]);
                if (eq_a < 0) continue;
                f_int(eq_a) += fe(a);
                if (!build_tangent) continue;
                const typename E::ElementMatrix& ke = ke_buf_[e];
                for (int b = 0; b < E::kDofCount; ++b) {
                    const int eq_b = dofs.equation(edofs[b]);
                    if (eq_b < 0) continue;
                    builder->add_entry(eq_a, eq_b, ke(a, b));
                }
            }
        }

        // --- Plate (structural wall/beam) contribution. 9-DOF Timoshenko beam; translational
        // DOFs (u_x,u_y) are SHARED with the soil nodes (assembly summation → soil-structure
        // interaction), the rotational DOF (φ) is plate-specific. props.plastic() off ⇒
        // ELASTIC f_int_plate = K·u_total (u_total = committed u_free + increment du_free;
        // OLD path bit-for-bit). On ⇒ M-N hinge: return mapping at the bending/axial Gauss
        // points, state committed/trial (plate.hpp §10). The soil path is UNCHANGED.
        for (size_t pli = 0; pli < structures.plates.size(); ++pli) {
            const auto& pe = structures.plates[pli];
            plate::NodeCoords Xe;
            std::array<int, 9> geq;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];  // geometry: mesh node coordinate
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                // Translational DOF: independent (wall) trans_dof if given, else share the mesh node.
                const int tx = pe.trans_dof[2 * k + 0], ty = pe.trans_dof[2 * k + 1];
                geq[3 * k + 0] = dofs.equation(tx >= 0 ? tx : dofs.global_dof(pe.nodes[k], 0));
                geq[3 * k + 1] = dofs.equation(ty >= 0 ? ty : dofs.global_dof(pe.nodes[k], 1));
                geq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
            }
            plate::Dof up = plate::Dof::Zero();
            for (int a = 0; a < 9; ++a)
                if (geq[a] >= 0) up(a) = u_free(geq[a]) + du_free(geq[a]);
            plate::Dof fp;
            plate::ElementMatrix Kbuf;
            const plate::ElementMatrix* Kp = nullptr;
            if (pe.props.plastic()) {
                const size_t off = pli * plate::kPlasticStateSize;
                plate::internal_force_plastic(Xe, pe.props, up, st.plate_c->data() + off,
                                              st.plate_t->data() + off, fp,
                                              build_tangent ? &Kbuf : nullptr);
                if (build_tangent) Kp = &Kbuf;
            } else {
                Kbuf = plate::stiffness(Xe, pe.props);
                fp = Kbuf * up;
                Kp = &Kbuf;
            }
            for (int a = 0; a < 9; ++a) {
                if (geq[a] < 0) continue;
                f_int(geq[a]) += fp(a);
                if (!build_tangent) continue;
                for (int b = 0; b < 9; ++b)
                    if (geq[b] >= 0) builder->add_entry(geq[a], geq[b], (*Kp)(a, b));
            }
        }

        // --- 5-node (quartic) plate contribution (tri15 edge). 15 DOFs; same logic as 3-node.
        for (size_t pli = 0; pli < structures.plates5.size(); ++pli) {
            const auto& pe = structures.plates5[pli];
            plate::NodeCoords5 Xe;
            std::array<int, 15> geq;
            for (int k = 0; k < 5; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                const int tx = pe.trans_dof[2 * k + 0], ty = pe.trans_dof[2 * k + 1];
                geq[3 * k + 0] = dofs.equation(tx >= 0 ? tx : dofs.global_dof(pe.nodes[k], 0));
                geq[3 * k + 1] = dofs.equation(ty >= 0 ? ty : dofs.global_dof(pe.nodes[k], 1));
                geq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
            }
            plate::Dof5 up = plate::Dof5::Zero();
            for (int a = 0; a < 15; ++a)
                if (geq[a] >= 0) up(a) = u_free(geq[a]) + du_free(geq[a]);
            plate::Dof5 fp;
            plate::ElementMatrix5 Kbuf;
            const plate::ElementMatrix5* Kp = nullptr;
            if (pe.props.plastic()) {
                const size_t off = pli * plate::kPlasticStateSize5;
                plate::internal_force_plastic5(Xe, pe.props, up, st.plate5_c->data() + off,
                                               st.plate5_t->data() + off, fp,
                                               build_tangent ? &Kbuf : nullptr);
                if (build_tangent) Kp = &Kbuf;
            } else {
                Kbuf = plate::stiffness5(Xe, pe.props);
                fp = Kbuf * up;
                Kp = &Kbuf;
            }
            for (int a = 0; a < 15; ++a) {
                if (geq[a] < 0) continue;
                f_int(geq[a]) += fp(a);
                if (!build_tangent) continue;
                for (int b = 0; b < 15; ++b)
                    if (geq[b] >= 0) builder->add_entry(geq[a], geq[b], (*Kp)(a, b));
            }
        }

        // --- Anchor (one-directional axial spring) contribution. N elastic (EA/L)·U,
        // U=(u_b−u_a)·dir; translational DOFs shared. fixed-end (node_b<0): far end fixed ⇒
        // only node_a (2 DOFs). ELASTOPLASTIC (§7a): 1D return mapping,
        // N∈[−Fmax_comp,+Fmax_tens]; at yield permanent U_p, tangent 0. Up is carried
        // committed/trial; the return map is a pure function of U_p^c (line-search safe).
        for (size_t ai = 0; ai < structures.anchors.size(); ++ai) {
            const auto& an = structures.anchors[ai];
            const Eigen::Vector2d Xa(mesh.x[an.node_a], mesh.y[an.node_a]);
            const Eigen::Vector2d Xb = an.node_b >= 0
                ? Eigen::Vector2d(mesh.x[an.node_b], mesh.y[an.node_b]) : an.fixed_point;
            const Eigen::Vector2d dvec = Xb - Xa;
            const double Lgeom = dvec.norm();
            if (Lgeom < 1e-30) continue;
            const Eigen::Vector2d dir = dvec / Lgeom;
            const double kk = an.EA / (an.L > 0.0 ? an.L : Lgeom);
            const int idx[4] = {dofs.equation(dofs.global_dof(an.node_a, 0)),
                                dofs.equation(dofs.global_dof(an.node_a, 1)),
                                an.node_b >= 0 ? dofs.equation(dofs.global_dof(an.node_b, 0)) : -1,
                                an.node_b >= 0 ? dofs.equation(dofs.global_dof(an.node_b, 1)) : -1};
            const double g[4] = {-dir(0), -dir(1), dir(0), dir(1)};  // ∂U/∂u
            double U = 0.0;
            for (int i = 0; i < 4; ++i)
                if (idx[i] >= 0) U += g[i] * (u_free(idx[i]) + du_free(idx[i]));
            const double Up_c = (*st.anchor_c)[ai];
            // The lock-off force rides on the elastic response: N = N0 + k(U - U_p). It is a
            // constant, so it enters the residual and not the tangent -- a prestressed anchor is
            // no stiffer than a slack one, it merely starts loaded. The capacity is checked on
            // the TOTAL force, which is why the plastic elongation below subtracts N0 as well:
            // yielding must leave N exactly at the cap, not at the cap plus the prestress.
            const double N0 = an.prestress;
            double N = N0 + kk * (U - Up_c), Dt = kk;
            const double Ft = an.Fmax_tens, Fc = an.Fmax_comp;
            if (Ft > 0.0 && N > Ft)        { N = Ft;  (*st.anchor_t)[ai] = U - (Ft - N0) / kk; Dt = 0.0; }
            else if (Fc > 0.0 && N < -Fc)  { N = -Fc; (*st.anchor_t)[ai] = U + (Fc + N0) / kk; Dt = 0.0; }
            else                           { (*st.anchor_t)[ai] = Up_c; }
            for (int i = 0; i < 4; ++i) {
                if (idx[i] < 0) continue;
                f_int(idx[i]) += N * g[i];
                if (!build_tangent) continue;
                for (int j = 0; j < 4; ++j)
                    if (idx[j] >= 0) builder->add_entry(idx[i], idx[j], Dt * g[i] * g[j]);
            }
        }

        // --- Geogrid (tension-only axial membrane; §8) contribution. 6 DOFs (translations;
        // no bending/rotation), shared with the soil. Per-Gauss-point tension-only + N_p
        // return mapping: ε=Be·u → N=clamp(EA(ε−ε_p),0,N_p); the compression cut-off is
        // reversible, N_p yield is permanent (ε_p).
        const auto gxi = geogrid::gauss_xi();
        for (size_t gi = 0; gi < structures.geogrids.size(); ++gi) {
            const auto& ge = structures.geogrids[gi];
            geogrid::NodeCoords Xe;
            int geq[6];
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[ge.nodes[k]];
                Xe(k, 1) = mesh.y[ge.nodes[k]];
                geq[2 * k + 0] = dofs.equation(dofs.global_dof(ge.nodes[k], 0));
                geq[2 * k + 1] = dofs.equation(dofs.global_dof(ge.nodes[k], 1));
            }
            geogrid::Dof ug = geogrid::Dof::Zero();
            for (int a = 0; a < 6; ++a)
                if (geq[a] >= 0) ug(a) = u_free(geq[a]) + du_free(geq[a]);
            for (int q = 0; q < geogrid::kGaussCount; ++q) {
                const auto kin = geogrid::axial_kin(Xe, gxi[q]);
                const double eps = (kin.Be * ug)(0);
                const size_t si = gi * geogrid::kGaussCount + q;
                const auto ret = geogrid::axial_return(ge.props, eps, (*st.geogrid_c)[si]);
                (*st.geogrid_t)[si] = ret.ep_new;
                const double w = kin.J;  // 2-point weight = 1
                for (int a = 0; a < 6; ++a) {
                    if (geq[a] < 0) continue;
                    f_int(geq[a]) += w * kin.Be(0, a) * ret.N;
                    if (!build_tangent) continue;
                    for (int b = 0; b < 6; ++b)
                        if (geq[b] >= 0)
                            builder->add_entry(geq[a], geq[b], w * kin.Be(0, a) * ret.Dt * kin.Be(0, b));
                }
            }
        }

        // --- Interface (zero-thickness soil-structure Coulomb interface; see
        // interface-formulation.md). Newton-Cotes (nodal, Day & Potts) integration ⇒ node
        // pairs decouple. At each pair: relative displacement [[u]]=u_struct−u_soil → local
        // (Δu_s,Δu_n) → Coulomb return mapping → 4-DOF block (soil x,y + structure x,y).
        // B(2×4): s-row=[−c,−s,+c,+s], n-row=[+s,−c,−s,+c].
        const auto ncpts = iface::nc_points();
        for (size_t ii = 0; ii < structures.interfaces.size(); ++ii) {
            const auto& ie = structures.interfaces[ii];
            iface::NodeCoords Xe;
            for (int k = 0; k < 3; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount; ++q) {
                const int nd = ncpts[q].node;  // local node index [A,B,middle]
                const auto fr = iface::edge_frame(Xe, ncpts[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts[q].w * fr.J;
                // Equation indices of the 4 DOFs: soil x,y (base) + structure x,y (extra DOFs).
                const int idx[4] = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                    dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                    dofs.equation(ie.struct_dof[2 * nd + 0]),
                                    dofs.equation(ie.struct_dof[2 * nd + 1])};
                const double a[4] = {-c, -s, c, s};   // ∂Δu_s/∂dof
                const double b[4] = {s, -c, -s, c};   // ∂Δu_n/∂dof
                double du_s = 0.0, du_n = 0.0;
                for (int i = 0; i < 4; ++i)
                    if (idx[i] >= 0) {
                        const double ui = u_free(idx[i]) + du_free(idx[i]);
                        du_s += a[i] * ui; du_n += b[i] * ui;
                    }
                const size_t si = ii * iface::kPointCount + q;
                const auto ret = iface::coulomb_return(ie.props, du_s, du_n, (*st.iface_c)[si],
                                                       ie.sigma_n0[q]);
                (*st.iface_t)[si] = ret.slip_p_new;
                for (int i = 0; i < 4; ++i) {
                    if (idx[i] < 0) continue;
                    f_int(idx[i]) += wJ * (a[i] * ret.tau + b[i] * ret.sigma_n);
                    if (!build_tangent) continue;
                    for (int j = 0; j < 4; ++j)
                        if (idx[j] >= 0)
                            builder->add_entry(idx[i], idx[j],
                                               wJ * (ret.Ds * a[i] * a[j] + ret.Dn * b[i] * b[j]));
                }
            }
        }

        // --- 5-node interface (tri15 edge, 5-point Newton-Cotes). Same logic as 3-node.
        const auto ncpts5 = iface::nc_points5();
        for (size_t ii = 0; ii < structures.interfaces5.size(); ++ii) {
            const auto& ie = structures.interfaces5[ii];
            iface::NodeCoords5 Xe;
            for (int k = 0; k < 5; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount5; ++q) {
                const int nd = ncpts5[q].node;
                const auto fr = iface::edge_frame5(Xe, ncpts5[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts5[q].w * fr.J;
                const int idx[4] = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                    dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                    dofs.equation(ie.struct_dof[2 * nd + 0]),
                                    dofs.equation(ie.struct_dof[2 * nd + 1])};
                const double a[4] = {-c, -s, c, s};
                const double b[4] = {s, -c, -s, c};
                double du_s = 0.0, du_n = 0.0;
                for (int i = 0; i < 4; ++i)
                    if (idx[i] >= 0) {
                        const double ui = u_free(idx[i]) + du_free(idx[i]);
                        du_s += a[i] * ui; du_n += b[i] * ui;
                    }
                const size_t si = ii * iface::kPointCount5 + q;
                const auto ret = iface::coulomb_return(ie.props, du_s, du_n, (*st.iface5_c)[si],
                                                       ie.sigma_n0[q]);
                (*st.iface5_t)[si] = ret.slip_p_new;
                for (int i = 0; i < 4; ++i) {
                    if (idx[i] < 0) continue;
                    f_int(idx[i]) += wJ * (a[i] * ret.tau + b[i] * ret.sigma_n);
                    if (!build_tangent) continue;
                    for (int j = 0; j < 4; ++j)
                        if (idx[j] >= 0)
                            builder->add_entry(idx[i], idx[j],
                                               wJ * (ret.Ds * a[i] * a[j] + ret.Dn * b[i] * b[j]));
                }
            }
        }

        // --- Embedded beam (pile row): beam (Timoshenko) stiffness + skin coupling (beam ↔
        // soil via N_s, mesh-nonconforming). Sci.Man §7.5. The beam lives on extra DOFs; the
        // skin splits each point with [N_b,−N_s]. Since E (the soil element) is known in this
        // assembler, N_s = E::shape_functions.
        // Shared helper: axial return mapping + scatter at one coupling point (eq,cx,cy lists).
        auto axial_couple = [&](const int* eqp, const double* cxp, const double* cyp, int nc,
                                const Eigen::Vector2d& tang, double k_a, double k_n, double cap,
                                double wJ, double slip_c, double& slip_t) {
            const Eigen::Vector2d nrm(-tang(1), tang(0));
            Eigen::Vector2d dur(0.0, 0.0);
            for (int d = 0; d < nc; ++d)
                if (eqp[d] >= 0) { const double ud = u_free(eqp[d]) + du_free(eqp[d]); dur(0) += cxp[d] * ud; dur(1) += cyp[d] * ud; }
            const double dua = dur.dot(tang), dun = dur.dot(nrm);
            double ta = k_a * (dua - slip_c), Da = k_a;
            if (cap > 0.0 && std::fabs(ta) > cap) { ta = std::copysign(cap, ta); slip_t = dua - ta / k_a; Da = 0.0; }
            else { slip_t = slip_c; }
            const double tn = k_n * dun;
            const Eigen::Vector2d tr = ta * tang + tn * nrm;
            const Eigen::Matrix2d D = Da * (tang * tang.transpose()) + k_n * (nrm * nrm.transpose());
            for (int d = 0; d < nc; ++d) {
                if (eqp[d] < 0) continue;
                f_int(eqp[d]) += wJ * (cxp[d] * tr(0) + cyp[d] * tr(1));
                if (!build_tangent) continue;
                const Eigen::RowVector2d cd(cxp[d], cyp[d]);
                for (int e2 = 0; e2 < nc; ++e2)
                    if (eqp[e2] >= 0) builder->add_entry(eqp[d], eqp[e2],
                                          wJ * (cd * D * Eigen::Vector2d(cxp[e2], cyp[e2]))(0, 0));
            }
        };

        size_t skin_off = 0;
        for (size_t bi = 0; bi < structures.embedded_beams.size(); ++bi) {
            const auto& eb = structures.embedded_beams[bi];
            for (const auto& el : eb.elements) {  // (1) beam (Timoshenko) stiffness
                plate::NodeCoords Xe;
                std::array<int, 9> geq;
                for (int k = 0; k < 3; ++k) {
                    Xe(k, 0) = eb.node_x[el[k]]; Xe(k, 1) = eb.node_y[el[k]];
                    geq[3 * k + 0] = ebeam::trans_eq(eb, el[k], 0, dofs);
                    geq[3 * k + 1] = ebeam::trans_eq(eb, el[k], 1, dofs);
                    geq[3 * k + 2] = dofs.equation(eb.dof_phi[el[k]]);
                }
                const plate::ElementMatrix Kp = plate::stiffness(Xe, eb.props);
                plate::Dof up = plate::Dof::Zero();
                for (int a = 0; a < 9; ++a) if (geq[a] >= 0) up(a) = u_free(geq[a]) + du_free(geq[a]);
                const plate::Dof fp = Kp * up;
                for (int a = 0; a < 9; ++a) {
                    if (geq[a] < 0) continue;
                    f_int(geq[a]) += fp(a);
                    if (!build_tangent) continue;
                    for (int b = 0; b < 9; ++b)
                        if (geq[b] >= 0) builder->add_entry(geq[a], geq[b], Kp(a, b));
                }
            }
            constexpr int NC = 6 + 2 * E::kNodeCount;
            for (size_t pi = 0; pi < eb.skin.size(); ++pi) {  // (2) skin (axial cap T_max)
                const auto& sp = eb.skin[pi];
                if (!sp.ok) continue;
                const auto Ns = E::shape_functions(sp.xi_s, sp.eta_s);
                std::array<int, NC> eq; std::array<double, NC> cx, cy; int nc = 0;
                for (int i = 0; i < 3; ++i) {
                    eq[nc] = ebeam::trans_eq(eb, sp.beam_node[i], 0, dofs); cx[nc] = sp.Nb(i); cy[nc] = 0.0; ++nc;
                    eq[nc] = ebeam::trans_eq(eb, sp.beam_node[i], 1, dofs); cx[nc] = 0.0; cy[nc] = sp.Nb(i); ++nc;
                }
                for (int j = 0; j < E::kNodeCount; ++j) {
                    const int sn = mesh.node_of(sp.soil_elem, j);
                    eq[nc] = dofs.equation(dofs.global_dof(sn, 0)); cx[nc] = -Ns(j); cy[nc] = 0.0; ++nc;
                    eq[nc] = dofs.equation(dofs.global_dof(sn, 1)); cx[nc] = 0.0; cy[nc] = -Ns(j); ++nc;
                }
                axial_couple(eq.data(), cx.data(), cy.data(), nc, sp.tang, sp.k_a, sp.k_n, sp.t_max,
                             sp.wJ, (*st.eskin_c)[skin_off + pi], (*st.eskin_t)[skin_off + pi]);
            }
            skin_off += eb.skin.size();
            if (eb.foot.D_foot > 0.0 && eb.foot.ok) {  // (3) foot (axial spring, cap F_max; wJ=1)
                const auto Ns = E::shape_functions(eb.foot.xi_s, eb.foot.eta_s);
                constexpr int NF = 2 + 2 * E::kNodeCount;
                std::array<int, NF> eq; std::array<double, NF> cx, cy; int nc = 0;
                eq[nc] = ebeam::trans_eq(eb, eb.foot.beam_node, 0, dofs); cx[nc] = 1.0; cy[nc] = 0.0; ++nc;
                eq[nc] = ebeam::trans_eq(eb, eb.foot.beam_node, 1, dofs); cx[nc] = 0.0; cy[nc] = 1.0; ++nc;
                for (int j = 0; j < E::kNodeCount; ++j) {
                    const int sn = mesh.node_of(eb.foot.soil_elem, j);
                    eq[nc] = dofs.equation(dofs.global_dof(sn, 0)); cx[nc] = -Ns(j); cy[nc] = 0.0; ++nc;
                    eq[nc] = dofs.equation(dofs.global_dof(sn, 1)); cx[nc] = 0.0; cy[nc] = -Ns(j); ++nc;
                }
                axial_couple(eq.data(), cx.data(), cy.data(), nc, eb.foot.tang, eb.foot.D_foot, 0.0,
                             eb.foot.f_max, 1.0, (*st.efoot_c)[bi], (*st.efoot_t)[bi]);
            }
        }
        if (timings) {
            if (build_tangent) { timings->assemble_tangent += elapsed(t_asm); ++timings->n_tangent; }
            else               { timings->assemble_residual += elapsed(t_asm); ++timings->n_residual; }
        }
        return f_int;
    }

private:
    const mesh::Mesh& mesh_;
    const DofMap& dofs_;
    const std::vector<MaterialModel>& materials_;
    const std::vector<char>& active_;
    const Structures& structures_;
    const std::vector<MaterialProfile>& profile_;
    int neq_;
    std::vector<typename E::NodeCoords> all_coords_;
    std::vector<std::array<int, E::kDofCount>> all_edofs_;
    std::vector<ElementVector> fe_buf_;
    std::vector<typename E::ElementMatrix> ke_buf_;
};

}  // namespace detail
}  // namespace katai::core
