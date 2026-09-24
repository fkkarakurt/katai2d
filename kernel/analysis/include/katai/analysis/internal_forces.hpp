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
    using Report = PointReportT<Tangent>;
    static void integrate(const MaterialModel& m, const GaussState& comm,
                          const Strain& de, GaussState& tr, Tangent& t, TangentMode mode,
                          double creep_dt, Report* rep = nullptr, double substep_tol = 0.0) {
        integrate_point(m, comm, de, tr, t, mode, creep_dt, rep, substep_tol);
    }
    static Eigen::Matrix<double, 3, 1> stress(const GaussState& s) { return s.stress; }
    // The stress with the OUT-OF-PLANE component put back: [xx, yy, xy, zz]. Plane strain
    // carries sigma_zz outside the 3-vector the element assembles with, so a local error
    // measured on the in-plane block alone would ignore a principal stress -- and under K0 it
    // is regularly the extreme one.
    static Eigen::Vector4d full_stress(const GaussState& s) {
        Eigen::Vector4d v;
        v << s.stress, s.stress_zz;
        return v;
    }
    // D^e acting on a strain increment, in the same four components. The zz row of an
    // ISOTROPIC elastic operator is lambda*(de_xx + de_yy), and lambda is the operator's own
    // off-diagonal -- so the fourth component is read out of D^e rather than rebuilt from an
    // (E, nu) pair the model may not even use.
    static Eigen::Vector4d elastic_step(const Tangent& De, const Strain& de) {
        Eigen::Vector4d v;
        v << De * de, De(0, 1) * (de(0) + de(1));
        return v;
    }
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
    using Report = PointReportT<Tangent>;
    static void integrate(const MaterialModel& m, const GaussState& comm,
                          const Strain& de, GaussState& tr, Tangent& t, TangentMode mode,
                          double creep_dt, Report* rep = nullptr, double substep_tol = 0.0) {
        integrate_point_axisym(m, comm, de, tr, t, mode, creep_dt, rep, substep_tol);
    }
    static Eigen::Matrix<double, 4, 1> stress(const GaussState& s) {
        Eigen::Matrix<double, 4, 1> v;
        v << s.stress, s.stress_zz;
        return v;
    }
    // Axisymmetry already carries all four components, hoop included.
    static Eigen::Vector4d full_stress(const GaussState& s) { return stress(s); }
    static Eigen::Vector4d elastic_step(const Tangent& De, const Strain& de) { return De * de; }
    // Pore-pressure direction m = [r, z, rz, theta]: the hoop is a real normal stress.
    static Strain pore_vector() {
        Strain m;
        m << 1.0, 1.0, 0.0, 1.0;
        return m;
    }
};

// One previous-iterate record for a STRUCTURAL stress point: what the law returned there, and
// the relative displacement it was returned for. The soil keeps the same pair in two parallel
// vectors because its stress is a tensor; a traction is one number, so it fits in a struct.
struct PrevPoint {
    double value = 0.0;   // tau (interface) or axial force (coupling spring / foot)
    double du = 0.0;      // the relative displacement that produced it
};

// LOCAL convergence measurement, gathered while the internal force is assembled because that
// is the only place where both stresses of a point exist at the same instant.
//
// The concept: a stress point carries TWO stresses per iteration. The CONSTITUTIVE stress is
// what the material law returns for the strain the point was given,
// sigma_c,j = sigma_0 + D^e (Delta-eps_j - Delta-eps_p,j). The EQUILIBRIUM stress is
// what the finite-element linearisation says the point carries, sigma_eq,j = sigma_c,j-1 +
// D^e delta-eps_j, built from the previous iterate with this iteration's displacement
// correction. They coincide only at convergence. A run that satisfies global equilibrium while
// those two disagree is a run whose nodal forces balance around stresses the material law
// would not produce -- and nothing in a global force residual can see it.
//
// The two vectors below are the iteration's memory: the caller hands the same pair back on
// every iterate and they are overwritten here with this iterate's values, so no other part of
// the solver has to know where the iteration was.
struct LocalErrorProbe {
    std::vector<GaussState>* prev_sigma_c = nullptr;  // sigma_c,j-1 in, sigma_c,j out
    std::vector<double>* prev_deps = nullptr;         // Delta-eps_{j-1}, kStrain per point
    // The structural counterparts, one entry per stress point of each kind, sized like the
    // committed-state vectors they shadow.
    std::vector<PrevPoint>* prev_iface = nullptr;     // 3-node interfaces
    std::vector<PrevPoint>* prev_iface5 = nullptr;    // 5-node interfaces
    std::vector<PrevPoint>* prev_skin = nullptr;      // embedded-beam skin coupling springs
    std::vector<PrevPoint>* prev_foot = nullptr;      // embedded-beam feet
    // The first iterate of an increment has applied no correction, so there is no local error to
    // measure there -- only a state to record. The soil gets this for free (its previous stress
    // is the committed one and its previous strain increment is zero); a traction has no
    // committed value stored anywhere, so it is recorded on the way past instead.
    bool first_iterate = false;
    double tolerated = 0.01;                          // ToleratedError, the bar for a local error

    // --- measured, over the ACTIVE soil stress points -----------------------------------
    int plastic = 0;                 // points with a yield or cap surface active
    int plastic_inaccurate = 0;      // ... of those, over the tolerated local error
    int elastic_total = 0;           // points with no surface active
    int nl_elastic = 0;              // ... of those, with a STRESS-DEPENDENT elastic stiffness
    int nl_elastic_inaccurate = 0;   // ... of those, over the tolerated local error
    double worst_plastic = 0.0;      // the largest local error of each kind, so that a count
    double worst_nl_elastic = 0.0;   //   of zero can still say how much room it had
    // Points whose CONSTITUTIVE INTEGRATION hit its substep guard on this iterate, i.e. walked
    // the material law without meeting the error tolerance it was given. Not a local convergence
    // criterion and not derivable from one: the equilibrium and constitutive stresses can agree
    // perfectly while both of them are wrong, because they are both computed from the same
    // cut-short walk.
    int integration_saturated = 0;

    // The current stiffness parameter (CSP) numerator and denominator, integrated over the
    // active volume: the work actually done against the stress increment, over the work the same
    // strain would have done had the response stayed elastic.
    double energy_total = 0.0;       // integral of Delta-eps . Delta-sigma
    double energy_elastic = 0.0;     // integral of Delta-eps . D^e Delta-eps

    // Moment criterion: the largest out-of-balance moment on a rotational equation, over m_ref.
    // m_ref is summed HERE, element by element, because it is a sum of ABSOLUTE moment
    // contributions: taken after assembly, the opposing contributions of two adjacent plate
    // elements cancel and the reference collapses towards zero exactly where the structure is in
    // equilibrium -- which is where the criterion has to work.
    double m_ref = 0.0;
    std::vector<int> moment_eq;      // equations carrying a rotational degree of freedom

    // Interfaces: the local error at a plastic interface point (see measure_traction_point). The
    // embedded beam's SKIN coupling springs are counted in the same tally as the soil-structure
    // interfaces -- no distinction is drawn between a standard interface and the special one a
    // pile skin is -- so they are counted here together and the record says so.
    int iface_plastic = 0, iface_plastic_inaccurate = 0;
    double worst_iface = 0.0;

    // Embedded-beam foot force error: the difference between equilibrium and constitutive foot
    // forces, relative to the foot forces themselves. Three sums rather than a count, because the
    // criterion is one ratio over all the feet in the model, not a per-point test.
    int feet = 0;
    double foot_num = 0.0;        // sum |F_foot,eq - F_foot,c|
    double foot_den_c = 0.0;      // sum |F_foot,c|
    double foot_den_max = 0.0;    // sum |F_foot,max|

    // Clear what one iterate measured, keeping what carries the iteration across iterates
    // (the previous-iterate vectors and the tolerance).
    void reset_counts() {
        plastic = plastic_inaccurate = 0;
        elastic_total = nl_elastic = nl_elastic_inaccurate = 0;
        worst_plastic = worst_nl_elastic = 0.0;
        integration_saturated = 0;
        energy_total = energy_elastic = 0.0;
        iface_plastic = iface_plastic_inaccurate = 0;
        worst_iface = 0.0;
        feet = 0;
        foot_num = foot_den_c = foot_den_max = 0.0;
    }
};

// One interface / coupling-spring stress point, measured the way a soil stress point is measured
// (LocalErrorProbe): the EQUILIBRIUM traction is the previous iterate's constitutive traction
// carried forward elastically through this iteration's slip increment, and the CONSTITUTIVE
// traction is what the Coulomb return actually gave. The local interface error normalises their
// difference by the point's own shear capacity, max(tau_max, c), with the same 1 kPa floor as
// the soil so that a point carrying nothing cannot report
// an enormous relative error on a difference that is numerically nothing.
//
// `prev` is read and then overwritten with this iterate's pair, which is what carries the
// measurement to the next iteration.
inline void measure_traction_point(LocalErrorProbe& probe, PrevPoint& prev, double tau_c,
                                   double du, double k_elastic, bool plastic, double tau_max,
                                   double cohesion) {
    if (!probe.first_iterate) {
        const double tau_eq = prev.value + k_elastic * (du - prev.du);
        const double denom = std::max(std::max(tau_max, cohesion), 1.0);
        const double err = std::fabs(tau_eq - tau_c) / denom;
        if (plastic) {
            ++probe.iface_plastic;
            if (err > probe.tolerated) ++probe.iface_plastic_inaccurate;
        }
        probe.worst_iface = std::max(probe.worst_iface, err);
    } else if (plastic) {
        ++probe.iface_plastic;
    }
    prev.value = tau_c;
    prev.du = du;
}

// Per-element partial of the above. The assembly runs in parallel over elements, and a global
// counter written from several threads is both a race and a source of run-to-run variation in
// the floating-point sums -- this tree's results are deterministic and stay that way: each
// element fills its own partial, and the reduction happens in the sequential scatter pass, in
// element order.
struct LocalErrorPartial {
    int plastic = 0, plastic_inaccurate = 0;
    int elastic_total = 0, nl_elastic = 0, nl_elastic_inaccurate = 0;
    double worst_plastic = 0.0, worst_nl_elastic = 0.0;
    int integration_saturated = 0;
    double energy_total = 0.0, energy_elastic = 0.0;
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
        // Local convergence measurement. nullptr (the default) = not measured, and then not a
        // single extra flop happens in the Gauss loop: the line search re-evaluates the
        // residual several times per iteration and none of those iterates is the one whose
        // local error means anything.
        LocalErrorProbe* local = nullptr;
    };

    // Nonzero-Dirichlet (prescribed displacement ū) ramp — used ONLY by the static path.
    // presc==null → no fixed-DOF displacement at all (the dynamic seismic path, in the
    // relative frame, does not use this). factor = (cur_target − cur_lambda): this step's
    // increment fraction.
    // CLOSED 2026-08-13. This block used to declare a hidden limit: the ramp entered only the
    // SOIL element loop, so a plate, geogrid, interface or embedded beam standing on a boundary
    // driven by a prescribed displacement did not see that motion at all (geq < 0 => the node was
    // read as ZERO). The declaration rested on a guard that had since expired -- "prescribed
    // u_bar is today a kernel/test path, not available as a deformation BC in the GUI" -- because
    // schema v2 added `disps`, which made it reachable from a .k2d file and from the GUI while
    // the limit stayed declared and unrevisited. A declared limit is worth exactly the sentence
    // that scopes it, and that sentence had stopped being true.
    //
    // Measured on the geogrid case (KV-STR-005) before the fix: a reinforcement pulled through a
    // block by a prescribed edge reported its own axial force as 7.6363 kN/m where EA*eps is
    // 2.0000, and the converged field was not affine when it had to be -- the node beside the
    // driven edge sat 3.79% BELOW its own left neighbour, because the element next to it believed
    // the driven node had never moved.
    //
    // The fix is the one the old comment prescribed: structural elements are TOTAL-displacement
    // formulations, so they need the total lambda*u_bar applied so far, not the increment the
    // soil loop uses. `factor` is this step's increment (cur_target - cur_lambda) and drives the
    // soil; `total` is cur_target and drives the structures.
    struct Ramp {
        const Eigen::VectorXd* presc = nullptr;   // total_dofs
        double factor = 0.0;                      // this step's INCREMENT fraction (soil)
        double total = 0.0;                       // fraction applied by the END of this step
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
    // The CONSTITUTIVE integration error tolerance for this solve (.k2d v15 / NewtonOptions).
    // 0 = the material class's own default, which is what every path did before it was writable.
    // It sits here beside dt_day because both are per-solve numerics the Gauss loop reads, not
    // material properties: a tolerance published as a soil parameter is the mistake this tree has
    // already made once.
    double substep_tol = 0.0;

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

        // The TOTAL displacement of one GLOBAL dof, for the total-displacement structural
        // formulations below. A fixed DOF is not simply zero: under a non-zero prescribed
        // displacement it carries the share of u_bar applied so far. Passing the global dof
        // rather than the equation is the whole point -- the equation of a fixed DOF is -1, and
        // that is exactly the information the old code threw away.
        const auto u_at = [&](int gdof) -> double {
            const int eq = dofs.equation(gdof);
            if (eq >= 0) return u_free(eq) + du_free(eq);
            return has_presc ? ramp.total * (*ramp.presc)(gdof) : 0.0;
        };
        // What a structural element reads: the total displacement less the datum of the phase it
        // was installed in (Structures::install_datum), on a free DOF. install < 0 is the zero
        // datum and returns u_at itself -- no subtraction is performed at all, so every element of
        // a chain whose structure set never changed computes exactly what it always did.
        const auto u_el = [&](int gdof, int install) -> double {
            if (install < 0) return u_at(gdof);
            const int eq = dofs.equation(gdof);
            return eq >= 0 ? u_at(gdof) - structures.install_datum[(size_t)install](eq)
                           : u_at(gdof);
        };

        Eigen::VectorXd f_int = Eigen::VectorXd::Zero(neq);

        // --- Soil elements: TWO-PHASE assembly. (1) Per-element computation in PARALLEL —
        // each element writes only to its OWN fe/ke buffer and its OWN Gauss trial states (no
        // data race; no shared accumulator → floating-point order unchanged). (2) Scatter
        // SEQUENTIAL — f_int accumulation and COO entry order are IDENTICAL to the serial
        // path = the result is deterministic, independent of the thread count.
        if (build_tangent && ke_buf_.size() != static_cast<size_t>(mesh.element_count))
            ke_buf_.resize(mesh.element_count);
        LocalErrorProbe* const probe = st.local;
        if (probe) {
            probe_buf_.assign(mesh.element_count, LocalErrorPartial{});
            probe->m_ref = 0.0;
            probe->moment_eq.clear();
            if (probe->prev_deps->size() !=
                static_cast<size_t>(mesh.element_count) * n_gp * Kin::kStrain)
                probe->prev_deps->assign(
                    static_cast<size_t>(mesh.element_count) * n_gp * Kin::kStrain, 0.0);
            // The structural previous-iterate vectors shadow the committed-state vectors one for
            // one, so they are sized FROM them rather than re-derived from the element counts --
            // one place decides how many stress points a structure has.
            const auto fit = [](std::vector<PrevPoint>* v, size_t n) {
                if (v && v->size() != n) v->assign(n, PrevPoint{});
            };
            fit(probe->prev_iface, st.iface_c ? st.iface_c->size() : 0);
            fit(probe->prev_iface5, st.iface5_c ? st.iface5_c->size() : 0);
            fit(probe->prev_skin, st.eskin_c ? st.eskin_c->size() : 0);
            fit(probe->prev_foot, st.efoot_c ? st.efoot_c->size() : 0);
        }
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
            // Hoisted out of the Gauss loop and reset only when it is going to be read: the
            // report holds an elastic matrix, so constructing one per stress point per assembly
            // would zero nine (or sixteen) doubles on the hottest loop in the program for the
            // benefit of a measurement that is off.
            typename Kin::Report rep;
            for (int g = 0; g < n_gp; ++g) {
                const auto grad =
                    Kin::template gradients<E>(coords, gauss[g].xi, gauss[g].eta);
                const typename Kin::Strain dstrain = grad.B * du_e;
                const int gi = e * n_gp + g;
                typename Kin::Tangent dt;
                const TangentMode tm = !build_tangent ? TangentMode::kNone : tmode;
                // Depth-varying stiffness / cohesion are evaluated HERE, at the stress point: the
                // whole point of a gradient is that it varies WITHIN an element,
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
                if (probe) rep = typename Kin::Report{};
                Kin::integrate(matg, committed[gi], dstrain, trial[gi], dt, tm, dt_day,
                               probe ? &rep : nullptr, substep_tol);
                typename Kin::Strain sigma = Kin::stress(trial[gi]);
                if (probe) {
                    LocalErrorPartial& acc = probe_buf_[e];
                    const double w_loc = gauss[g].weight * grad.weight;
                    // --- CSP, the current stiffness parameter: the work this increment actually
                    // did, over the work the same strain would have done had the response stayed
                    // elastic. Unity while elastic, falling towards zero as the body plastifies.
                    const typename Kin::Strain dsig = sigma - Kin::stress(committed[gi]);
                    acc.energy_total += w_loc * dstrain.dot(dsig);
                    acc.energy_elastic += w_loc * dstrain.dot(rep.elastic * dstrain);
                    // --- the two stresses of the point: constitutive (s_c) and equilibrium
                    // (s_eq = sigma_c,j-1 + D^e delta-eps_j).
                    const size_t base = static_cast<size_t>(gi) * Kin::kStrain;
                    typename Kin::Strain deps_prev;
                    for (int k = 0; k < Kin::kStrain; ++k)
                        deps_prev(k) = (*probe->prev_deps)[base + k];
                    const typename Kin::Strain ddeps = dstrain - deps_prev;  // delta-eps_j
                    const Eigen::Vector4d s_c = Kin::full_stress(trial[gi]);
                    const Eigen::Vector4d s_eq =
                        Kin::full_stress((*probe->prev_sigma_c)[gi]) +
                        Kin::elastic_step(rep.elastic, ddeps);
                    if (rep.integration_saturated) ++acc.integration_saturated;
                    const double diff = (s_eq - s_c).norm();
                    const double tmax = tau_max_of(s_c(0), s_c(1), s_c(2), s_c(3));
                    const double coh = cohesion_of(matg);
                    if (rep.plastic) {
                        // The local error at a plastic stress point: |sigma_eq - sigma_c| over
                        // max(tau_max, c, 1 kPa). The 1 kPa floor stops a point that carries
                        // almost no stress at all from reporting an enormous relative error on a
                        // difference that is numerically nothing.
                        const double err = diff / std::max(std::max(tmax, coh), 1.0);
                        ++acc.plastic;
                        if (err > probe->tolerated) ++acc.plastic_inaccurate;
                        acc.worst_plastic = std::max(acc.worst_plastic, err);
                    } else {
                        ++acc.elastic_total;
                        if (rep.stress_dependent) {
                            // The local error at a stress-dependent elastic point. Same
                            // difference, a different floor: a point whose ELASTIC stiffness
                            // moves with stress can be inaccurate while carrying no plasticity
                            // at all, and p_ref/200 is the scale that is measured against.
                            const double err = diff / std::max(std::max(tmax, coh),
                                                               p_ref_of(matg) / 200.0);
                            ++acc.nl_elastic;
                            if (err > probe->tolerated) ++acc.nl_elastic_inaccurate;
                            acc.worst_nl_elastic = std::max(acc.worst_nl_elastic, err);
                        }
                    }
                    (*probe->prev_sigma_c)[gi] = trial[gi];
                    for (int k = 0; k < Kin::kStrain; ++k)
                        (*probe->prev_deps)[base + k] = dstrain(k);
                }
                // Undrained (A): the constitutive model returns EFFECTIVE stress and its effective
                // tangent; add the pore fluid's volumetric (bulk) contribution. Total = sigma' +
                // (Kw/n) eps_v m and tangent += (Kw/n) m m^T = D_u. The excess pore pressure
                // u = -(Kw/n) eps_v is carried in the Gauss state (eps_vol_und, pw_carried) from
                // phase to phase. (See effective-stress-formulation.md.)
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
                // The excess pore pressure is the point's own state (GaussState::eps_vol_und,
                // pw_carried): generated here only while the material is undrained and the phase
                // does not ignore it, and carried -- never dropped, never re-derived from the
                // volume change of an earlier drained stage -- otherwise.
                trial[gi].eps_vol_und = committed[gi].eps_vol_und;
                trial[gi].pw_carried = committed[gi].pw_carried;
                if (matg.undrained) {
                    // Undrained (A)/(B): the constitutive model returned EFFECTIVE stress; add the
                    // pore fluid's volumetric contribution. u = -(Kw/n) eps_v.
                    const double kwn = matg.kw_over_n(matg.undrained_poisson);
                    if (!matg.ignore_undrained) {
                        trial[gi].eps_vol_und = committed[gi].eps_vol_und + mvec.dot(dstrain);
                        sigma += kwn * trial[gi].eps_vol_und * mvec;
                        if (build_tangent) dt += kwn * mvec * mvec.transpose();
                    } else if (trial[gi].eps_vol_und != 0.0) {
                        // Ignoring undrained behaviour: no water stiffness and nothing new, but
                        // the pressure generated before this phase stays where it is.
                        sigma += kwn * trial[gi].eps_vol_und * mvec;
                    }
                }
                if (trial[gi].pw_carried != 0.0) sigma += trial[gi].pw_carried * mvec;
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
            if (probe) {  // reduction in element order -> the sums do not depend on the threads
                const LocalErrorPartial& a = probe_buf_[e];
                probe->plastic += a.plastic;
                probe->plastic_inaccurate += a.plastic_inaccurate;
                probe->elastic_total += a.elastic_total;
                probe->nl_elastic += a.nl_elastic;
                probe->nl_elastic_inaccurate += a.nl_elastic_inaccurate;
                probe->worst_plastic = std::max(probe->worst_plastic, a.worst_plastic);
                probe->worst_nl_elastic = std::max(probe->worst_nl_elastic, a.worst_nl_elastic);
                probe->integration_saturated += a.integration_saturated;
                probe->energy_total += a.energy_total;
                probe->energy_elastic += a.energy_elastic;
            }
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
            std::array<int, 9> geq, gd;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];  // geometry: mesh node coordinate
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                // Translational DOF: independent (wall) trans_dof if given, else share the mesh node.
                const int tx = pe.trans_dof[2 * k + 0], ty = pe.trans_dof[2 * k + 1];
                gd[3 * k + 0] = tx >= 0 ? tx : dofs.global_dof(pe.nodes[k], 0);
                gd[3 * k + 1] = ty >= 0 ? ty : dofs.global_dof(pe.nodes[k], 1);
                gd[3 * k + 2] = pe.rot_dof[k];
                for (int c = 0; c < 3; ++c) geq[3 * k + c] = dofs.equation(gd[3 * k + c]);
            }
            plate::Dof up = plate::Dof::Zero();
            for (int a = 0; a < 9; ++a) up(a) = u_el(gd[a], pe.install);
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
            if (probe)
                for (int k = 0; k < 3; ++k) {
                    probe->m_ref += std::fabs(fp(3 * k + 2));
                    if (geq[3 * k + 2] >= 0) probe->moment_eq.push_back(geq[3 * k + 2]);
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
            std::array<int, 15> geq, gd;
            for (int k = 0; k < 5; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                const int tx = pe.trans_dof[2 * k + 0], ty = pe.trans_dof[2 * k + 1];
                gd[3 * k + 0] = tx >= 0 ? tx : dofs.global_dof(pe.nodes[k], 0);
                gd[3 * k + 1] = ty >= 0 ? ty : dofs.global_dof(pe.nodes[k], 1);
                gd[3 * k + 2] = pe.rot_dof[k];
                for (int c = 0; c < 3; ++c) geq[3 * k + c] = dofs.equation(gd[3 * k + c]);
            }
            plate::Dof5 up = plate::Dof5::Zero();
            for (int a = 0; a < 15; ++a) up(a) = u_el(gd[a], pe.install);
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
            if (probe)
                for (int k = 0; k < 5; ++k) {
                    probe->m_ref += std::fabs(fp(3 * k + 2));
                    if (geq[3 * k + 2] >= 0) probe->moment_eq.push_back(geq[3 * k + 2]);
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
            const int gdx[4] = {anchor_dof(an, dofs, 0), anchor_dof(an, dofs, 1),
                                anchor_dof(an, dofs, 2), anchor_dof(an, dofs, 3)};
            const int idx[4] = {dofs.equation(gdx[0]), dofs.equation(gdx[1]),
                                gdx[2] >= 0 ? dofs.equation(gdx[2]) : -1,
                                gdx[3] >= 0 ? dofs.equation(gdx[3]) : -1};
            const double g[4] = {-dir(0), -dir(1), dir(0), dir(1)};  // ∂U/∂u
            double U = 0.0;
            for (int i = 0; i < 4; ++i)
                if (gdx[i] >= 0) U += g[i] * u_el(gdx[i], an.install);
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
            int geq[6], gd[6];
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[ge.nodes[k]];
                Xe(k, 1) = mesh.y[ge.nodes[k]];
                gd[2 * k + 0] = geogrid_dof(ge, dofs, k, 0);
                gd[2 * k + 1] = geogrid_dof(ge, dofs, k, 1);
                geq[2 * k + 0] = dofs.equation(gd[2 * k + 0]);
                geq[2 * k + 1] = dofs.equation(gd[2 * k + 1]);
            }
            geogrid::Dof ug = geogrid::Dof::Zero();
            for (int a = 0; a < 6; ++a) ug(a) = u_el(gd[a], ge.install);
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
            if (!ie.active) {   // the ground on one side is not there: the joint carries nothing
                for (int q = 0; q < iface::kPointCount; ++q)
                    (*st.iface_t)[ii * iface::kPointCount + q] = (*st.iface_c)[ii * iface::kPointCount + q];
                continue;
            }
            iface::NodeCoords Xe;
            for (int k = 0; k < 3; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount; ++q) {
                const int nd = ncpts[q].node;  // local node index [A,B,middle]
                const auto fr = iface::edge_frame(Xe, ncpts[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts[q].w * fr.J;
                // Equation indices of the 4 DOFs: soil x,y (base) + structure x,y (extra DOFs).
                const int gdx[4] = {dofs.global_dof(ie.soil_nodes[nd], 0),
                                    dofs.global_dof(ie.soil_nodes[nd], 1),
                                    ie.struct_dof[2 * nd + 0], ie.struct_dof[2 * nd + 1]};
                const int idx[4] = {dofs.equation(gdx[0]), dofs.equation(gdx[1]),
                                    dofs.equation(gdx[2]), dofs.equation(gdx[3])};
                const double a[4] = {-c, -s, c, s};   // ∂Δu_s/∂dof
                const double b[4] = {s, -c, -s, c};   // ∂Δu_n/∂dof
                double du_s = 0.0, du_n = 0.0;
                for (int i = 0; i < 4; ++i) {
                    const double ui = u_el(gdx[i], ie.install);
                    du_s += a[i] * ui; du_n += b[i] * ui;
                }
                const size_t si = ii * iface::kPointCount + q;
                const auto ret = iface::coulomb_return(ie.props, du_s, du_n, (*st.iface_c)[si],
                                                       ie.sigma_n0[q]);
                (*st.iface_t)[si] = ret.slip_p_new;
                if (probe && probe->prev_iface)
                    measure_traction_point(*probe, (*probe->prev_iface)[si], ret.tau, du_s,
                                           ie.props.ks, ret.Ds == 0.0,
                                           std::max(0.0, ie.props.c_i -
                                                             ret.sigma_n * std::tan(ie.props.phi_i)),
                                           ie.props.c_i);
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
            if (!ie.active) {
                for (int q = 0; q < iface::kPointCount5; ++q)
                    (*st.iface5_t)[ii * iface::kPointCount5 + q] = (*st.iface5_c)[ii * iface::kPointCount5 + q];
                continue;
            }
            iface::NodeCoords5 Xe;
            for (int k = 0; k < 5; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount5; ++q) {
                const int nd = ncpts5[q].node;
                const auto fr = iface::edge_frame5(Xe, ncpts5[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts5[q].w * fr.J;
                const int gdx[4] = {dofs.global_dof(ie.soil_nodes[nd], 0),
                                    dofs.global_dof(ie.soil_nodes[nd], 1),
                                    ie.struct_dof[2 * nd + 0], ie.struct_dof[2 * nd + 1]};
                const int idx[4] = {dofs.equation(gdx[0]), dofs.equation(gdx[1]),
                                    dofs.equation(gdx[2]), dofs.equation(gdx[3])};
                const double a[4] = {-c, -s, c, s};
                const double b[4] = {s, -c, -s, c};
                double du_s = 0.0, du_n = 0.0;
                for (int i = 0; i < 4; ++i) {
                    const double ui = u_el(gdx[i], ie.install);
                    du_s += a[i] * ui; du_n += b[i] * ui;
                }
                const size_t si = ii * iface::kPointCount5 + q;
                const auto ret = iface::coulomb_return(ie.props, du_s, du_n, (*st.iface5_c)[si],
                                                       ie.sigma_n0[q]);
                (*st.iface5_t)[si] = ret.slip_p_new;
                if (probe && probe->prev_iface5)
                    measure_traction_point(*probe, (*probe->prev_iface5)[si], ret.tau, du_s,
                                           ie.props.ks, ret.Ds == 0.0,
                                           std::max(0.0, ie.props.c_i -
                                                             ret.sigma_n * std::tan(ie.props.phi_i)),
                                           ie.props.c_i);
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
        // soil via N_s, mesh-nonconforming). The beam lives on extra DOFs; the
        // skin splits each point with [N_b,−N_s]. Since E (the soil element) is known in this
        // assembler, N_s = E::shape_functions.
        // Shared helper: axial return mapping + scatter at one coupling point (eq,cx,cy lists).
        // out_f / out_du hand back the axial force the law returned and the relative axial
        // displacement it was returned for -- the pair the local criteria need (the local
        // interface error for the skin springs, the foot force error for the foot). Both default
        // to null, so the measurement costs
        // nothing where it is not asked for.
        auto axial_couple = [&](const int* eqp, const int* gdp, const double* cxp,
                                const double* cyp, int nc, int install,
                                const Eigen::Vector2d& tang, double k_a, double k_n, double cap,
                                double wJ, double slip_c, double& slip_t,
                                double* out_f = nullptr, double* out_du = nullptr,
                                bool* out_capped = nullptr, bool no_tension = false) {
            const Eigen::Vector2d nrm(-tang(1), tang(0));
            Eigen::Vector2d dur(0.0, 0.0);
            for (int d = 0; d < nc; ++d)
                if (gdp[d] >= 0) { const double ud = u_el(gdp[d], install); dur(0) += cxp[d] * ud; dur(1) += cyp[d] * ud; }
            const double dua = dur.dot(tang), dun = dur.dot(nrm);
            double ta = k_a * (dua - slip_c), Da = k_a;
            // A pile BASE bears, it does not hold: pulled away from the soil under it the tip lifts
            // off and carries nothing (embedded-beam-formulation.md sec 4, Ref. sec 6.6.3.3). The
            // gap is reversible -- the plastic state is kept -- and tension is positive along
            // `tang`, which points from the tip up the beam.
            if (no_tension && ta > 0.0) { ta = 0.0; slip_t = slip_c; Da = 0.0; }
            else if (cap > 0.0 && std::fabs(ta) > cap) { ta = std::copysign(cap, ta); slip_t = dua - ta / k_a; Da = 0.0; }
            else { slip_t = slip_c; }
            if (out_f) *out_f = ta;
            if (out_du) *out_du = dua;
            if (out_capped) *out_capped = (Da == 0.0);
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
                std::array<int, 9> geq, gd;
                for (int k = 0; k < 3; ++k) {
                    Xe(k, 0) = eb.node_x[el[k]]; Xe(k, 1) = eb.node_y[el[k]];
                    gd[3 * k + 0] = ebeam::trans_gdof(eb, el[k], 0, dofs);
                    gd[3 * k + 1] = ebeam::trans_gdof(eb, el[k], 1, dofs);
                    gd[3 * k + 2] = eb.dof_phi[el[k]];
                    for (int c = 0; c < 3; ++c) geq[3 * k + c] = dofs.equation(gd[3 * k + c]);
                }
                const plate::ElementMatrix Kp = plate::stiffness(Xe, eb.props);
                plate::Dof up = plate::Dof::Zero();
                for (int a = 0; a < 9; ++a) up(a) = u_el(gd[a], eb.install);
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
                std::array<int, NC> eq, gdx; std::array<double, NC> cx, cy; int nc = 0;
                for (int i = 0; i < 3; ++i) {
                    for (int c = 0; c < 2; ++c) {
                        gdx[nc] = ebeam::trans_gdof(eb, sp.beam_node[i], c, dofs);
                        eq[nc] = gdx[nc] >= 0 ? dofs.equation(gdx[nc]) : -1;
                        cx[nc] = c == 0 ? sp.Nb(i) : 0.0; cy[nc] = c == 0 ? 0.0 : sp.Nb(i); ++nc;
                    }
                }
                for (int j = 0; j < E::kNodeCount; ++j) {
                    const int sn = mesh.node_of(sp.soil_elem, j);
                    for (int c = 0; c < 2; ++c) {
                        gdx[nc] = dofs.global_dof(sn, c);
                        eq[nc] = dofs.equation(gdx[nc]);
                        cx[nc] = c == 0 ? -Ns(j) : 0.0; cy[nc] = c == 0 ? 0.0 : -Ns(j); ++nc;
                    }
                }
                double f_sp = 0.0, du_sp = 0.0; bool capped = false;
                const bool want = probe && probe->prev_skin;
                axial_couple(eq.data(), gdx.data(), cx.data(), cy.data(), nc, eb.install, sp.tang, sp.k_a, sp.k_n, sp.t_max,
                             sp.wJ, (*st.eskin_c)[skin_off + pi], (*st.eskin_t)[skin_off + pi],
                             want ? &f_sp : nullptr, want ? &du_sp : nullptr,
                             want ? &capped : nullptr);
                // An embedded beam's skin springs are counted with the interface plastic points
                // rather than separately, so they land in the same tally (the local interface
                // error).
                if (want)
                    measure_traction_point(*probe, (*probe->prev_skin)[skin_off + pi], f_sp, du_sp,
                                           sp.k_a, capped, sp.t_max, 0.0);
            }
            skin_off += eb.skin.size();
            if (eb.foot.D_foot > 0.0 && eb.foot.ok) {  // (3) foot (axial spring, cap F_max; wJ=1)
                const auto Ns = E::shape_functions(eb.foot.xi_s, eb.foot.eta_s);
                constexpr int NF = 2 + 2 * E::kNodeCount;
                std::array<int, NF> eq, gdx; std::array<double, NF> cx, cy; int nc = 0;
                for (int c = 0; c < 2; ++c) {
                    gdx[nc] = ebeam::trans_gdof(eb, eb.foot.beam_node, c, dofs);
                    eq[nc] = gdx[nc] >= 0 ? dofs.equation(gdx[nc]) : -1;
                    cx[nc] = c == 0 ? 1.0 : 0.0; cy[nc] = c == 0 ? 0.0 : 1.0; ++nc;
                }
                for (int j = 0; j < E::kNodeCount; ++j) {
                    const int sn = mesh.node_of(eb.foot.soil_elem, j);
                    for (int c = 0; c < 2; ++c) {
                        gdx[nc] = dofs.global_dof(sn, c);
                        eq[nc] = dofs.equation(gdx[nc]);
                        cx[nc] = c == 0 ? -Ns(j) : 0.0; cy[nc] = c == 0 ? 0.0 : -Ns(j); ++nc;
                    }
                }
                double f_ft = 0.0, du_ft = 0.0;
                const bool want_f = probe && probe->prev_foot;
                axial_couple(eq.data(), gdx.data(), cx.data(), cy.data(), nc, eb.install, eb.foot.tang, eb.foot.D_foot, 0.0,
                             eb.foot.f_max, 1.0, (*st.efoot_c)[bi], (*st.efoot_t)[bi],
                             want_f ? &f_ft : nullptr, want_f ? &du_ft : nullptr, nullptr,
                             /*no_tension=*/true);
                // The foot force error. The foot is not counted as a point: its criterion is ONE
                // ratio over every foot in the model, so what is gathered here are the three sums
                // that ratio is formed from. The equilibrium force is the previous iterate's
                // carried forward elastically -- the same equilibrium-versus-constitutive
                // construction as a soil stress point, applied in the same way to the pile tip.
                if (want_f) {
                    PrevPoint& pv = (*probe->prev_foot)[bi];
                    if (!probe->first_iterate) {
                        // A base that has lifted off carries nothing whichever way it is reached.
                        const double f_eq = std::min(0.0, pv.value + eb.foot.D_foot * (du_ft - pv.du));
                        probe->foot_num += std::fabs(f_eq - f_ft);
                    }
                    probe->foot_den_c += std::fabs(f_ft);
                    probe->foot_den_max += std::fabs(eb.foot.f_max > 0.0 ? eb.foot.f_max : 0.0);
                    ++probe->feet;
                    pv.value = f_ft;
                    pv.du = du_ft;
                }
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
    std::vector<LocalErrorPartial> probe_buf_;  // per-element partials, reduced sequentially
};

}  // namespace detail
}  // namespace katai::core
