#include <katai/materials/material_model.hpp>

namespace katai::core {

HsReturnCore hs_return_core(const HardeningSoilParams& pe, double Eur,
                            const Eigen::Vector3d& comm_in_plane, double comm_zz,
                            const Eigen::Vector3d& trial_in_plane, double trial_zz,
                            double gamma_p_n, double pp_n,
                            const HsSubstepPlan* plan_in,
                            double sigma_t_cap,
                            HsSubstepPlan* plan_out,
                            // The constitutive integration error tolerance (STOL), from
                            // the phase (.k2d v15) or 0 for the material class's default.
                            double substep_tol) {
    const double cxx = comm_in_plane(0), cyy = comm_in_plane(1), cxy = comm_in_plane(2);
    const double cmean = 0.5 * (cxx + cyy);
    const double cR = std::sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);

    // Spectral decomposition of the trial stress (tension-positive), tracking sources.
    const double sxx = trial_in_plane(0), syy = trial_in_plane(1), sxy = trial_in_plane(2);
    const double mean = 0.5 * (sxx + syy), hd = 0.5 * (sxx - syy);
    const double radius = std::sqrt(hd * hd + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = hd / radius; sin2t = sxy / radius; }
    struct PV { double v; int src; };
    PV pv[3] = {{mean + radius, 0}, {mean - radius, 1}, {trial_zz, 2}};
    std::sort(pv, pv + 3, [](const PV& a, const PV& b) { return a.v > b.v; });

    // A POINT AT THE TENSION CORNER THAT IS PULLED FURTHER STAYS THERE. With all three principal
    // stresses at the tension cut-off, an elastic trial whose own cut-off return is the corner
    // again -- the hydrostatic apex of the cut-off -- has the corner as its answer whenever the
    // corner is admissible for the model's own surfaces: it is the closest admissible point of
    // the smaller set too. The walk below would get there as well, but from a stress-free state
    // it runs to its ceiling of substeps on the way (the stiffness there is the floor of a power
    // law) and reports an integration that did not meet its tolerance for a point that carries
    // nothing: measured on a strip load on sand, 0.6% of the integrations and 37% of all the
    // substeps taken. No strain hardens anything there, so the state variables stand.
    if (sigma_t_cap < 0.5 * kNoTensionCap) {
        const double tolc = 1e-10 * (1.0 + std::fabs(sigma_t_cap));
        const double cmin = std::min(cmean - cR, comm_zz), cmax = std::max(cmean + cR, comm_zz);
        const double sphi = std::sin(pe.friction);
        const bool in_cone = !(sphi > 1e-12) ||
                             pe.cohesion * std::cos(pe.friction) / sphi >= sigma_t_cap;
        const bool in_cap = !(pe.cap_beta > 0.0) ||
                            sigma_t_cap * sigma_t_cap <= pp_n * pp_n + 1e-20;
        if (cmin >= sigma_t_cap - tolc && cmax <= sigma_t_cap + tolc && in_cone && in_cap) {
            double rs[3] = {pv[0].v, pv[1].v, pv[2].v};
            const LameConstants lc = lame_from(Eur, pe.nu_ur);
            if (apply_rankine_cap(rs, sigma_t_cap, lc.lambda, lc.mu, nullptr) &&
                rs[0] == sigma_t_cap && rs[1] == sigma_t_cap && rs[2] == sigma_t_cap) {
                HsReturnCore out;
                out.in_plane << sigma_t_cap, sigma_t_cap, 0.0;
                out.zz = sigma_t_cap;
                out.gamma_p = gamma_p_n;
                out.pp = pp_n;
                out.tan.tangent.setZero();        // the corner does not depend on the strain
                out.tan.algo_jacobian.setZero();
                out.plastic = true;
                out.nsub = 0;
                if (plan_out) plan_out->clear();
                return out;
            }
        }
    }

    // Robust strain-driven substepping return (hs_integrate). What is handed over is a STRAIN
    // increment: C_e (sigma_tr - sigma_n) inverts the elastic predictor the FE just formed with
    // this same Eur, so it recovers exactly the strain increment the element applied -- that is
    // the physical input, and it is exact whatever the integrator does with it afterwards. The
    // earlier nested-Newton projection diverged on shear-dominated, low-confinement BVPs (strip
    // footings, free surfaces).
    // (Before 2026-08-20 the integrator held E_ur fixed at this same value, so a step that
    // turned out elastic reproduced the predictor bit for bit. It no longer does: E_ur follows
    // sigma3 through the substeps, which makes the unloading-reloading response nonlinear, which
    // is what the model says it is. The Jacobian chain below is unaffected -- deps/dsigma_tr is
    // C_e exactly, by construction of the predictor.)
    const Eigen::Vector3d sigHS(-pv[2].v, -pv[1].v, -pv[0].v);   // trial, comp-pos desc
    // The committed stress IN THE TRIAL'S PRINCIPAL FRAME, position by position: the normal
    // components along the trial's eigen-directions. With isotropic elasticity the difference
    // from the trial is then exactly D_e applied to the strain increment's components in that
    // frame, so C_e below recovers the increment the element applied. Until 2026-09 the committed
    // principals were sorted on their own and subtracted from the sorted trial: wherever the
    // ORDER of the principal values changed between the two -- the axial stress of a triaxial
    // extension passing from major to minor, or two nearly equal stresses swapping -- components
    // of different directions were subtracted from each other.
    const double chd = 0.5 * (cxx - cyy);
    const double c_A = cmean + chd * cos2t + cxy * sin2t;   // along the trial's direction A
    const double c_B = cmean - chd * cos2t - cxy * sin2t;   // along B
    auto comm_along = [&](int src) { return src == 0 ? c_A : (src == 1 ? c_B : comm_zz); };
    const Eigen::Vector3d sig_n_cp(-comm_along(pv[2].src), -comm_along(pv[1].src),
                                   -comm_along(pv[0].src));   // committed, same positions
    const double nu = pe.nu_ur;
    Eigen::Matrix3d Ce;  // principal elastic compliance at the frozen Eur
    Ce << 1.0, -nu, -nu, -nu, 1.0, -nu, -nu, -nu, 1.0;
    Ce /= Eur;
    const Eigen::Vector3d deps_p = Ce * (sigHS - sig_n_cp);
    const HsIntegrated ret = hs_integrate(pe, sig_n_cp, gamma_p_n, pp_n, deps_p, substep_tol,
                                         plan_out, plan_in);
    double r[3] = {-ret.stress(2), -ret.stress(1), -ret.stress(0)};  // tension, by position
    // Tension cut-off (sigma_i <= sigma_t), applied to the principals the model's own return
    // produced. The positions are the trial's directions; the returned VALUES need not be in the
    // trial's order any more, and the cut-off reads its argument as descending -- so it is
    // applied to the sorted values and the result put back in place.
    bool capped = false;
    Eigen::Matrix3d Jcap = Eigen::Matrix3d::Identity();   // d r_capped / d r, by position
    {
        int perm[3] = {0, 1, 2};
        std::sort(perm, perm + 3, [&](int a, int b) { return r[a] > r[b]; });
        double rs[3] = {r[perm[0]], r[perm[1]], r[perm[2]]};
        const LameConstants lc = lame_from(Eur, pe.nu_ur);
        Eigen::Matrix3d Js;   // in the sorted order the cut-off works in
        capped = apply_rankine_cap(rs, sigma_t_cap, lc.lambda, lc.mu, &Js);
        if (capped)
            for (int i = 0; i < 3; ++i) {
                r[perm[i]] = rs[i];
                for (int j = 0; j < 3; ++j) Jcap(perm[i], perm[j]) = Js(i, j);
            }
    }

    // Analytic consistent (continuum) tangent: the shared principal-space assembly (spin +
    // out-of-plane coupling) used by Mohr-Coulomb, fed J(i,j)=J_comp(2-i,2-j) where J_comp =
    // D_pp C_e is the comp-positive stress-to-stress Jacobian (sign flip leaves J invariant,
    // the sort order reverses). Returns both the 3x3 plane-strain D_T and the 4x4 algo_jacobian.
    const Eigen::Matrix3d Jcomp = ret.tangent * Ce;
    Eigen::Matrix3d Jten;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Jten(i, j) = Jcomp(2 - i, 2 - j);
    // The cut-off acts after the model's own return, so its Jacobian comes after it in the chain.
    if (capped) Jten = Jcap * Jten;
    const int src[3] = {pv[0].src, pv[1].src, pv[2].src};
    const LameConstants lame_ur = lame_from(Eur, nu);

    HsReturnCore out;
    out.tan = principal_consistent_tangent(cos2t, sin2t, radius, src, r, Jten, lame_ur);

    double pa = 0.0, pb = 0.0, pz = 0.0;  // coaxial reconstruction on the trial eigenframe
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == 0) pa = r[i];
        else if (pv[i].src == 1) pb = r[i];
        else pz = r[i];
    }
    const double m2 = 0.5 * (pa + pb), r2 = 0.5 * (pa - pb);
    out.in_plane(0) = m2 + r2 * cos2t;
    out.in_plane(1) = m2 - r2 * cos2t;
    out.in_plane(2) = r2 * sin2t;
    out.zz = pz;
    out.gamma_p = ret.gamma_p;
    out.pp = ret.pp;
    out.plastic = ret.plastic || capped;
    out.nsub = ret.nsub;
    out.saturated = ret.saturated;
    return out;
}

// G0 / G_t of this increment: how much stiffer than the stiffness the increment was taken with
// the point answers an increment that turns against its strain history. 1 without small-strain
// stiffness. (Normal strains exx, eyy, ezz and the ENGINEERING shear gxy.)
static double hs_small_strain_upper(const HardeningSoilParams& p, const GaussState& committed,
                                    double exx, double eyy, double ezz, double gxy) {
    if (!(p.G0_ref > 0.0)) return 1.0;
    const double gfac = 2.0 * (1.0 + p.nu_ur) * p.G0_ref / p.Eur_ref;
    if (!(gfac > 1.0)) return 1.0;
    const double em = (exx + eyy + ezz) / 3.0;
    const HsSmallStrainStep ss =
        hs_small_strain_step(p, committed, exx - em, eyy - em, ezz - em, 0.5 * gxy);
    return gfac / ss.pe.elastic_scale;
}

void hs_forward(const MaterialModel& m, const GaussState& committed,
                const Eigen::Vector3d& de, GaussState& trial,
                Eigen::Matrix3d* tangent_out,
                const HsSubstepPlan* plan_in,
                bool* plastic_out, HsSubstepPlan* plan_out,
                double* Eur_out, double substep_tol,
                int* saturated_out, const Eigen::Vector3d* de_small) {
    // The small-strain step: the deviatoric increment (plane strain, eps_zz = 0) against the
    // strain history, and the stiffness and hardening rate this increment is taken with.
    // de_small, where given, is the increment that step is read from instead of de: the
    // perturbed runs of the consistent tangent pass the base increment, so that they difference
    // the integration at a fixed stiffness and the stiffness' own slope is added once, by
    // hs_small_strain_tangent, with a step that suits it.
    const Eigen::Vector3d& des = de_small ? *de_small : de;
    const double em_ss = (des(0) + des(1)) / 3.0;
    const HsSmallStrainStep ss =
        hs_small_strain_step(m.hs, committed, des(0) - em_ss, des(1) - em_ss, -em_ss, 0.5 * des(2));
    HardeningSoilParams pe = ss.pe;
    // Dilatancy cut-off: psi = 0 clamps the mobilised dilatancy sin(psi_m) to [0, 0] inside the
    // return core, which IS the cut-off rule -- the rule enters at the mobilised dilatancy, and
    // nothing else in the HS machinery has to know about void ratios. The flag carries the SAME
    // fact a second way because for HSsmall the two stopped being equivalent: the Li & Dafalias
    // rule reads psi only through phi_cv, so zeroing psi alone would move phi_cv up to phi and
    // switch the Li & Dafalias contraction on everywhere below failure -- a "stop dilating"
    // option that starts producing volume loss. The cut-off sets psi_m itself back to zero, so
    // it is set back to zero.
    if (dilatancy_cut(m, committed)) { pe.dilatancy = 0.0; pe.dilatancy_cut = true; }
    const double Eur = hs_frozen_Eur(pe, committed.stress, committed.stress_zz);
    const LameConstants lame_ur = lame_from(Eur, pe.nu_ur);
    const PlaneStrainStress comm{committed.stress, committed.stress_zz};
    const PlaneStrainStress pred = elastic_predictor(comm, de, lame_ur);

    const HsReturnCore c = hs_return_core(pe, Eur, committed.stress, committed.stress_zz,
                                          pred.in_plane, pred.zz, committed.gamma_p,
                                          committed.pp, plan_in, tension_cap_of(m), plan_out,
                                          substep_tol);
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.gamma_p = c.gamma_p;
    trial.pp = c.pp;
    trial.eps_vol = committed.eps_vol;
    if (tangent_out) *tangent_out = c.tan.tangent;
    if (plastic_out) *plastic_out = c.plastic;
    if (Eur_out) *Eur_out = Eur;
    if (saturated_out) *saturated_out = c.saturated;

    // HSsmall: the strain history as the increment leaves it -- the increment the stress was
    // integrated over, whichever one the stiffness was read from.
    if (de_small) {
        const double em = (de(0) + de(1)) / 3.0;
        const HsSmallStrainStep sh =
            hs_small_strain_step(m.hs, committed, de(0) - em, de(1) - em, -em, 0.5 * de(2));
        trial.gamma_hist = sh.gamma_hist;
        trial.hs_hist = sh.hist;
        trial.hs_g_mono = sh.g_mono;
    } else {
        trial.gamma_hist = ss.gamma_hist;
        trial.hs_hist = ss.hist;
        trial.hs_g_mono = ss.g_mono;
    }
}

// THE STIFFNESS OF THE INCREMENT DEPENDS ON THE INCREMENT. With the strain history a tensor, the
// small-strain multiplier r = E_t / E_ur is the secant of its curve over this strain increment
// and rises towards G0 where the increment turns against the history; the hardening rate h_i
// moves with it. A tangent taken at fixed r and h_i is not the derivative of the stress:
// d sigma = (that tangent) d eps + (d sigma / d r)(d r / d eps) . d eps + the same for h_i.
// Without the second part the iteration works on an operator that is up to G0/G_ur times too
// soft or too stiff at a reversal -- measured on a tunnel lining installed beside a building on
// piles (HSsmall sands), the phase stopped at 8% of its load with no step along the Newton
// direction reducing the out-of-balance force.
// The multiplier of a small increment is a function of its DIRECTION, not of its size (the
// reversal test reads signs), so its slope is of order 1 / |d eps| and a fixed step of 1e-8 is
// as large as the increment it differences: the step is a fraction of the increment, taken both
// ways, and with no increment there is no direction and no slope. d sigma / d r and
// d sigma / d h_i are one more integration each, replaying the base run's subdivision, and only
// where the gradient is not zero.
Eigen::Matrix3d hs_small_strain_tangent(const MaterialModel& m, const GaussState& committed,
                                        const Eigen::Vector3d& de, const GaussState& trial,
                                        bool plastic, const HsSubstepPlan& plan,
                                        double substep_tol, const Eigen::Matrix3d& tangent) {
    Eigen::Matrix3d out = Eigen::Matrix3d::Zero();
    // An increment below 1e-8 has no slope worth having: its stress change is under a
    // thousandth of a kPa, and the slope, of order 1 / |d eps|, multiplies the rounding
    // of the two integrations below into numbers as large as the stiffness itself. (Measured: a
    // strip load on sand, the linear solve returned NaN from such a tangent.)
    constexpr double kSlopeFloor = 1.0e-8;
    const double h = 1.0e-4 * de.norm();
    if (!(m.hs.G0_ref > 0.0) || !(de.norm() > kSlopeFloor)) return out;
    const auto step_at = [&](const Eigen::Vector3d& d) {
        const double em = (d(0) + d(1)) / 3.0;
        return hs_small_strain_step(m.hs, committed, d(0) - em, d(1) - em, -em, 0.5 * d(2));
    };
    const HsSmallStrainStep ss = step_at(de);
    Eigen::RowVector3d grad_r = Eigen::RowVector3d::Zero(), grad_h = Eigen::RowVector3d::Zero();
    for (int j = 0; j < 3; ++j) {
        Eigen::Vector3d da = de, db = de;
        da(j) += h; db(j) -= h;
        const HardeningSoilParams pa = step_at(da).pe, pb = step_at(db).pe;
        grad_r(j) = (pa.elastic_scale - pb.elastic_scale) / (2.0 * h * ss.pe.elastic_scale);
        grad_h(j) = (pa.hardening_scale - pb.hardening_scale) / (2.0 * h * ss.pe.hardening_scale);
    }
    HardeningSoilParams pe = ss.pe;
    if (dilatancy_cut(m, committed)) { pe.dilatancy = 0.0; pe.dilatancy_cut = true; }
    const PlaneStrainStress comm{committed.stress, committed.stress_zz};
    // The two integrations are differenced over 0.1% of the stiffness, not over a millionth of it:
    // the stress is a smooth function of the multiplier (the subdivision is replayed), and what
    // the difference is multiplied by, the slope of the multiplier, is of order 1 / |d eps|. At
    // 1e-6 the difference of an increment of 1e-7 was the rounding of the stress, and its
    // product with that slope put numbers a hundred times the stiffness into the tangent
    // (measured: the first time step of a consolidation phase that unloads the ground went
    // from 0.17 of its scale to 50 on a step of a ten-thousandth of the Newton direction).
    constexpr double kRel = 1.0e-3;
    const auto stress_with = [&](const HardeningSoilParams& pp) {
        const double Eur2 = hs_frozen_Eur(pp, committed.stress, committed.stress_zz);
        const PlaneStrainStress pred2 = elastic_predictor(comm, de, lame_from(Eur2, pp.nu_ur));
        return hs_return_core(pp, Eur2, committed.stress, committed.stress_zz, pred2.in_plane,
                              pred2.zz, committed.gamma_p, committed.pp, &plan, tension_cap_of(m),
                              nullptr, substep_tol).in_plane;
    };
    if (grad_r.squaredNorm() > 0.0) {
        if (!plastic) {
            // An elastic step is the same stress path walked (1 + k) times as far when the
            // stiffness is scaled by (1 + k): d sigma / d ln r is the tangent applied to the
            // increment, and no integration is needed for it.
            out += (tangent * de) * grad_r;
        } else {
            HardeningSoilParams pp = pe;
            pp.elastic_scale *= 1.0 + kRel;
            out += ((stress_with(pp) - trial.stress) / kRel) * grad_r;
        }
    }
    if (plastic && grad_h.squaredNorm() > 0.0) {
        HardeningSoilParams pp = pe;
        pp.hardening_scale *= 1.0 + kRel;
        out += ((stress_with(pp) - trial.stress) / kRel) * grad_h;
    }
    return out;
}

SsReturnCore ss_return_core(const softsoil::Params& P,
                            const Eigen::Vector3d& comm_in_plane, double comm_zz,
                            const Eigen::Vector4d& de4, double pp_n,
                            int nsub_fixed,
                            const softsoilcreep::Params* ssc,
                            double dt_day,
                            double sigma_t_cap) {
    const double cxx = comm_in_plane(0), cyy = comm_in_plane(1), cxy = comm_in_plane(2);

    // SS elastic trial in Voigt (tension-positive): exponentially EXACT mean + linear
    // deviatoric with G at the trial p — the SAME formula pair as ss_step's predictor
    // (guarantees the inversion below matches to round-off).
    const double mean_c = (cxx + cyy + comm_zz) / 3.0;
    const double p_c = -mean_c;                                       // compression-positive
    const double dv = -(de4(0) + de4(1) + de4(3));                    // compression-positive Δε_v
    const double shift = softsoil::cohesion_shift(P);   // the law is in p' + c cot(phi)
    const double p_tr = softsoil::mean_after(p_c, dv, P.kap_star, shift);
    const double K_tr = softsoil::bulk_at(p_tr, P.kap_star, shift);
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const double em = (de4(0) + de4(1) + de4(3)) / 3.0;
    const double sxx = -p_tr + (cxx - mean_c) + 2.0 * G * (de4(0) - em);
    const double syy = -p_tr + (cyy - mean_c) + 2.0 * G * (de4(1) - em);
    const double szz = -p_tr + (comm_zz - mean_c) + 2.0 * G * (de4(3) - em);
    const double sxy = cxy + G * de4(2);

    // Spectral decomposition of the trial (tension-positive) — the hs_return_core machinery verbatim.
    const double mean2 = 0.5 * (sxx + syy), hd = 0.5 * (sxx - syy);
    const double radius = std::sqrt(hd * hd + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = hd / radius; sin2t = sxy / radius; }
    struct PV { double v; int src; };
    PV pv[3] = {{mean2 + radius, 0}, {mean2 - radius, 1}, {szz, 2}};
    std::sort(pv, pv + 3, [](const PV& a, const PV& b) { return a.v > b.v; });

    // The committed stress IN THE TRIAL'S PRINCIPAL FRAME, position by position -- the same
    // construction as the Hardening Soil wrapper. The committed principals used to be sorted on
    // their own and subtracted from the sorted trial, which pairs components of different
    // directions wherever the ORDER of the values changes within the step. On a Soft Soil path
    // that crosses over, the two values are equal where they cross and the error measured 1.3e-13
    // (test_soft_soil g5); the projection makes the increment exact by construction.
    const double cmean = 0.5 * (cxx + cyy), chd = 0.5 * (cxx - cyy);
    const double c_A = cmean + chd * cos2t + cxy * sin2t;
    const double c_B = cmean - chd * cos2t - cxy * sin2t;
    auto comm_along = [&](int src) { return src == 0 ? c_A : (src == 1 ? c_B : comm_zz); };
    const Eigen::Vector3d sig_n_cp(-comm_along(pv[2].src), -comm_along(pv[1].src),
                                   -comm_along(pv[0].src));   // committed, same positions
    const Eigen::Vector3d sig_tr_cp(-pv[2].v, -pv[1].v, -pv[0].v);  // trial, comp-pos descending

    // Principal Δε = the exact inverse of the elastic law between the two states (the
    // principal mean = the tensor mean, so the p's match the ones above to round-off).
    const double p_n_p = sig_n_cp.mean();
    const double p_tr_p = sig_tr_cp.mean();
    const double dv_p = softsoil::strain_between(p_n_p, p_tr_p, P.kap_star, shift);
    const Eigen::Vector3d ddev =
        ((sig_tr_cp - Eigen::Vector3d::Constant(sig_tr_cp.mean())) -
         (sig_n_cp - Eigen::Vector3d::Constant(sig_n_cp.mean()))) / (2.0 * G);
    const Eigen::Vector3d deps_p = ddev + Eigen::Vector3d::Constant(dv_p / 3.0);

    Eigen::Vector3d ret_sig;
    double ret_pp;
    bool ret_plastic;
    int ret_nsub;
    if (ssc) {
        const softsoilcreep::StepResult rc =
            softsoilcreep::ssc_step(*ssc, sig_n_cp, pp_n, deps_p, dt_day, nsub_fixed);
        ret_sig = rc.sig;
        ret_pp = rc.pp;
        ret_plastic = rc.mc_active || rc.devc > 1e-12;
        ret_nsub = rc.nsub;
    } else {
        const softsoil::StepResult ret = softsoil::ss_step(P, sig_n_cp, pp_n, deps_p, nsub_fixed);
        ret_sig = ret.sig;
        ret_pp = ret.pp;
        ret_plastic = ret.cap_active || ret.mc_active;
        ret_nsub = ret.nsub;
    }
    double r[3] = {-ret_sig(2), -ret_sig(1), -ret_sig(0)};   // tension-positive, by position
    // Tension cut-off (sigma_i <= sigma_t); see the Hardening Soil branch. K_tr and G are this
    // step's secant moduli, so the cap is returned on the same elasticity the step was taken with.
    // The returned values need not be in the trial's order any more, and the cut-off reads its
    // argument as descending, so it is applied to the sorted values and put back in place.
    bool ss_capped = false;
    {
        int perm[3] = {0, 1, 2};
        std::sort(perm, perm + 3, [&](int a, int b) { return r[a] > r[b]; });
        double rs[3] = {r[perm[0]], r[perm[1]], r[perm[2]]};
        ss_capped = apply_rankine_cap(rs, sigma_t_cap, K_tr - 2.0 * G / 3.0, G);
        if (ss_capped)
            for (int i = 0; i < 3; ++i) r[perm[i]] = rs[i];
    }

    SsReturnCore out;
    double pa = 0.0, pb = 0.0, pz = 0.0;   // coaxial reconstruction (in the trial eigenframe)
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == 0) pa = r[i];
        else if (pv[i].src == 1) pb = r[i];
        else pz = r[i];
    }
    const double m2 = 0.5 * (pa + pb), r2 = 0.5 * (pa - pb);
    out.in_plane(0) = m2 + r2 * cos2t;
    out.in_plane(1) = m2 - r2 * cos2t;
    out.in_plane(2) = r2 * sin2t;
    out.zz = pz;
    out.pp = ret_pp;
    out.K = K_tr;
    out.G = G;
    out.plastic = ret_plastic || ss_capped;
    out.nsub = ret_nsub;
    return out;
}

void ss_forward(const MaterialModel& m, const GaussState& committed,
                const Eigen::Vector3d& de, GaussState& trial,
                bool* plastic_out, double* K_out,
                double* G_out, int nsub_fixed, int* nsub_out) {
    Eigen::Vector4d de4;
    de4 << de(0), de(1), de(2), 0.0;
    const SsReturnCore c = ss_return_core(m.ssoil, committed.stress, committed.stress_zz, de4,
                                          committed.pp, nsub_fixed, nullptr, 0.0,
                                          tension_cap_of(m));
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.pp = c.pp;
    trial.gamma_p = committed.gamma_p;
    trial.gamma_hist = committed.gamma_hist;
    trial.eps_vol = committed.eps_vol;
    if (plastic_out) *plastic_out = c.plastic;
    if (K_out) *K_out = c.K;
    if (G_out) *G_out = c.G;
    if (nsub_out) *nsub_out = c.nsub;
}

void ssc_forward(const MaterialModel& m, const GaussState& committed,
                 const Eigen::Vector3d& de, double dt_day, GaussState& trial,
                 bool* plastic_out, double* K_out,
                 double* G_out, int nsub_fixed, int* nsub_out) {
    Eigen::Vector4d de4;
    de4 << de(0), de(1), de(2), 0.0;
    const softsoil::Params S = m.ssc.ss();
    const SsReturnCore c = ss_return_core(S, committed.stress, committed.stress_zz, de4,
                                          committed.pp, nsub_fixed, &m.ssc, dt_day,
                                          tension_cap_of(m));
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.pp = c.pp;
    trial.gamma_p = committed.gamma_p;
    trial.gamma_hist = committed.gamma_hist;
    trial.eps_vol = committed.eps_vol;
    if (plastic_out) *plastic_out = c.plastic;
    if (K_out) *K_out = c.K;
    if (G_out) *G_out = c.G;
    if (nsub_out) *nsub_out = c.nsub;
}

void integrate_point(const MaterialModel& m, const GaussState& committed,
                     const Eigen::Vector3d& strain_increment,
                     GaussState& trial, Eigen::Matrix3d& tangent,
                     TangentMode mode, double dt_day,
                     PointReport* report, double substep_tol,
                     const Eigen::Vector3d* small_strain_increment) {
    const LameConstants lame = lame_from(m.youngs_modulus, m.poisson_ratio);
    const PlaneStrainStress previous{committed.stress, committed.stress_zz};
    const PlaneStrainStress predictor =
        elastic_predictor(previous, strain_increment, lame);

    // PLASTICITY EXCLUDED (the elastic zone at the foot of a plate with end bearing): the point
    // answers with the elastic part of its own model -- for Hardening Soil the unloading
    // stiffness at its stress, for the others E and nu -- and its hardening state stays where it
    // is. The compression-index models have no such branch here and are integrated as usual.
    if (committed.plasticity_excluded && m.type != MaterialType::SoftSoil &&
        m.type != MaterialType::SoftSoilCreep) {
        trial = committed;
        if (m.type == MaterialType::HardeningSoil) {
            const HardeningSoilParams pe = hs_small_strain_params(m.hs, committed.gamma_hist);
            const double Eur = hs_frozen_Eur(pe, committed.stress, committed.stress_zz);
            const PlaneStrainStress pred =
                elastic_predictor(previous, strain_increment, lame_from(Eur, pe.nu_ur));
            trial.stress = pred.in_plane;
            trial.stress_zz = pred.zz;
            tangent = elastic_plane_strain_of(Eur, pe.nu_ur);
            if (report) report->stress_dependent = true;
        } else {
            trial.stress = predictor.in_plane;
            trial.stress_zz = predictor.zz;
            tangent = m.elastic_plane_strain();
        }
        if (report) report->elastic = tangent;
        return;
    }

    switch (m.type) {
        case MaterialType::LinearElastic: {
            trial.stress = predictor.in_plane;
            trial.stress_zz = predictor.zz;
            tangent = m.elastic_plane_strain();
            if (report) report->elastic = tangent;
            break;
        }
        case MaterialType::MohrCoulomb: {
            const MohrCoulombParams params{m.youngs_modulus, m.poisson_ratio,
                                           m.cohesion, m.friction_angle,
                                           effective_dilatancy(m, committed), m.tension_cutoff,
                                           m.tensile_strength};
            const McReturn base = mc_return_mapping(predictor, params);
            trial.stress = base.stress.in_plane;
            trial.stress_zz = base.stress.zz;
            // Elastic step -> elastic operator; plastic step -> the closed-form
            // consistent (algorithmic) tangent assembled inside the return mapping
            // (the active region's Jacobian carried through the eigen-decomposition,
            // incl. in-plane frame rotation and the sigma_zz coupling). A region-
            // consistent analytic tangent is essential for phi > 0: a forward
            // finite difference straddles region boundaries near the edges
            // (sigma_2 = sigma_3, ubiquitous under gravity) and breaks global
            // Newton convergence.
            tangent = base.plastic ? base.tangent : m.elastic_plane_strain();
            if (report) {
                report->plastic = base.plastic;
                report->elastic = m.elastic_plane_strain();
            }
            break;
        }
        case MaterialType::HoekBrown: {
            // Rock: Hooke plus a non-linear strength criterion, so the trial IS the elastic
            // predictor the caller already built and the only work is the return. The tangent is
            // the elastic operator (hoek_brown.hpp says why); the equilibrium iteration pays for
            // that in steps rather than in accuracy, which is the trade the Hardening Soil branch
            // makes for its own reasons.
            hoekbrown::Params hp = m.hb;
            hp.E = m.youngs_modulus; hp.nu = m.poisson_ratio;
            const hoekbrown::Constants hc = hoekbrown::constants_of(hp);
            const hoekbrown::PlaneReturn r = hoekbrown::plane_return(predictor, hp, hc);
            trial.stress = r.stress.in_plane;
            trial.stress_zz = r.stress.zz;
            tangent = m.elastic_plane_strain();
            if (report) {
                report->plastic = r.plastic;
                report->elastic = m.elastic_plane_strain();
            }
            // THE CONSISTENT TANGENT, BY FINITE DIFFERENCE -- the same device the Hardening
            // Soil and Soft Soil branches use, for the same reason and at a smaller cost (this
            // return is a bisection on one scalar, not a substepped walk). The elastic operator
            // alone was enough to integrate the material correctly, and every material-point
            // test passes with it, because those tests ask what stress comes back. It is NOT
            // enough to iterate an ill-conditioned boundary value problem to equilibrium:
            // measured on a circular tunnel unloaded into a Hoek-Brown rock mass, the elastic
            // tangent stalled as soon as the plastic annulus passed about 4% of the tunnel
            // radius, and no number of load steps bought more than a little (150 steps per
            // stage reached the same wall as 40). Newton with a wrong-but-symmetric operator
            // converges linearly, and linearly is not always convergent.
            if (mode == TangentMode::kConsistent && r.plastic) {
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j));
                    dep(j) += h;
                    const PlaneStrainStress pp = elastic_predictor(previous, dep, lame);
                    const hoekbrown::PlaneReturn rp = hoekbrown::plane_return(pp, hp, hc);
                    tangent(0, j) = (rp.stress.in_plane(0) - r.stress.in_plane(0)) / h;
                    tangent(1, j) = (rp.stress.in_plane(1) - r.stress.in_plane(1)) / h;
                    tangent(2, j) = (rp.stress.in_plane(2) - r.stress.in_plane(2)) / h;
                }
            }
            break;
        }
        case MaterialType::HardeningSoil: {
            bool plastic = false;
            double Eur = 0.0;
            int saturated = 0;
            HsSubstepPlan plan;   // the base run's subdivision, replayed by the perturbed runs
            hs_forward(m, committed, strain_increment, trial, &tangent, nullptr, &plastic, &plan,
                       &Eur, substep_tol, &saturated, small_strain_increment);
            if (report) {
                report->plastic = plastic;
                report->stress_dependent = true;
                report->integration_saturated = saturated != 0;
                report->elastic = elastic_plane_strain_of(Eur, m.hs.nu_ur);
                const double up = hs_small_strain_upper(m.hs, committed, strain_increment(0),
                                                        strain_increment(1), 0.0,
                                                        strain_increment(2));
                if (up > 1.0) {
                    report->has_stiffest = true;
                    report->stiffest = elastic_plane_strain_of(up * Eur, m.hs.nu_ur);
                }
            }
            // kConsistent + plastic step → numerical consistent tangent (TangentMode
            // block): 3 perturbed forward runs, substep count pinned to the base run.
            if (mode == TangentMode::kConsistent && plastic) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j), strain_increment.norm());
                    dep(j) += h;
                    hs_forward(m, committed, dep, pert, nullptr, &plan, nullptr, nullptr,
                               nullptr, substep_tol, nullptr,
                               small_strain_increment ? small_strain_increment
                                                      : &strain_increment);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            // (No slope of the stiffness where it is read from a fixed increment.)
            if (mode != TangentMode::kNone && !small_strain_increment)
                tangent += hs_small_strain_tangent(m, committed, strain_increment, trial, plastic,
                                                   plan, substep_tol, tangent);
            break;
        }
        case MaterialType::SoftSoil: {
            bool plastic = false;
            double K = 0.0, G = 0.0;
            int nsub = 0;
            ss_forward(m, committed, strain_increment, trial, &plastic, &K, &G, 0, &nsub);
            tangent = elastic_plane_strain_kg(K, G);
            if (report) {
                report->plastic = plastic;
                report->stress_dependent = true;
                report->elastic = tangent;
            }
            // Plastic step + tangent wanted → NUMERICAL (forward-difference) tangent —
            // unlike HS, in BOTH modes (kContinuum/kConsistent): SS has no "cheap
            // continuum", the elastic operator is λ*/κ* times too stiff on the cap and
            // breaks Newton. kNone: the tangent is unused, skip (the elastic operator
            // stays, harmless). The perturbed runs are pinned to the base run's SUBSTEP
            // count (nsub — the HS lesson).
            if (plastic && mode != TangentMode::kNone) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    ss_forward(m, committed, dep, pert, nullptr, nullptr, nullptr, nsub);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            break;
        }
        case MaterialType::SoftSoilCreep: {
            // Same skeleton as the SS branch; time enters via dt_day and the FD columns run
            // with the SAME dt + SAME nsub (the tangent ∂σ/∂ε is taken at fixed time —
            // which is exactly what Newton needs).
            bool plastic = false;
            double K = 0.0, G = 0.0;
            int nsub = 0;
            ssc_forward(m, committed, strain_increment, dt_day, trial, &plastic, &K, &G, 0, &nsub);
            tangent = elastic_plane_strain_kg(K, G);
            if (report) {
                report->plastic = plastic;
                report->stress_dependent = true;
                report->elastic = tangent;
            }
            if (plastic && mode != TangentMode::kNone) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    ssc_forward(m, committed, dep, dt_day, pert, nullptr, nullptr, nullptr, nsub);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            break;
        }
    }
}

void integrate_point_axisym(const MaterialModel& m,
                            const GaussState& committed,
                            const Eigen::Vector4d& strain_increment,
                            GaussState& trial, Eigen::Matrix4d& tangent,
                            TangentMode mode,
                            double dt_day,
                            PointReport4* report,
                            double substep_tol,
                            const Eigen::Vector4d* small_strain_increment) {
    const Eigen::Matrix4d De = m.elastic_axisym();
    Eigen::Vector4d s_n;
    s_n << committed.stress, committed.stress_zz;       // [r, z, rz, theta]
    const Eigen::Vector4d s_tr = s_n + De * strain_increment;  // elastic predictor

    switch (m.type) {
        case MaterialType::LinearElastic: {
            trial.stress = s_tr.head<3>();
            trial.stress_zz = s_tr(3);
            tangent = De;
            if (report) report->elastic = De;
            break;
        }
        case MaterialType::MohrCoulomb: {
            PlaneStrainStress predictor;
            predictor.in_plane = s_tr.head<3>();
            predictor.zz = s_tr(3);
            const MohrCoulombParams params{m.youngs_modulus, m.poisson_ratio,
                                           m.cohesion, m.friction_angle,
                                           effective_dilatancy(m, committed), m.tension_cutoff,
                                           m.tensile_strength};
            const McReturn base = mc_return_mapping(predictor, params);
            trial.stress = base.stress.in_plane;
            trial.stress_zz = base.stress.zz;
            tangent = base.plastic ? (base.algo_jacobian * De) : De;
            if (report) { report->plastic = base.plastic; report->elastic = De; }
            break;
        }
        case MaterialType::HoekBrown: {
            // Axisymmetric Hoek-Brown. The (r,z) block is the in-plane Mohr circle and the hoop
            // sigma_theta is the third principal -- shear on the theta planes vanishes in
            // axisymmetry, exactly as sigma_zz is the third principal in plane strain -- so
            // plane_return is reused verbatim and only the elastic predictor differs: it uses
            // the AXISYMMETRIC operator, because the hoop strain is a real strain here and not
            // a constraint. This is the same argument the Hardening Soil block below makes.
            //
            // THIS CASE WAS MISSING when the model shipped, and the switch has no default, so a
            // Hoek-Brown material in axisymmetry left `trial` and `tangent` exactly as the
            // caller passed them -- the previous iterate, returned as if it had been integrated.
            // Nothing failed loudly: the equilibrium iteration simply never converged, and the
            // run reported a stalled search. A missing branch in a switch over a closed enum is
            // the cheapest silent-wrong there is, which is why test_material_axisym_coverage now
            // walks the registry and poisons the outputs first.
            hoekbrown::Params hp = m.hb;
            hp.E = m.youngs_modulus; hp.nu = m.poisson_ratio;
            const hoekbrown::Constants hc = hoekbrown::constants_of(hp);
            PlaneStrainStress predictor;
            predictor.in_plane = s_tr.head<3>();
            predictor.zz = s_tr(3);
            const hoekbrown::PlaneReturn r = hoekbrown::plane_return(predictor, hp, hc);
            trial.stress = r.stress.in_plane;
            trial.stress_zz = r.stress.zz;
            tangent = De;
            if (report) { report->plastic = r.plastic; report->elastic = De; }
            // The consistent 4x4 by finite difference, as in the plane-strain branch; the
            // fourth column and row are the HOOP, which is a real strain here.
            if (mode == TangentMode::kConsistent && r.plastic) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j));
                    dep(j) += h;
                    const Eigen::Vector4d s_tr_p = s_n + De * dep;
                    PlaneStrainStress pp;
                    pp.in_plane = s_tr_p.head<3>();
                    pp.zz = s_tr_p(3);
                    const hoekbrown::PlaneReturn rp = hoekbrown::plane_return(pp, hp, hc);
                    tangent(0, j) = (rp.stress.in_plane(0) - r.stress.in_plane(0)) / h;
                    tangent(1, j) = (rp.stress.in_plane(1) - r.stress.in_plane(1)) / h;
                    tangent(2, j) = (rp.stress.in_plane(2) - r.stress.in_plane(2)) / h;
                    tangent(3, j) = (rp.stress.zz - r.stress.zz) / h;
                }
            }
            break;
        }
        case MaterialType::HardeningSoil: {
            // Axisymmetric HS: the (r,z) block is the in-plane Mohr circle and the hoop
            // sigma_theta is the third principal -- the SAME principal structure as plane
            // strain, so hs_return_core is reused verbatim. Only the elastic predictor /
            // operator differ: the trial uses the AXISYMMETRIC elastic matrix at the frozen
            // Eur (NOT m.youngs_modulus; the hoop is a real strain), and the consistent 4x4
            // tangent is D_T = algo_jacobian * D_e_axisym (Psi from the shared principal
            // assembly), exactly mirroring Mohr-Coulomb.
            // The small-strain step, with the real hoop strain in the deviator.
            const auto small_step = [&](const Eigen::Vector4d& d) {
                const double em = (d(0) + d(1) + d(3)) / 3.0;
                return hs_small_strain_step(m.hs, committed, d(0) - em, d(1) - em, d(3) - em,
                                            0.5 * d(2));
            };
            // The history is advanced by the increment; the stiffness is read from it too unless
            // the caller fixed another one (see small_strain_increment in the header).
            const HsSmallStrainStep ss = small_step(strain_increment);
            HardeningSoilParams pe =
                small_strain_increment ? small_step(*small_strain_increment).pe : ss.pe;
            // Dilatancy cut-off: see the plane-strain block (hs_forward).
            if (dilatancy_cut(m, committed)) { pe.dilatancy = 0.0; pe.dilatancy_cut = true; }
            const double Eur = hs_frozen_Eur(pe, committed.stress, committed.stress_zz);
            const double nu = pe.nu_ur;
            const Eigen::Matrix4d De_ur = elastic_axisym_of(Eur, nu);
            const Eigen::Vector4d s_tr_ur = s_n + De_ur * strain_increment;
            HsSubstepPlan plan;   // as in the plane-strain block: the base run's subdivision
            const HsReturnCore c = hs_return_core(pe, Eur, committed.stress, committed.stress_zz,
                                                  s_tr_ur.head<3>(), s_tr_ur(3),
                                                  committed.gamma_p, committed.pp, nullptr,
                                                  tension_cap_of(m), &plan, substep_tol);
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.gamma_p = c.gamma_p;
            trial.pp = c.pp;
            trial.eps_vol = committed.eps_vol;
            tangent = c.tan.algo_jacobian * De_ur;  // continuum 4x4 (Psi * D_e)
            if (report) {
                report->plastic = c.plastic;
                report->stress_dependent = true;
                report->integration_saturated = c.saturated != 0;
                report->elastic = De_ur;
                const double up = hs_small_strain_upper(m.hs, committed, strain_increment(0),
                                                        strain_increment(1), strain_increment(3),
                                                        strain_increment(2));
                if (up > 1.0) {
                    report->has_stiffest = true;
                    report->stiffest = elastic_axisym_of(up * Eur, m.hs.nu_ur);
                }
            }
            // kConsistent + plastic step → numerical consistent 4x4 tangent (same rationale
            // as the plane-strain block; the perturbed runs replay the base subdivision).
            if (mode == TangentMode::kConsistent && c.plastic) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j), strain_increment.norm());
                    dep(j) += h;
                    const Eigen::Vector4d s_tr_p = s_n + De_ur * dep;
                    const HsReturnCore cp = hs_return_core(
                        pe, Eur, committed.stress, committed.stress_zz,
                        s_tr_p.head<3>(), s_tr_p(3), committed.gamma_p, committed.pp,
                        &plan, tension_cap_of(m), nullptr, substep_tol);
                    tangent(0, j) = (cp.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cp.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cp.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cp.zz - c.zz) / h;
                }
            }
            // The stiffness of the increment depends on the increment: the same two terms as in
            // the plane-strain routine (hs_small_strain_tangent), with the hoop strain as a fourth component. Added
            // after the consistent tangent above, whose perturbed runs keep the stiffness fixed.
            if (mode != TangentMode::kNone && m.hs.G0_ref > 0.0 && !small_strain_increment) {
                Eigen::RowVector4d grad_r = Eigen::RowVector4d::Zero(),
                                   grad_h = Eigen::RowVector4d::Zero();
                // see hs_small_strain_tangent, also for the floor
                const double h = 1.0e-4 * strain_increment.norm();
                for (int j = 0; j < 4 && strain_increment.norm() > 1.0e-8; ++j) {
                    const auto at = [&](double sgn) {
                        Eigen::Vector4d dp = strain_increment;
                        dp(j) += sgn * h;
                        const double emp = (dp(0) + dp(1) + dp(3)) / 3.0;
                        return hs_small_strain_step(m.hs, committed, dp(0) - emp, dp(1) - emp,
                                                    dp(3) - emp, 0.5 * dp(2)).pe;
                    };
                    const HardeningSoilParams pa = at(1.0), pb = at(-1.0);
                    grad_r(j) = (pa.elastic_scale - pb.elastic_scale) / (2.0 * h * ss.pe.elastic_scale);
                    grad_h(j) = (pa.hardening_scale - pb.hardening_scale) /
                                (2.0 * h * ss.pe.hardening_scale);
                }
                constexpr double kRel = 1.0e-3;   // see hs_small_strain_tangent
                const auto stress_with = [&](const HardeningSoilParams& pp) {
                    const double Eur2 = hs_frozen_Eur(pp, committed.stress, committed.stress_zz);
                    const Eigen::Vector4d s_tr2 =
                        s_n + elastic_axisym_of(Eur2, pp.nu_ur) * strain_increment;
                    const HsReturnCore c2 = hs_return_core(
                        pp, Eur2, committed.stress, committed.stress_zz, s_tr2.head<3>(), s_tr2(3),
                        committed.gamma_p, committed.pp, &plan, tension_cap_of(m), nullptr,
                        substep_tol);
                    Eigen::Vector4d s;
                    s << c2.in_plane, c2.zz;
                    return s;
                };
                Eigen::Vector4d s_base;
                s_base << c.in_plane, c.zz;
                if (grad_r.squaredNorm() > 0.0) {
                    if (!c.plastic) {   // see hs_small_strain_tangent
                        const Eigen::Vector4d along = tangent * strain_increment;
                        tangent += along * grad_r;
                    } else {
                        HardeningSoilParams pp = pe;
                        pp.elastic_scale *= 1.0 + kRel;
                        tangent += ((stress_with(pp) - s_base) / kRel) * grad_r;
                    }
                }
                if (c.plastic && grad_h.squaredNorm() > 0.0) {
                    HardeningSoilParams pp = pe;
                    pp.hardening_scale *= 1.0 + kRel;
                    tangent += ((stress_with(pp) - s_base) / kRel) * grad_h;
                }
            }
            // HSsmall: the strain history as the increment leaves it.
            trial.gamma_hist = ss.gamma_hist;
            trial.hs_hist = ss.hist;
            trial.hs_g_mono = ss.g_mono;
            break;
        }
        case MaterialType::SoftSoil: {
            // Axisymmetric SS: the (r,z) block is the in-plane Mohr circle, the hoop σ_θ is
            // the third principal — ss_return_core takes de4 DIRECTLY as [Δεr, Δεz, Δγrz,
            // Δεθ] (the hoop is a real strain; the volumetric term contains all three).
            // Tangent: elastic (K_tr,G) isotropic 4×4; plastic by FD, same rationale as the
            // plane-strain branch.
            const SsReturnCore c = ss_return_core(m.ssoil, committed.stress, committed.stress_zz,
                                                  strain_increment, committed.pp, 0, nullptr,
                                                  0.0, tension_cap_of(m));
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.pp = c.pp;
            trial.gamma_p = committed.gamma_p;
            trial.gamma_hist = committed.gamma_hist;
            trial.eps_vol = committed.eps_vol;
            tangent = elastic_axisym_kg(c.K, c.G);
            if (report) {
                report->plastic = c.plastic;
                report->stress_dependent = true;
                report->elastic = tangent;
            }
            if (c.plastic && mode != TangentMode::kNone) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    const SsReturnCore cq =
                        ss_return_core(m.ssoil, committed.stress, committed.stress_zz, dep,
                                       committed.pp, c.nsub, nullptr, 0.0, tension_cap_of(m));
                    tangent(0, j) = (cq.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cq.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cq.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cq.zz - c.zz) / h;
                }
            }
            break;
        }
        case MaterialType::SoftSoilCreep: {
            // Same skeleton as the SS axisymmetric branch + time (dt_day); FD columns with
            // the same dt + the base nsub.
            const softsoil::Params S = m.ssc.ss();
            const SsReturnCore c = ss_return_core(S, committed.stress, committed.stress_zz,
                                                  strain_increment, committed.pp, 0, &m.ssc,
                                                  dt_day, tension_cap_of(m));
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.pp = c.pp;
            trial.gamma_p = committed.gamma_p;
            trial.gamma_hist = committed.gamma_hist;
            trial.eps_vol = committed.eps_vol;
            tangent = elastic_axisym_kg(c.K, c.G);
            if (report) {
                report->plastic = c.plastic;
                report->stress_dependent = true;
                report->elastic = tangent;
            }
            if (c.plastic && mode != TangentMode::kNone) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    const SsReturnCore cq =
                        ss_return_core(S, committed.stress, committed.stress_zz, dep,
                                       committed.pp, c.nsub, &m.ssc, dt_day,
                                       tension_cap_of(m));
                    tangent(0, j) = (cq.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cq.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cq.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cq.zz - c.zz) / h;
                }
            }
            break;
        }
    }
}

}  // namespace katai::core
