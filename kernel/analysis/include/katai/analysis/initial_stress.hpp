#pragma once
// Initial stress — the K0 procedure (P2.2). Produces the geostatic (horizontal
// surface/layers) initial stress field DIRECTLY (without solving equations): the vertical
// effective stress from the overburden, the lateral stress σ'_h = K0·σ'_v. This is handed
// to the nonlinear solver as the committed INITIAL Gauss state (pre-stress); in the first
// phase Δε=0 ⇒ zero displacement (the pre-stressed body is already in equilibrium with its
// self-weight). The hook for staged construction (phase-to-phase σ carry-over).
//
// Math + equilibrium proof + sources: docs/references/initial-stress-k0.md
// (Jaky 1944; PLAXIS Reference/Scientific Manual — K0 procedure / Gravity loading).

#include <algorithm>
#include <functional>
#include <vector>

#include <Eigen/Core>

#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/soft_soil.hpp>   // ss_initial_pp (preconsolidation seeding)
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Homogeneous (single-layer) geostatic initial stress. Below the water table unit_weight
// must be the buoyant γ' = γ_sat − γ_w (effective stress; see
// effective-stress-formulation.md); in dry soil the total γ. k0 = σ'_h/σ'_v (Jaky 1−sinφ'
// or a user value).
struct K0Options {
    double surface_elevation = 0.0;  // z_surf (topmost elevation)
    double unit_weight = 0.0;        // γ' (effective/buoyant for the overburden)
    double k0 = 0.0;                 // coefficient of lateral earth pressure
};

namespace detail {
template <class E>
std::vector<GaussState> k0_initial_stress_impl(const mesh::Mesh& mesh,
                                               const K0Options& opt) {
    std::vector<GaussState> states(static_cast<size_t>(mesh.element_count) *
                                   E::kGaussCount);
    const auto gauss = E::gauss_points();
    typename E::NodeCoords coords;
    for (int e = 0; e < mesh.element_count; ++e) {
        for (int k = 0; k < E::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
        }
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto N = E::shape_functions(gauss[g].xi, gauss[g].eta);
            double y = 0.0;
            for (int k = 0; k < E::kNodeCount; ++k) y += N(k) * coords(k, 1);
            // tension-positive: compression negative. σ'_v more negative with depth.
            const double sv = -opt.unit_weight * (opt.surface_elevation - y);
            const double sh = opt.k0 * sv;
            GaussState& s = states[static_cast<size_t>(e) * E::kGaussCount + g];
            s.stress << sh, sv, 0.0;  // [σ_xx=σ_h, σ_yy=σ_v, σ_xy=0]
            s.stress_zz = sh;          // the out-of-plane horizontal direction is K0·σ_v too
            s.eps_vol = 0.0;
        }
    }
    return states;
}
}  // namespace detail

// Returns the K0 initial stress field as Gauss-point states (solver seed). The element
// type is picked from mesh.nodes_per_element (tri6/tri15); the Gauss order matches the
// solver (e*kGaussCount + g).
inline std::vector<GaussState> compute_k0_initial_stress(const mesh::Mesh& mesh,
                                                         const K0Options& opt) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::k0_initial_stress_impl<Tri15Element>(mesh, opt);
    return detail::k0_initial_stress_impl<Tri6Element>(mesh, opt);
}

// --- Generalized (layered + water-aware) K0 initial stress ---------------------------------------
// The homogeneous K0Options assumes a horizontal single layer; a real GUI model can be
// layered (material per polygon), have a water table and an irregular surface. Here σ'_v is
// found at each Gauss point by the VERTICAL integral of the effective unit weight ABOVE
// that point (the classic K0-procedure assumption: column-based overburden — exact for
// horizontally layered soil, approximate on a sloped surface; PLAXIS gives the same
// warning). γ_eff must return the buoyant γ'=γ_sat−γ_w below the water table (effective
// stress). σ'_h = K0·σ'_v, K0 by material id (Jaky 1−sinφ' or a user value).
struct K0LayeredOptions {
    std::function<double(double x, double y)> eff_unit_weight;  // effective (buoyant) γ' at a point
    std::function<double(double x)> ground_surface;            // ground surface elevation y_surf(x)
    std::vector<double> k0;                                     // K0 by material id
    int integration_steps = 400;                               // vertical-ray midpoint steps
    // Discontinuity elevations of γ'(x,·) (layer boundaries, water table) in an x column.
    // If given, the column integral is taken segment by segment: for piecewise-constant γ'
    // the midpoint rule is exact WITHIN a segment, so σ'_v is EXACT (the horizontal-layer K0
    // identity f_int(σ_K0)=f_gravity holds to round-off). If absent, the old single-piece
    // midpoint behaviour.
    std::function<std::vector<double>(double x)> strata_breaks;
    // NonPorous target correction (PLAXIS: a non-porous material sees NEITHER initial NOR
    // excess pore pressure): the eff_unit_weight slice function yields the effective σ'_v
    // for POROUS targets; if the target Gauss point's material is non-porous, the seed must
    // be the TOTAL stress — σ_v = σ'_v − u(x,y) (tension-positive; u ≥ 0 is the pressure
    // magnitude). σ_h = K0·σ_v on the same total. Both empty = old behaviour bit-for-bit.
    std::vector<char> nonporous;                                // by material id (1 = non-porous)
    std::function<double(double x, double y)> pore;             // u(x,y) ≥ 0 hydrostatic magnitude
};

namespace detail {
template <class E>
std::vector<GaussState> k0_layered_impl(const mesh::Mesh& mesh, const K0LayeredOptions& opt) {
    std::vector<GaussState> states(static_cast<size_t>(mesh.element_count) * E::kGaussCount);
    const auto gauss = E::gauss_points();
    typename E::NodeCoords coords;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int mat = mesh.element_material[e];
        const double k0 = (mat >= 0 && mat < static_cast<int>(opt.k0.size())) ? opt.k0[mat] : 0.0;
        for (int k = 0; k < E::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n]; coords(k, 1) = mesh.y[n];
        }
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto N = E::shape_functions(gauss[g].xi, gauss[g].eta);
            double x = 0.0, y = 0.0;
            for (int k = 0; k < E::kNodeCount; ++k) { x += N(k) * coords(k, 0); y += N(k) * coords(k, 1); }
            const double y_surf = opt.ground_surface(x);
            // σ'_v = −∫_y^{y_surf} γ'(x, t) dt (midpoint rule; compression negative). If the
            // discontinuity elevations (strata_breaks) are known, integrate segment by
            // segment — exact for piecewise-constant γ' (removes the O(Δγ·dt) error of a
            // step straddling a layer boundary).
            double sv = 0.0;
            if (y_surf > y) {
                std::vector<double> seg{y, y_surf};
                if (opt.strata_breaks) {
                    for (double b : opt.strata_breaks(x))
                        if (b > y + 1e-12 && b < y_surf - 1e-12) seg.push_back(b);
                    std::sort(seg.begin(), seg.end());
                }
                const int ns_total = std::max(1, opt.integration_steps);
                for (size_t sgi = 0; sgi + 1 < seg.size(); ++sgi) {
                    const double a = seg[sgi], b = seg[sgi + 1];
                    const int ns = std::max(1, (int)((b - a) / (y_surf - y) * ns_total + 0.5));
                    const double dt = (b - a) / ns;
                    for (int i = 0; i < ns; ++i) sv += opt.eff_unit_weight(x, a + (i + 0.5) * dt) * dt;
                }
            }
            double sigv = -sv;
            // Non-porous target: the seed is the total stress (header note above).
            if (mat >= 0 && mat < static_cast<int>(opt.nonporous.size()) &&
                opt.nonporous[mat] && opt.pore)
                sigv -= opt.pore(x, y);
            const double sigh = k0 * sigv;
            GaussState& s = states[static_cast<size_t>(e) * E::kGaussCount + g];
            s.stress << sigh, sigv, 0.0;
            s.stress_zz = sigh;
            s.eps_vol = 0.0;
        }
    }
    return states;
}
}  // namespace detail

inline std::vector<GaussState> compute_k0_initial_stress_layered(const mesh::Mesh& mesh,
                                                                 const K0LayeredOptions& opt) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::k0_layered_impl<Tri15Element>(mesh, opt);
    return detail::k0_layered_impl<Tri6Element>(mesh, opt);
}

// Overconsolidation raises the automatic K0 (PLAXIS Reference, elastic unloading):
//   K0 = K0nc OCR - nu/(1 - nu) (OCR - 1),
// clamped to [0, Kp = (1 + sin phi)/(1 - sin phi)] as a coarse passive-limit
// safeguard. nu is nu_ur for the advanced models and nu for LE/MC (the caller
// chooses); POP's contribution to K0 is depth-dependent (a sigma'_v ratio) and
// is deliberately NOT folded in here in v1 -- it enters the cap seed only,
// stated in the GUI help. Before this correction OCR reached only the HS/SS
// cap seed and was SILENTLY ignored by the MC/LE geostatic field (audit
// finding).
inline double k0_overconsolidated(double k0nc, double ocr, double nu, double sin_phi) {
    const double nu_c = std::clamp(nu, 0.0, 0.49);
    const double k0 = k0nc * ocr - nu_c / (1.0 - nu_c) * (ocr - 1.0);
    return std::clamp(k0, 0.0, (1.0 + sin_phi) / std::max(1e-6, 1.0 - sin_phi));
}

// Per-material overconsolidation input for the preconsolidation seed, resolved
// from the schema at the caller's seam. mode: 0 = normally consolidated,
// 1 = OCR given, 2 = POP given.
struct Overconsolidation {
    int mode = 0;
    double OCR = 1.0;
    double POP = 0.0;
};

// Hardening Soil / Soft Soil (Creep): seed the cap preconsolidation pp and the
// shear hardening gamma_p so the geostatic K0 state is admissible (else the
// first stress return is a large inconsistent correction and the solver
// diverges -- hardening-soil-formulation.md sec 7). `states` is the K0 seed in
// solver order (element * gauss); inactive elements are skipped (their seed is
// zero by construction). Model-family selection reads the engine's models
// table, never a schema enum.
inline void seed_preconsolidation(const mesh::Mesh& mesh,
                                  const std::vector<MaterialModel>& models,
                                  const std::vector<Overconsolidation>& oc,
                                  const std::vector<char>& active,
                                  std::vector<GaussState>& states) {
    const int ng = mesh.element_count > 0 ? (int)(states.size() / mesh.element_count) : 0;
    for (int e = 0; e < mesh.element_count && ng > 0; ++e) {
        if (!active.empty() && !active[e]) continue;
        const int mat = mesh.element_material[e];
        const bool is_hs = models[mat].type == MaterialType::HardeningSoil;
        const bool is_ss = models[mat].type == MaterialType::SoftSoil ||
                           models[mat].type == MaterialType::SoftSoilCreep;
        if (!is_hs && !is_ss) continue;
        for (int g = 0; g < ng; ++g) {
            auto& gs = states[(size_t)e * ng + g];
            const Eigen::Vector3d sig_cp(-gs.stress(0), -gs.stress(1), -gs.stress_zz);
            // OCR_eq: OCR mode directly; POP mode as the equivalent ratio
            // (sigma'_v0 + POP)/sigma'_v0 (vertical = yy). POP used to be
            // SILENTLY ignored (only OCR was wired) -- hs_initial_pp is
            // first-order homogeneous in sigma, so the ratio-scaling is the
            // exact counterpart of PLAXIS's POP definition; for SS the same
            // ratio applies to f-bar (see the ss_initial_pp note).
            double ocr_eq = 1.0;
            if (oc[mat].mode == 1) {
                ocr_eq = std::max(1.0, oc[mat].OCR);
            } else if (oc[mat].mode == 2) {
                const double sv0 = sig_cp(1);
                if (sv0 > 1e-9)
                    ocr_eq = std::max(1.0, (sv0 + oc[mat].POP) / sv0);
            }
            if (is_ss) {
                // SSC: the elastic law and the p_eq measure are IDENTICAL to SS
                // (the ssc.ss() view), so the same seed is correct. Creep itself
                // does not change the seed (the t = 0 instant; ageing is extra).
                const softsoil::Params sp =
                    models[mat].type == MaterialType::SoftSoil
                        ? models[mat].ssoil
                        : models[mat].ssc.ss();
                gs.pp = softsoil::ss_initial_pp(sp, sig_cp, ocr_eq);
                continue;
            }
            const double ocr = ocr_eq;
            // Seed with the parameter set the FIRST integrate_point will actually use. The whole
            // job of hs_initial_gamma_p is to put the pre-stress exactly ON the shear surface
            // (gamma_p = fbar(q0)), and fbar depends on E_ur. But hs_forward evaluates the model
            // through hs_small_strain_params(hs, gamma_hist), which for HSsmall (G0_ref > 0)
            // REPLACES E_ur,ref by E_0,ref = 2(1+nu_ur) G_0,ref -- 288 vs 90 MPa on the defaults.
            // Seeding with the raw params therefore put the state OFF the surface the model then
            // evaluated: measured f_s = +1.45e-3 and |dsigma| = 9.9 kPa (5.5% of sigma_v) for a ZERO
            // strain increment, so the K0 phase started with a residual it could never clear and
            // reported "equilibrated 0%". gamma_hist = 0 is the seed's own state (no strain history
            // yet), so pe is exactly what the first evaluation uses. Plain HS (G0_ref = 0) gets its
            // params back unchanged -> bit-for-bit.
            const auto pe_seed = hs_small_strain_params(models[mat].hs, 0.0);
            gs.pp = hs_initial_pp(pe_seed, sig_cp, ocr);
            gs.gamma_p = hs_initial_gamma_p(pe_seed, sig_cp);
        }
    }
}

}  // namespace katai::core
