#pragma once
// .k2d input validation -- the executable half of the input contract (Stage C).
// validate_project() walks a parsed model::Project and reports every violation
// as an Issue whose path names the exact JSON field (`materials[2].Rf`) and
// whose message states the rule the way an engineer would, echoing the
// offending value. The validator REPORTS and never repairs: silent correction
// is how silently-wrong numbers are born.
//
// Scope -- what belongs here and what deliberately does not:
//   * IN:  the schema's structural consistency (array sizes, index
//          cross-references, enum ranges the reader cannot check without
//          losing information) and parameter bounds that make a solve
//          indefensible on sight (E <= 0, Rf >= 1, kappa* >= lambda*).
//   * OUT: model x drainage x phase feasibility (e.g. Soft Soil Creep under a
//          Consolidation phase, NonPorous in a flow phase). The ENGINE refuses
//          those at solve time with its own tested messages -- one source of
//          truth; duplicating the rules here would let the two drift apart.
//   * OUT: geometry that only the mesher can judge (polygon self-intersection,
//          region overlap). The mesher's failure is the honest report there.
//
// Parameter-bound references: PLAXIS 2D Material Models Manual -- HS ranges and
// the Rf < 1 asymptote (par. 6.4), Soft Soil lambda* > kappa* (par. 10), SSC mu*
// (par. 11); van Genuchten (1980) requires n > 1; Rayleigh target ordering
// 0 < f1 < f2 per docs/references/dynamic-seismic-formulation.md.
//
// Every rule is pinned by tests/test_k2d_validator.cpp, which also prints the
// full message catalogue so engineer-readability stays reviewable.

#include <cstdio>
#include <string>
#include <vector>

#include <katai/io/issue.hpp>
#include <katai/model/project.hpp>

namespace katai::io {

namespace detail {

inline std::string at(const char* list, size_t i, const char* field = nullptr) {
    char b[128];
    if (field) std::snprintf(b, sizeof(b), "%s[%zu].%s", list, i, field);
    else std::snprintf(b, sizeof(b), "%s[%zu]", list, i);
    return b;
}
inline std::string num(double v) {   // echo a value the way the file stores it
    char b[32];
    std::snprintf(b, sizeof(b), "%g", v);
    return b;
}

// Phase rules are shared between `initial` and `phases[i]`; `base` is the JSON
// path prefix and the counts tie the activity vectors to the project's lists.
inline void check_phase(ValidationReport& r, const model::Phase& ph, const std::string& base,
                        bool is_initial, size_t n_poly, size_t n_struct, size_t n_load,
                        size_t n_disp) {
    using model::PhaseType;
    using model::SeismicWave;
    const auto path = [&base](const char* f) { return base + "." + f; };

    // Per-phase water conditions (staged dewatering): the polyline must be usable, or the
    // phase would silently fall back to the project's level and run an excavation that was
    // never dewatered -- the exact failure the version bump exists to prevent.
    if (ph.water_override) {
        if (ph.wx.size() != ph.wy.size())
            r.add(Severity::Error, path("wy"),
                  "the phase water line has " + std::to_string(ph.wx.size()) +
                      " x-value(s) but " + std::to_string(ph.wy.size()) + " y-value(s)");
        else if (ph.wx.size() < 2)
            r.add(Severity::Error, path("wx"),
                  "water_override is set but the phase water line has " +
                      std::to_string(ph.wx.size()) +
                      " point(s); at least two are needed to define a phreatic surface");
    } else if (!ph.wx.empty() || !ph.wy.empty()) {
        r.add(Severity::Warning, path("water_override"),
              "the phase carries a water line but water_override is false, so the project's "
              "water level is used and these points are ignored");
    }

    // Numerical controls (0 = let the program choose). A value that is set must be usable,
    // because the alternative is a run that quietly ignores the numerics the file was written
    // to carry -- which is the whole reason they are in the file.
    if (ph.tolerance < 0.0)
        r.add(Severity::Error, path("tol"),
              "a tolerated error cannot be negative (got " + num(ph.tolerance) + ")");
    else if (ph.tolerance >= 1.0)
        r.add(Severity::Error, path("tol"),
              "a tolerated error of " + num(ph.tolerance) +
                  " accepts a residual as large as the load itself; it is a relative error, so "
                  "it must be below 1 (PLAXIS's default is 0.01)");
    else if (ph.tolerance > 0.05)
        r.add(Severity::Warning, path("tol"),
              "a tolerated error of " + num(ph.tolerance) + " (" + num(100.0 * ph.tolerance) +
                  "%) is far looser than any published default (PLAXIS: 1%); the equilibrium "
                  "this phase reports may be a long way from equilibrium");
    if (ph.load_steps < 0)
        r.add(Severity::Error, path("loadsteps"),
              "the number of load increments cannot be negative (got " +
                  std::to_string(ph.load_steps) + "; 0 means the program chooses)");
    if (ph.max_iterations < 0)
        r.add(Severity::Error, path("maxiter"),
              "the iteration limit cannot be negative (got " +
                  std::to_string(ph.max_iterations) + "; 0 means the program chooses)");
    else if (ph.max_iterations > 0 && ph.max_iterations < 4)
        r.add(Severity::Warning, path("maxiter"),
              "an iteration limit of " + std::to_string(ph.max_iterations) +
                  " will make the solver cut the load step back almost every time; the run will "
                  "be slow rather than wrong");

    const int type = (int)ph.type;
    if (type < 0 || type >= 6) {
        // The driver treats an unknown type as Plastic, so an unchecked value
        // would not fail loudly -- it would solve the WRONG analysis.
        r.add(Severity::Error, path("type"),
              "unknown phase type " + std::to_string(type) +
                  "; this build knows 0..5 (Plastic, Safety, Consolidation, TransientFlow, "
                  "FullyCoupled, Dynamic)");
        return;   // the remaining rules key off the type
    }
    if (is_initial && ph.type != PhaseType::Plastic)
        r.add(Severity::Warning, path("type"),
              "the initial phase establishes the initial stress state (K0 procedure / gravity); "
              "its phase type is ignored");

    const bool timed = ph.type == PhaseType::Consolidation || ph.type == PhaseType::TransientFlow ||
                       ph.type == PhaseType::FullyCoupled || ph.type == PhaseType::Dynamic;
    if (timed && !(ph.duration > 0.0))
        r.add(Severity::Error, path("duration"),
              "a time-dependent phase needs a positive time interval (got " + num(ph.duration) +
                  (ph.type == PhaseType::Dynamic ? " s)" : " day)"));
    if (timed && ph.time_steps < 1)
        r.add(Severity::Error, path("steps"),
              "a time-dependent phase needs at least one time step (got " +
                  std::to_string(ph.time_steps) + ")");

    if (ph.type == PhaseType::Dynamic) {
        const int wave = (int)ph.seismic_wave;
        if (wave < 0 || wave >= 3) {
            r.add(Severity::Error, path("seiswave"),
                  "unknown seismic waveform " + std::to_string(wave) +
                      "; this build knows 0..2 (Harmonic, Ricker, Record)");
        } else if (ph.seismic_wave == SeismicWave::Record) {
            if (ph.accel_record.empty())
                r.add(Severity::Error, path("rec"),
                      "the waveform is an accelerogram record but the record is empty");
            if (!(ph.record_dt > 0.0))
                r.add(Severity::Error, path("recdt"),
                      "the accelerogram sampling interval must be positive (got " +
                          num(ph.record_dt) + " s)");
        } else if (!(ph.seismic_freq > 0.0)) {
            r.add(Severity::Error, path("seisfreq"),
                  "a Harmonic/Ricker input needs a positive dominant frequency (got " +
                      num(ph.seismic_freq) + " Hz)");
        }
        if (ph.damping_ratio < 0.0)
            r.add(Severity::Error, path("damp"),
                  "the damping ratio cannot be negative (got " + num(ph.damping_ratio) + ")");
        else if (ph.damping_ratio >= 1.0)
            r.add(Severity::Warning, path("damp"),
                  "a damping ratio of " + num(ph.damping_ratio) +
                      " is at or above critical damping; typical soil values are 0.01-0.10");
        if (!(ph.rayleigh_f1 > 0.0) || !(ph.rayleigh_f2 > ph.rayleigh_f1))
            r.add(Severity::Error, path("rayf1"),
                  "the Rayleigh damping targets must satisfy 0 < f1 < f2 (got f1 = " +
                      num(ph.rayleigh_f1) + " Hz, f2 = " + num(ph.rayleigh_f2) + " Hz)");
        if (ph.site_class < 0 || ph.site_class > 4)
            r.add(Severity::Error, path("siteclass"),
                  "the TBDY site class must be 0..4 (ZA..ZE); got " +
                      std::to_string(ph.site_class));
        if (!(ph.tbdy_ss > 0.0) || !(ph.tbdy_s1 > 0.0))
            r.add(Severity::Error, path("tbdyss"),
                  "the TBDY spectral coefficients S_S and S_1 must be positive (got S_S = " +
                      num(ph.tbdy_ss) + ", S_1 = " + num(ph.tbdy_s1) + ")");
        if (ph.ec8_enabled) {
            if (ph.ec8_ground < 0 || ph.ec8_ground > 4)
                r.add(Severity::Error, path("ec8gnd"),
                      "the EC8 ground type must be 0..4 (A..E); got " +
                          std::to_string(ph.ec8_ground));
            if (ph.ec8_type != 0 && ph.ec8_type != 1)
                r.add(Severity::Error, path("ec8typ"),
                      "the EC8 spectrum type must be 0 (Type 1) or 1 (Type 2); got " +
                          std::to_string(ph.ec8_type));
            if (ph.ec8_agr < 0.0)
                r.add(Severity::Error, path("ec8agr"),
                      "the EC8 reference peak ground acceleration cannot be negative (got " +
                          num(ph.ec8_agr) + " g)");
            if (!(ph.ec8_gamma > 0.0))
                r.add(Severity::Error, path("ec8gi"),
                      "the EC8 importance factor must be positive (got " + num(ph.ec8_gamma) + ")");
        }
    }

    const int design = (int)ph.design_approach;
    if (design < 0 || design >= 7)
        r.add(Severity::Error, path("design"),
              "unknown design approach " + std::to_string(design) +
                  "; this build knows 0..6 (None, EC7 DA1-C1/DA1-C2/DA2/DA3, TBDY static/seismic)");

    // A longer-than-the-list activity vector is silently ignored entry by entry
    // (the "missing = active" rule reads only the first N) -- surface it.
    const auto check_activity = [&](const char* field, size_t have, size_t want,
                                    const char* what) {
        if (have > want)
            r.add(Severity::Warning, path(field),
                  "the phase lists " + std::to_string(have) + " " + what +
                      " activity flags but the project has " + std::to_string(want) + "; the " +
                      std::to_string(have - want) + " extra flag(s) are ignored");
    };
    check_activity("poly", ph.poly_active.size(), n_poly, "soil-region");
    check_activity("struct", ph.struct_active.size(), n_struct, "structural-element");
    check_activity("load", ph.load_active.size(), n_load, "load");
    check_activity("disp", ph.disp_active.size(), n_disp, "prescribed-displacement");

    // Prescribed displacements are a STATIC (Plastic) capability in this build. A
    // time-dependent or Safety phase with one active would silently need semantics this
    // build does not have (creep under held displacement, phi-c with imposed motion),
    // so it is refused by name rather than half-applied.
    bool any_disp = false;
    for (size_t i = 0; i < n_disp; ++i) any_disp = any_disp || ph.active_disp(i);
    if (any_disp && !is_initial && ph.type != PhaseType::Plastic)
        r.add(Severity::Error, path("disp"),
              "prescribed displacements are only supported in Plastic (staged construction) "
              "phases in this build; deactivate them in this phase or change the phase type");
}

inline void check_material(ValidationReport& r, const model::Material& m, size_t i) {
    using model::Drainage;
    using model::SoilModel;
    const auto path = [i](const char* f) { return at("materials", i, f); };
    const std::string who = "\"" + m.name + "\": ";

    // -- General ---------------------------------------------------------------
    // Zero unit weight is legal on purpose: weightless verification materials
    // (Prandtl, Boussinesq oracles) are part of the validated record.
    if (m.gamma_unsat < 0.0)
        r.add(Severity::Error, path("gamma_unsat"),
              who + "the unit weight cannot be negative (got " + num(m.gamma_unsat) + " kN/m3)");
    if (m.gamma_sat < 0.0)
        r.add(Severity::Error, path("gamma_sat"),
              who + "the unit weight cannot be negative (got " + num(m.gamma_sat) + " kN/m3)");
    if (m.drainage != Drainage::NonPorous && m.gamma_sat < m.gamma_unsat)
        r.add(Severity::Warning, path("gamma_sat"),
              who + "the saturated unit weight (" + num(m.gamma_sat) +
                  " kN/m3) is below the unsaturated one (" + num(m.gamma_unsat) +
                  " kN/m3), which is physically unusual");
    if (m.e_init < 0.0)
        r.add(Severity::Error, path("e_init"),
              who + "the initial void ratio cannot be negative (got " + num(m.e_init) + ")");

    // -- Stiffness and strength, by the fields the chosen model actually reads --
    const bool le = m.model == SoilModel::LinearElastic;
    const bool mc = m.model == SoilModel::MohrCoulomb;
    const bool hs = m.model == SoilModel::HardeningSoil || m.model == SoilModel::HSsmall;
    const bool ss = m.model == SoilModel::SoftSoil || m.model == SoilModel::SoftSoilCreep;

    if ((le || mc) && !(m.E > 0.0))
        r.add(Severity::Error, path("E"),
              who + "Young's modulus must be positive (got " + num(m.E) + " kN/m2)");
    if ((le || mc) && !(m.nu > -1.0 && m.nu < 0.5))
        r.add(Severity::Error, path("nu"),
              who + "Poisson's ratio must lie in (-1, 0.5); at 0.5 the material is "
                    "incompressible and the stiffness matrix is singular (got " +
                  num(m.nu) + ")");

    if (!le) {   // every plastic model carries a Mohr-Coulomb strength envelope
        if (m.phi < 0.0 || m.phi >= 90.0)
            r.add(Severity::Error, path("phi"),
                  who + "the friction angle must lie in [0, 90) degrees (got " + num(m.phi) + ")");
        if (m.c < 0.0)
            r.add(Severity::Error, path("c"),
                  who + "cohesion cannot be negative (got " + num(m.c) + " kN/m2)");
        if (m.psi > m.phi)
            r.add(Severity::Error, path("psi"),
                  who + "the dilatancy angle cannot exceed the friction angle (psi = " +
                      num(m.psi) + " > phi = " + num(m.phi) + " degrees)");
        if (m.c == 0.0 && m.phi == 0.0)
            r.add(Severity::Error, path("c"),
                  who + "the material has no shear strength (c = 0 and phi = 0)");
        // The dilatancy cut-off needs room to act: a critical void ratio at or below the in-situ
        // one means the soil is already at critical density and the cut-off would fire on the
        // first increment -- almost always a typo rather than an intention.
        if (m.dilatancy_cutoff && !(m.e_max > m.e_init))
            r.add(Severity::Error, path("e_max"),
                  who + "the dilatancy cut-off needs a critical void ratio above the in-situ one "
                        "(got e_max = " + num(m.e_max) + ", e_init = " + num(m.e_init) + ")");
        if (m.dilatancy_cutoff && !(m.psi > 0.0))
            r.add(Severity::Warning, path("dilatancy_cutoff"),
                  who + "the dilatancy cut-off is set but the dilatancy angle is " + num(m.psi) +
                      " deg, so there is no dilation to cut off");
        if (m.tensile_strength < 0.0)
            r.add(Severity::Error, path("tensile_strength"),
                  who + "the tensile strength cannot be negative (got " +
                      num(m.tensile_strength) + " kN/m2)");
        if (m.drainage == Drainage::UndrainedB) {
            if (!(m.c > 0.0))
                r.add(Severity::Error, path("c"),
                      who + "Undrained (B) reads c as the undrained shear strength su, "
                            "which must be positive (got " + num(m.c) + " kN/m2)");
            if (m.phi != 0.0)
                r.add(Severity::Warning, path("phi"),
                      who + "Undrained (B) is a Tresca envelope: phi is forced to 0 and the "
                            "entered value (" + num(m.phi) + " degrees) is ignored");
        }
    }

    if (hs) {
        if (!(m.E50ref > 0.0))
            r.add(Severity::Error, path("E50ref"),
                  who + "E50_ref must be positive (got " + num(m.E50ref) + " kN/m2)");
        if (!(m.Eoedref > 0.0))
            r.add(Severity::Error, path("Eoedref"),
                  who + "Eoed_ref must be positive (got " + num(m.Eoedref) + " kN/m2)");
        if (!(m.Eurref > 0.0))
            r.add(Severity::Error, path("Eurref"),
                  who + "Eur_ref must be positive (got " + num(m.Eurref) + " kN/m2)");
        else if (m.Eurref <= m.E50ref)
            r.add(Severity::Error, path("Eurref"),
                  who + "the unload-reload stiffness must exceed E50 (Eur = " + num(m.Eurref) +
                      " <= E50 = " + num(m.E50ref) + " kN/m2)");
        else if (m.Eurref < 2.0 * m.E50ref)
            r.add(Severity::Warning, path("Eurref"),
                  who + "Eur = " + num(m.Eurref) + " kN/m2 is below the recommended Eur >= 2 E50 "
                        "(E50 = " + num(m.E50ref) + " kN/m2); the hyperbolic law may misbehave");
        if (m.m < 0.0 || m.m > 1.0)
            r.add(Severity::Error, path("m"),
                  who + "the stiffness power m must lie in [0, 1] (got " + num(m.m) + ")");
        if (!(m.nu_ur > -1.0 && m.nu_ur < 0.5))
            r.add(Severity::Error, path("nu_ur"),
                  who + "the unload-reload Poisson's ratio must lie in (-1, 0.5) (got " +
                      num(m.nu_ur) + ")");
        if (!(m.p_ref > 0.0))
            r.add(Severity::Error, path("p_ref"),
                  who + "the reference pressure must be positive (got " + num(m.p_ref) +
                      " kN/m2)");
        if (!(m.Rf > 0.0 && m.Rf < 1.0))
            r.add(Severity::Error, path("Rf"),
                  who + "the failure ratio must lie in (0, 1); at Rf = 1 the hyperbola's "
                        "asymptote equals the failure load and q_a becomes singular (got " +
                      num(m.Rf) + ")");
        if (!m.k0nc_auto && !(m.k0nc > 0.0))
            r.add(Severity::Error, path("k0nc"),
                  who + "K0_nc must be positive (got " + num(m.k0nc) + ")");
        if (m.model == SoilModel::HSsmall) {
            if (!(m.G0ref > 0.0))
                r.add(Severity::Error, path("G0ref"),
                      who + "G0_ref must be positive (got " + num(m.G0ref) + " kN/m2)");
            if (!(m.gamma07 > 0.0))
                r.add(Severity::Error, path("gamma07"),
                      who + "gamma_0.7 must be positive (got " + num(m.gamma07) + ")");
            const double gur = m.Eurref / (2.0 * (1.0 + m.nu_ur));
            if (m.G0ref > 0.0 && m.Eurref > 0.0 && m.G0ref < gur)
                r.add(Severity::Warning, path("G0ref"),
                      who + "the small-strain modulus G0 = " + num(m.G0ref) +
                          " kN/m2 is below the unload-reload shear modulus Gur = Eur/(2(1+nu_ur))"
                          " = " + num(gur) + " kN/m2; the small-strain overlay would soften "
                          "rather than stiffen");
        }
    }

    if (ss) {
        if (!(m.lam_star > 0.0))
            r.add(Severity::Error, path("lamstar"),
                  who + "the modified compression index lambda* must be positive (got " +
                      num(m.lam_star) + ")");
        if (!(m.kap_star > 0.0))
            r.add(Severity::Error, path("kapstar"),
                  who + "the modified swelling index kappa* must be positive (got " +
                      num(m.kap_star) + ")");
        else if (m.lam_star > 0.0 && m.kap_star >= m.lam_star)
            r.add(Severity::Error, path("kapstar"),
                  who + "the swelling index must be below the compression index (kappa* = " +
                      num(m.kap_star) + " >= lambda* = " + num(m.lam_star) + ")");
        if (!(m.nu_ur > -1.0 && m.nu_ur < 0.5))
            r.add(Severity::Error, path("nu_ur"),
                  who + "the unload-reload Poisson's ratio must lie in (-1, 0.5) (got " +
                      num(m.nu_ur) + ")");
        if (m.model == SoilModel::SoftSoilCreep && !(m.mu_star > 0.0))
            r.add(Severity::Error, path("mustar"),
                  who + "the modified creep index mu* must be positive (got " + num(m.mu_star) +
                      ")");
    }

    // -- Groundwater -----------------------------------------------------------
    if (m.kx < 0.0 || m.ky < 0.0)
        r.add(Severity::Error, path(m.kx < 0.0 ? "kx" : "ky"),
              who + "permeability cannot be negative (got kx = " + num(m.kx) + ", ky = " +
                  num(m.ky) + " m/day); use 0 for impermeable");
    if (!(m.gw_ga > 0.0))
        r.add(Severity::Error, path("gw_ga"),
              who + "the van Genuchten g_a (inverse air-entry) must be positive (got " +
                  num(m.gw_ga) + " 1/m)");
    if (!(m.gw_gn > 1.0))
        r.add(Severity::Error, path("gw_gn"),
              who + "the van Genuchten g_n must exceed 1 (got " + num(m.gw_gn) +
                  "); at n <= 1 the retention curve is undefined");
    if (m.gw_Sres < 0.0 || m.gw_Sres >= 1.0)
        r.add(Severity::Error, path("gw_Sres"),
              who + "the residual saturation must lie in [0, 1) (got " + num(m.gw_Sres) + ")");

    // -- Interfaces / initial state --------------------------------------------
    if (!m.rinter_rigid && !(m.Rinter > 0.0 && m.Rinter <= 1.0))
        r.add(Severity::Error, path("Rinter"),
              who + "the interface strength factor must lie in (0, 1] (got " + num(m.Rinter) +
                  ")");
    if (!m.k0_auto && !(m.k0 > 0.0))
        r.add(Severity::Error, path("k0"),
              who + "K0 must be positive (got " + num(m.k0) + ")");
    if (m.oc_mode < 0 || m.oc_mode > 2)
        r.add(Severity::Error, path("oc_mode"),
              who + "the stress-history mode must be 0 (none), 1 (OCR) or 2 (POP); got " +
                  std::to_string(m.oc_mode));
    if (m.oc_mode == 1 && m.OCR < 1.0)
        r.add(Severity::Error, path("OCR"),
              who + "OCR below 1 means under-consolidated soil, which the K0 procedure cannot "
                    "seed; enter 1 for normally consolidated (got " + num(m.OCR) + ")");
    if (m.oc_mode == 2 && m.POP < 0.0)
        r.add(Severity::Error, path("POP"),
              who + "the pre-overburden pressure cannot be negative (got " + num(m.POP) +
                  " kN/m2)");
}

}  // namespace detail

// Validate a parsed project against the input contract. Pure and total: never
// throws, never mutates, reports everything it finds in file order.
inline ValidationReport validate_project(const model::Project& p) {
    using model::LoadKind;
    using model::StructKind;
    using detail::at;
    using detail::num;

    ValidationReport r;

    // -- Domain ----------------------------------------------------------------
    if (!(p.x_max > p.x_min))
        r.add(Severity::Error, "x_max",
              "the model must span x_min < x_max (got x_min = " + num(p.x_min) +
                  ", x_max = " + num(p.x_max) + " m)");
    if (!(p.y_max > p.y_min))
        r.add(Severity::Error, "y_max",
              "the model must span y_min < y_max (got y_min = " + num(p.y_min) +
                  ", y_max = " + num(p.y_max) + " m)");
    if (p.axisymmetric && p.x_min < 0.0)
        r.add(Severity::Error, "x_min",
              "an axisymmetric model measures x as the radius from the symmetry axis, so "
              "x_min must be >= 0 (got " + num(p.x_min) + " m)");

    // -- Mesh and initial procedure --------------------------------------------
    if (!(p.mesh.elem_size > 0.0))
        r.add(Severity::Error, "mesh.elem_size",
              "the target element size must be positive (got " + num(p.mesh.elem_size) + " m)");
    if (p.mesh.order != 6 && p.mesh.order != 15)
        r.add(Severity::Error, "mesh.order",
              "the element order must be 6 (quadratic tri6) or 15 (quartic tri15); got " +
                  std::to_string(p.mesh.order));
    const int proc = (int)p.initial_procedure;
    if (proc < 0 || proc >= 3)
        r.add(Severity::Error, "initial_procedure",
              "unknown initial procedure " + std::to_string(proc) +
                  "; this build knows 0..2 (K0 procedure, Gravity loading, Safety)");
    else if (p.initial_procedure == model::InitialProcedure::Safety && !p.phases.empty())
        r.add(Severity::Warning, "initial_procedure",
              "Safety as the initial procedure computes the factor of safety of the initial "
              "state only; with staged phases present, add a Safety phase at the point of "
              "interest instead");

    // -- Water -----------------------------------------------------------------
    if (p.has_water && p.wx.size() != p.wy.size())
        r.add(Severity::Error, "wy",
              "the water line has " + std::to_string(p.wx.size()) + " x-value(s) but " +
                  std::to_string(p.wy.size()) + " y-value(s)");
    else if (p.has_water && p.wx.size() < 2)
        r.add(Severity::Warning, "wx",
              "has_water is set but the water line has " + std::to_string(p.wx.size()) +
                  " point(s); at least two are needed, so the model is treated as DRY");

    // -- Material data sets ----------------------------------------------------
    for (size_t i = 0; i < p.materials.size(); ++i) detail::check_material(r, p.materials[i], i);

    for (size_t i = 0; i < p.plates.size(); ++i) {
        const auto& m = p.plates[i];
        const std::string who = "\"" + m.name + "\": ";
        if (!(m.EA > 0.0))
            r.add(Severity::Error, at("plates", i, "EA"),
                  who + "the axial stiffness must be positive (got " + num(m.EA) + " kN/m)");
        if (!(m.EI > 0.0))
            r.add(Severity::Error, at("plates", i, "EI"),
                  who + "the bending stiffness must be positive (got " + num(m.EI) +
                      " kN m2/m)");
        if (m.w < 0.0)
            r.add(Severity::Error, at("plates", i, "w"),
                  who + "the plate weight cannot be negative (got " + num(m.w) + " kN/m/m)");
        if (m.nu < 0.0 || m.nu >= 0.5)
            r.add(Severity::Error, at("plates", i, "nu"),
                  who + "Poisson's ratio must lie in [0, 0.5) (got " + num(m.nu) + ")");
        if (m.Mp < 0.0 || m.Np < 0.0)
            r.add(Severity::Error, at("plates", i, m.Mp < 0.0 ? "Mp" : "Np"),
                  who + "a plastic capacity cannot be negative (got Mp = " + num(m.Mp) +
                      " kN m/m, Np = " + num(m.Np) + " kN/m); use 0 for elastic");
        if (m.elastoplastic && m.Mp == 0.0 && m.Np == 0.0)
            r.add(Severity::Warning, at("plates", i, "elastoplastic"),
                  who + "elastoplastic is set but Mp = Np = 0 means no capacity is checked; "
                        "the plate behaves elastically");
    }
    for (size_t i = 0; i < p.anchors.size(); ++i) {
        const auto& m = p.anchors[i];
        const std::string who = "\"" + m.name + "\": ";
        if (!(m.EA > 0.0))
            r.add(Severity::Error, at("anchors", i, "EA"),
                  who + "the axial stiffness must be positive (got " + num(m.EA) + " kN)");
        if (m.Fmax_tens < 0.0 || m.Fmax_comp < 0.0)
            r.add(Severity::Error, at("anchors", i, m.Fmax_tens < 0.0 ? "Fmax_tens" : "Fmax_comp"),
                  who + "an anchor capacity cannot be negative (got tension " + num(m.Fmax_tens) +
                      ", compression " + num(m.Fmax_comp) + " kN); use 0 for unlimited");
        if (!(m.Lspacing > 0.0))
            r.add(Severity::Error, at("anchors", i, "Lspacing"),
                  who + "the out-of-plane spacing must be positive (got " + num(m.Lspacing) +
                      " m)");
        if (m.elastoplastic && m.Fmax_tens == 0.0 && m.Fmax_comp == 0.0)
            r.add(Severity::Warning, at("anchors", i, "elastoplastic"),
                  who + "elastoplastic is set but both capacities are 0 (unlimited); the anchor "
                        "behaves elastically");
        // A lock-off force is a TENSION: a negative one would mean the anchor was installed
        // pushing the wall into the soil, which is not what a prestressed anchor or strut does.
        if (m.prestress < 0.0)
            r.add(Severity::Error, at("anchors", i, "prestress"),
                  who + "the lock-off force is tension-positive and cannot be negative (got " +
                      num(m.prestress) + " kN)");
        // Locking off beyond the capacity is not a modelling preference: the anchor would yield
        // on installation, so the force the file asks for cannot exist.
        if (m.elastoplastic && m.Fmax_tens > 0.0 && m.prestress > m.Fmax_tens)
            r.add(Severity::Error, at("anchors", i, "prestress"),
                  who + "the lock-off force " + num(m.prestress) +
                      " kN exceeds the tension capacity " + num(m.Fmax_tens) +
                      " kN; the anchor would yield as it is installed");
    }
    for (size_t i = 0; i < p.geogrids.size(); ++i) {
        const auto& m = p.geogrids[i];
        const std::string who = "\"" + m.name + "\": ";
        if (!(m.EA > 0.0))
            r.add(Severity::Error, at("geogrids", i, "EA"),
                  who + "the axial stiffness must be positive (got " + num(m.EA) + " kN/m)");
        if (m.Np < 0.0)
            r.add(Severity::Error, at("geogrids", i, "Np"),
                  who + "the tension cap cannot be negative (got " + num(m.Np) +
                      " kN/m); use 0 for unlimited");
        if (m.elastoplastic && m.Np == 0.0)
            r.add(Severity::Warning, at("geogrids", i, "elastoplastic"),
                  who + "elastoplastic is set but Np = 0 (unlimited); the geogrid behaves "
                        "elastically");
    }
    for (size_t i = 0; i < p.embedded.size(); ++i) {
        const auto& m = p.embedded[i];
        const std::string who = "\"" + m.name + "\": ";
        if (!(m.E > 0.0))
            r.add(Severity::Error, at("embedded", i, "E"),
                  who + "Young's modulus must be positive (got " + num(m.E) + " kN/m2)");
        if (m.gamma < 0.0)
            r.add(Severity::Error, at("embedded", i, "gamma"),
                  who + "the unit weight cannot be negative (got " + num(m.gamma) + " kN/m3)");
        if (!(m.diameter > 0.0))
            r.add(Severity::Error, at("embedded", i, "diameter"),
                  who + "the pile diameter must be positive (got " + num(m.diameter) + " m)");
        if (!(m.Lspacing > 0.0))
            r.add(Severity::Error, at("embedded", i, "Lspacing"),
                  who + "the out-of-plane spacing must be positive (got " + num(m.Lspacing) +
                      " m)");
        if (m.Tskin_max < 0.0 || m.Fmax_base < 0.0)
            r.add(Severity::Error,
                  at("embedded", i, m.Tskin_max < 0.0 ? "Tskin_max" : "Fmax_base"),
                  who + "a resistance cannot be negative (got skin " + num(m.Tskin_max) +
                      " kN/m, base " + num(m.Fmax_base) + " kN); use 0 for unlimited");
    }

    // -- Soil regions ----------------------------------------------------------
    if (p.polygons.empty())
        r.add(Severity::Error, "polygons",
              "the project has no soil regions; there is nothing to mesh");
    for (size_t i = 0; i < p.polygons.size(); ++i) {
        const auto& P = p.polygons[i];
        const std::string who = "\"" + P.name + "\": ";
        const size_t nv = P.x.size();
        if (P.x.size() != P.y.size()) {
            r.add(Severity::Error, at("polygons", i, "y"),
                  who + "the region has " + std::to_string(P.x.size()) + " x-value(s) but " +
                      std::to_string(P.y.size()) + " y-value(s)");
        } else if (nv < 3) {
            r.add(Severity::Error, at("polygons", i, "x"),
                  who + "a soil region needs at least three vertices (got " +
                      std::to_string(nv) + ")");
        }
        if (P.material < 0)
            r.add(Severity::Error, at("polygons", i, "material"),
                  who + "no soil material assigned");
        else if ((size_t)P.material >= p.materials.size())
            r.add(Severity::Error, at("polygons", i, "material"),
                  who + "material index " + std::to_string(P.material) +
                      " does not exist; the project has " + std::to_string(p.materials.size()) +
                      " soil material(s)");
        if (!(P.coarseness > 0.0))
            r.add(Severity::Error, at("polygons", i, "coarseness"),
                  who + "the mesh coarseness factor must be positive (got " +
                      num(P.coarseness) + ")");
        const auto check_edges = [&](const char* field, const std::vector<int>& v, int lo,
                                     int hi, const char* what) {
            if (v.empty()) return;
            if (v.size() != nv) {
                r.add(Severity::Error, at("polygons", i, field),
                      who + "the region has " + std::to_string(nv) + " edge(s) but " +
                          std::to_string(v.size()) + " " + what + " value(s)");
                return;
            }
            for (size_t j = 0; j < v.size(); ++j)
                if (v[j] < lo || v[j] >= hi) {
                    r.add(Severity::Error,
                          at("polygons", i, field) + "[" + std::to_string(j) + "]",
                          who + std::string(what) + " value " + std::to_string(v[j]) +
                              " is outside " + std::to_string(lo) + ".." +
                              std::to_string(hi - 1));
                    break;   // one report per array keeps the output readable
                }
        };
        check_edges("edge_bc", P.edge_bc, 0, 5, "boundary-condition");
        check_edges("edge_flow", P.edge_flow, 0, 4, "flow-boundary");
        // A prescribed-flux edge needs a value to prescribe: without the array the edge would
        // silently behave as closed, which is a different problem from the one asked for.
        {
            bool has_flux_edge = false;
            for (size_t j = 0; j < P.edge_flow.size(); ++j)
                if (P.edge_flow[j] == (int)model::FlowBCType::Flux) has_flux_edge = true;
            if (has_flux_edge && P.edge_flux.size() != nv)
                r.add(Severity::Error, at("polygons", i, "edge_flux"),
                      who + "an edge prescribes a flux, so edge_flux must list one value per edge "
                            "(got " + std::to_string(P.edge_flux.size()) + " of " +
                            std::to_string(nv) + ")");
        }
        if (!P.edge_flow.empty() && P.edge_head.size() != nv)
            r.add(Severity::Error, at("polygons", i, "edge_head"),
                  who + "edge_flow is present, so edge_head must list one head per edge (got " +
                      std::to_string(P.edge_head.size()) + " of " + std::to_string(nv) + ")");
    }

    // -- Structural elements ---------------------------------------------------
    for (size_t i = 0; i < p.structs.size(); ++i) {
        const auto& s = p.structs[i];
        const std::string who = "\"" + s.name + "\": ";
        const int kind = (int)s.kind;
        if (kind < 0 || kind >= 5) {
            r.add(Severity::Error, at("structs", i, "kind"),
                  who + "unknown element kind " + std::to_string(kind) +
                      "; this build knows 0..4 (Plate, Anchor, Geogrid, Embedded beam, "
                      "Interface)");
            continue;   // the material rule keys off the kind
        }
        if (s.x1 == s.x2 && s.y1 == s.y2)
            r.add(Severity::Error, at("structs", i, "x2"),
                  who + "the element has zero length (both endpoints at (" + num(s.x1) + ", " +
                      num(s.y1) + "))");
        const auto need_material = [&](size_t nsets, const char* list) {
            if (s.material < 0)
                r.add(Severity::Error, at("structs", i, "material"),
                      who + "no " + std::string(list) + " material set assigned");
            else if ((size_t)s.material >= nsets)
                r.add(Severity::Error, at("structs", i, "material"),
                      who + std::string(list) + " material index " +
                          std::to_string(s.material) + " does not exist; the project has " +
                          std::to_string(nsets) + " set(s)");
        };
        switch (s.kind) {
            case StructKind::Plate: need_material(p.plates.size(), "plate"); break;
            case StructKind::Anchor: need_material(p.anchors.size(), "anchor"); break;
            case StructKind::Geogrid: need_material(p.geogrids.size(), "geogrid"); break;
            case StructKind::EmbeddedBeam: need_material(p.embedded.size(), "embedded-beam"); break;
            case StructKind::Interface:
                if (s.material >= 0)
                    r.add(Severity::Warning, at("structs", i, "material"),
                          who + "an interface takes its strength from the adjacent soil (or "
                                "iface_material); the material field is ignored");
                break;
        }
        if (s.iface_material >= 0 && (size_t)s.iface_material >= p.materials.size())
            r.add(Severity::Error, at("structs", i, "iface_material"),
                  who + "interface material index " + std::to_string(s.iface_material) +
                      " does not exist; the project has " + std::to_string(p.materials.size()) +
                      " soil material(s)");
        if (!(s.coarseness > 0.0))
            r.add(Severity::Error, at("structs", i, "coarseness"),
                  who + "the mesh coarseness factor must be positive (got " + num(s.coarseness) +
                      ")");
    }

    // -- Loads -----------------------------------------------------------------
    for (size_t i = 0; i < p.loads.size(); ++i) {
        const auto& L = p.loads[i];
        const std::string who = "\"" + L.name + "\": ";
        const int kind = (int)L.kind;
        if (kind < 0 || kind >= 2) {
            r.add(Severity::Error, at("loads", i, "kind"),
                  who + "unknown load kind " + std::to_string(kind) +
                      "; this build knows 0 (Point) and 1 (Distributed)");
            continue;
        }
        if (L.kind == LoadKind::Distributed && L.x1 == L.x2 && L.y1 == L.y2)
            r.add(Severity::Error, at("loads", i, "x2"),
                  who + "a distributed load needs a line of positive length (both endpoints at "
                        "(" + num(L.x1) + ", " + num(L.y1) + "))");
        if (!(L.coarseness > 0.0))
            r.add(Severity::Error, at("loads", i, "coarseness"),
                  who + "the mesh coarseness factor must be positive (got " + num(L.coarseness) +
                      ")");
    }

    // -- Prescribed displacements ----------------------------------------------
    for (size_t i = 0; i < p.disps.size(); ++i) {
        const auto& D = p.disps[i];
        const std::string who = "\"" + D.name + "\": ";
        if (D.x1 == D.x2 && D.y1 == D.y2)
            r.add(Severity::Error, at("disps", i, "x2"),
                  who + "a prescribed displacement needs a line of positive length (both "
                        "endpoints at (" + num(D.x1) + ", " + num(D.y1) + "))");
        if (!D.set_ux && !D.set_uy)
            r.add(Severity::Error, at("disps", i, "set_uy"),
                  who + "neither component is prescribed; set set_ux and/or set_uy (a set "
                        "component with value 0 is a support line)");
        if (!(D.coarseness > 0.0))
            r.add(Severity::Error, at("disps", i, "coarseness"),
                  who + "the mesh coarseness factor must be positive (got " + num(D.coarseness) +
                      ")");
    }
    // The initial phase is the K0 / gravity / Safety PROCEDURE, not a staged solve: only
    // gravity loading actually equilibrates imposed displacements. Refuse the other two by
    // name -- K0 would set stresses AROUND the constraint and Safety would reduce phi-c
    // under an imposed motion neither procedure defines.
    if (p.initial_procedure != model::InitialProcedure::GravityLoading) {
        bool any0 = false;
        for (size_t i = 0; i < p.disps.size(); ++i) any0 = any0 || p.initial.active_disp(i);
        if (any0 && !p.disps.empty())
            r.add(Severity::Error, "initial.disp",
                  "prescribed displacements active in the initial phase need the gravity-loading "
                  "initial procedure; the K0 procedure and Safety do not solve for them "
                  "(deactivate them initially and activate them in a staged phase)");
    }

    // -- Phases ----------------------------------------------------------------
    detail::check_phase(r, p.initial, "initial", true, p.polygons.size(), p.structs.size(),
                        p.loads.size(), p.disps.size());
    for (size_t i = 0; i < p.phases.size(); ++i)
        detail::check_phase(r, p.phases[i], at("phases", i), false, p.polygons.size(),
                            p.structs.size(), p.loads.size(), p.disps.size());

    return r;
}

}  // namespace katai::io
