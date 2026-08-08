// katai._core -- the nanobind binding of the published facade (roadmap section 7.2,
// layer 4). This translation unit includes exactly ONE project header,
// <katai/api/katai.hpp>: the binding is a front end, and the section 7.3 boundary
// (total over the input contract, closed over the implementation) is inherited from
// the facade rather than re-argued here. Solver internals, element state and Gauss
// points do not appear because the facade does not export them.
//
// Layering (decision 2026-08-02): `_core` is the mechanical 1:1 layer -- schema
// structs, enums, file IO, validation, the Job -- and the engineer-facing `katai`
// package (pure Python, katai/__init__.py) builds its DSL on top. Units are the
// schema's own, fixed: kN, m, day (values pass through verbatim; every docstring
// states the unit).
//
// Completeness is mechanical, not remembered: test_python_schema_coverage walks
// every key the .k2d writer emits and fails unless the binding reads it
// faithfully and writes it back -- a schema field without a Python attribute
// cannot land silently.

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <katai/api/katai.hpp>

namespace nb = nanobind;
using namespace nb::literals;
namespace api = katai::api;

namespace {

// The validator's report rendered the way the engineer reads it -- reused by the
// exception text so a refusal in Python carries the same words as the CLI.
std::string report_text(const api::ValidationReport& r) { return r.to_string(); }

// Every material family carries a float[3] display colour; a raw array does not
// bind, so each class gets the same list<->array property.
template <class T>
void bind_color(nb::class_<T>& cls) {
    cls.def_prop_rw(
        "color",
        [](const T& v) { return std::vector<double>(v.color, v.color + 3); },
        [](T& v, const std::vector<double>& c) {
            if (c.size() != 3) throw nb::value_error("color must be [r, g, b] with 0..1 components");
            for (int i = 0; i < 3; ++i) v.color[i] = (float)c[i];
        },
        "display colour [r, g, b], components 0..1");
}

}  // namespace

NB_MODULE(_core, m) {
    m.doc() = "KATAI 2D core bindings: the published facade, 1:1. Units are fixed: "
              "kN, m, day. The engineer-facing API lives in the `katai` package.";

    // ----------------------------------------------------------------- enums --
    nb::enum_<api::SoilModel>(m, "SoilModel", nb::is_arithmetic())
        .value("LinearElastic", api::SoilModel::LinearElastic)
        .value("MohrCoulomb", api::SoilModel::MohrCoulomb)
        .value("HardeningSoil", api::SoilModel::HardeningSoil)
        .value("HSsmall", api::SoilModel::HSsmall)
        .value("SoftSoil", api::SoilModel::SoftSoil)
        .value("SoftSoilCreep", api::SoilModel::SoftSoilCreep);
    nb::enum_<api::Drainage>(m, "Drainage", nb::is_arithmetic())
        .value("Drained", api::Drainage::Drained)
        .value("Undrained", api::Drainage::Undrained)
        .value("NonPorous", api::Drainage::NonPorous)
        .value("UndrainedB", api::Drainage::UndrainedB);
    nb::enum_<api::PhaseType>(m, "PhaseType", nb::is_arithmetic())
        .value("Plastic", api::PhaseType::Plastic)
        .value("Safety", api::PhaseType::Safety)
        .value("Consolidation", api::PhaseType::Consolidation)
        .value("TransientFlow", api::PhaseType::TransientFlow)
        .value("FullyCoupled", api::PhaseType::FullyCoupled)
        .value("Dynamic", api::PhaseType::Dynamic);
    nb::enum_<api::BCType>(m, "BCType", nb::is_arithmetic())
        .value("Free", api::BCType::Free)
        .value("NormallyFixed", api::BCType::NormallyFixed)
        .value("HorizontallyFixed", api::BCType::HorizontallyFixed)
        .value("VerticallyFixed", api::BCType::VerticallyFixed)
        .value("FullyFixed", api::BCType::FullyFixed);
    nb::enum_<api::FlowBCType>(m, "FlowBCType", nb::is_arithmetic())
        .value("Closed", api::FlowBCType::Closed)
        .value("Head", api::FlowBCType::Head)
        .value("Seepage", api::FlowBCType::Seepage)
        .value("Flux", api::FlowBCType::Flux);
    nb::enum_<api::StructKind>(m, "StructKind", nb::is_arithmetic())
        .value("Plate", api::StructKind::Plate)
        .value("Anchor", api::StructKind::Anchor)
        .value("Geogrid", api::StructKind::Geogrid)
        .value("EmbeddedBeam", api::StructKind::EmbeddedBeam)
        .value("Interface", api::StructKind::Interface);
    nb::enum_<api::LoadKind>(m, "LoadKind", nb::is_arithmetic())
        .value("Point", api::LoadKind::Point)
        .value("Distributed", api::LoadKind::Distributed);
    nb::enum_<api::InitialProcedure>(m, "InitialProcedure", nb::is_arithmetic())
        .value("K0Procedure", api::InitialProcedure::K0Procedure)
        .value("GravityLoading", api::InitialProcedure::GravityLoading)
        .value("Safety", api::InitialProcedure::Safety);
    nb::enum_<api::SeismicWave>(m, "SeismicWave", nb::is_arithmetic())
        .value("Harmonic", api::SeismicWave::Harmonic)
        .value("Ricker", api::SeismicWave::Ricker)
        .value("Record", api::SeismicWave::Record);
    nb::enum_<api::DesignApproach>(m, "DesignApproach", nb::is_arithmetic())
        .value("None_", api::DesignApproach::None)
        .value("EC7_DA1_C1", api::DesignApproach::EC7_DA1_C1)
        .value("EC7_DA1_C2", api::DesignApproach::EC7_DA1_C2)
        .value("EC7_DA2", api::DesignApproach::EC7_DA2)
        .value("EC7_DA3", api::DesignApproach::EC7_DA3)
        .value("TBDY2018_Static", api::DesignApproach::TBDY2018_Static)
        .value("TBDY2018_Seismic", api::DesignApproach::TBDY2018_Seismic);
    nb::enum_<api::Severity>(m, "Severity", nb::is_arithmetic())
        .value("Warning", api::Severity::Warning)
        .value("Error", api::Severity::Error);
    nb::enum_<katai::jobs::JobState>(m, "JobState", nb::is_arithmetic())
        .value("Pending", katai::jobs::JobState::Pending)
        .value("Validating", katai::jobs::JobState::Validating)
        .value("Meshing", katai::jobs::JobState::Meshing)
        .value("Solving", katai::jobs::JobState::Solving)
        .value("Done", katai::jobs::JobState::Done)
        .value("Failed", katai::jobs::JobState::Failed)
        .value("Cancelled", katai::jobs::JobState::Cancelled);

    // ---------------------------------------------------------- schema types --
    auto material_cls = nb::class_<api::Material>(m, "Material");
    material_cls
        .def(nb::init<>())
        .def_rw("name", &api::Material::name)
        .def_rw("model", &api::Material::model)
        .def_rw("drainage", &api::Material::drainage)
        .def_rw("gamma_unsat", &api::Material::gamma_unsat, "unit weight above water [kN/m3]")
        .def_rw("gamma_sat", &api::Material::gamma_sat, "unit weight below water [kN/m3]")
        .def_rw("e_init", &api::Material::e_init)
        .def_rw("E", &api::Material::E, "Young's modulus [kN/m2]")
        .def_rw("nu", &api::Material::nu)
        .def_rw("c", &api::Material::c, "cohesion (Undrained B: su) [kN/m2]")
        .def_rw("phi", &api::Material::phi, "friction angle [deg]")
        .def_rw("psi", &api::Material::psi, "dilatancy angle [deg]")
        .def_rw("E_inc", &api::Material::E_inc)
        .def_rw("c_inc", &api::Material::c_inc)
        .def_rw("y_ref", &api::Material::y_ref)
        .def_rw("tension_cutoff", &api::Material::tension_cutoff)
        .def_rw("dilatancy_cutoff", &api::Material::dilatancy_cutoff,
                "stop dilatancy at the critical void ratio e_max (PLAXIS MMM Eq. 5.16b)")
        .def_rw("e_max", &api::Material::e_max, "critical (maximum) void ratio")
        .def_rw("tensile_strength", &api::Material::tensile_strength)
        .def_rw("E50ref", &api::Material::E50ref)
        .def_rw("Eoedref", &api::Material::Eoedref)
        .def_rw("Eurref", &api::Material::Eurref)
        .def_rw("m", &api::Material::m)
        .def_rw("nu_ur", &api::Material::nu_ur)
        .def_rw("p_ref", &api::Material::p_ref)
        .def_rw("Rf", &api::Material::Rf)
        .def_rw("k0nc_auto", &api::Material::k0nc_auto)
        .def_rw("k0nc", &api::Material::k0nc)
        .def_rw("G0ref", &api::Material::G0ref)
        .def_rw("gamma07", &api::Material::gamma07)
        .def_rw("lam_star", &api::Material::lam_star)
        .def_rw("kap_star", &api::Material::kap_star)
        .def_rw("mu_star", &api::Material::mu_star)
        .def_rw("kx", &api::Material::kx, "permeability x [m/day]")
        .def_rw("ky", &api::Material::ky, "permeability y [m/day]")
        .def_rw("gw_ga", &api::Material::gw_ga)
        .def_rw("gw_gn", &api::Material::gw_gn)
        .def_rw("gw_gl", &api::Material::gw_gl)
        .def_rw("gw_Sres", &api::Material::gw_Sres)
        .def_rw("rinter_rigid", &api::Material::rinter_rigid)
        .def_rw("Rinter", &api::Material::Rinter)
        .def_rw("k0_auto", &api::Material::k0_auto)
        .def_rw("k0", &api::Material::k0)
        .def_rw("oc_mode", &api::Material::oc_mode)
        .def_rw("OCR", &api::Material::OCR)
        .def_rw("POP", &api::Material::POP);
    bind_color(material_cls);

    auto plate_cls = nb::class_<api::PlateMaterial>(m, "PlateMaterial");
    plate_cls
        .def(nb::init<>())
        .def_rw("name", &api::PlateMaterial::name)
        .def_rw("EA", &api::PlateMaterial::EA, "[kN/m]")
        .def_rw("EI", &api::PlateMaterial::EI, "[kN m2/m]")
        .def_rw("w", &api::PlateMaterial::w, "weight [kN/m/m]")
        .def_rw("nu", &api::PlateMaterial::nu)
        .def_rw("elastoplastic", &api::PlateMaterial::elastoplastic)
        .def_rw("Mp", &api::PlateMaterial::Mp)
        .def_rw("Np", &api::PlateMaterial::Np);
    bind_color(plate_cls);

    auto anchor_cls = nb::class_<api::AnchorMaterial>(m, "AnchorMaterial");
    anchor_cls
        .def(nb::init<>())
        .def_rw("name", &api::AnchorMaterial::name)
        .def_rw("EA", &api::AnchorMaterial::EA, "[kN]")
        .def_rw("Lspacing", &api::AnchorMaterial::Lspacing, "out-of-plane spacing [m]")
        .def_rw("elastoplastic", &api::AnchorMaterial::elastoplastic)
        .def_rw("Fmax_tens", &api::AnchorMaterial::Fmax_tens)
        .def_rw("Fmax_comp", &api::AnchorMaterial::Fmax_comp)
        .def_rw("prestress", &api::AnchorMaterial::prestress,
                "lock-off force of one anchor [kN], tension-positive; 0 = installed slack");
    bind_color(anchor_cls);

    auto geogrid_cls = nb::class_<api::GeogridMaterial>(m, "GeogridMaterial");
    geogrid_cls
        .def(nb::init<>())
        .def_rw("name", &api::GeogridMaterial::name)
        .def_rw("EA", &api::GeogridMaterial::EA, "[kN/m]")
        .def_rw("elastoplastic", &api::GeogridMaterial::elastoplastic)
        .def_rw("Np", &api::GeogridMaterial::Np);
    bind_color(geogrid_cls);

    auto embedded_cls = nb::class_<api::EmbeddedBeamMaterial>(m, "EmbeddedBeamMaterial");
    embedded_cls
        .def(nb::init<>())
        .def_rw("name", &api::EmbeddedBeamMaterial::name)
        .def_rw("E", &api::EmbeddedBeamMaterial::E)
        .def_rw("gamma", &api::EmbeddedBeamMaterial::gamma)
        .def_rw("diameter", &api::EmbeddedBeamMaterial::diameter)
        .def_rw("Lspacing", &api::EmbeddedBeamMaterial::Lspacing)
        .def_rw("Tskin_max", &api::EmbeddedBeamMaterial::Tskin_max)
        .def_rw("Fmax_base", &api::EmbeddedBeamMaterial::Fmax_base);
    bind_color(embedded_cls);

    nb::class_<api::SoilPolygon>(m, "SoilPolygon")
        .def(nb::init<>())
        .def_rw("name", &api::SoilPolygon::name)
        .def_rw("material", &api::SoilPolygon::material)
        .def_rw("x", &api::SoilPolygon::x)
        .def_rw("y", &api::SoilPolygon::y)
        .def_rw("edge_bc", &api::SoilPolygon::edge_bc)
        .def_rw("edge_flow", &api::SoilPolygon::edge_flow)
        .def_rw("edge_head", &api::SoilPolygon::edge_head)
        .def_rw("edge_flux", &api::SoilPolygon::edge_flux,
                "prescribed boundary flux per edge [m/day], inflow positive (FlowBCType.Flux)")
        .def_rw("coarseness", &api::SoilPolygon::coarseness);

    nb::class_<api::StructElement>(m, "StructElement")
        .def(nb::init<>())
        .def_rw("kind", &api::StructElement::kind)
        .def_rw("name", &api::StructElement::name)
        .def_rw("x1", &api::StructElement::x1)
        .def_rw("y1", &api::StructElement::y1)
        .def_rw("x2", &api::StructElement::x2)
        .def_rw("y2", &api::StructElement::y2)
        .def_rw("material", &api::StructElement::material)
        .def_rw("iface_pos", &api::StructElement::iface_pos)
        .def_rw("iface_neg", &api::StructElement::iface_neg)
        .def_rw("iface_material", &api::StructElement::iface_material)
        .def_rw("coarseness", &api::StructElement::coarseness);

    nb::class_<api::PrescribedDisp>(m, "PrescribedDisp")
        .def(nb::init<>())
        .def_rw("name", &api::PrescribedDisp::name)
        .def_rw("x1", &api::PrescribedDisp::x1)
        .def_rw("y1", &api::PrescribedDisp::y1)
        .def_rw("x2", &api::PrescribedDisp::x2)
        .def_rw("y2", &api::PrescribedDisp::y2)
        .def_rw("set_ux", &api::PrescribedDisp::set_ux, "the horizontal component is prescribed")
        .def_rw("ux", &api::PrescribedDisp::ux, "horizontal displacement [m] (when set_ux)")
        .def_rw("set_uy", &api::PrescribedDisp::set_uy, "the vertical component is prescribed")
        .def_rw("uy", &api::PrescribedDisp::uy, "vertical displacement [m] (when set_uy)")
        .def_rw("coarseness", &api::PrescribedDisp::coarseness);

    nb::class_<api::Load>(m, "Load")
        .def(nb::init<>())
        .def_rw("kind", &api::Load::kind)
        .def_rw("name", &api::Load::name)
        .def_rw("x1", &api::Load::x1)
        .def_rw("y1", &api::Load::y1)
        .def_rw("x2", &api::Load::x2)
        .def_rw("y2", &api::Load::y2)
        .def_rw("qx1", &api::Load::qx1, "[kN/m] (point: kN/m out-of-plane)")
        .def_rw("qy1", &api::Load::qy1)
        .def_rw("qx2", &api::Load::qx2)
        .def_rw("qy2", &api::Load::qy2)
        .def_rw("coarseness", &api::Load::coarseness);

    nb::class_<api::Phase>(m, "Phase")
        .def(nb::init<>())
        .def_rw("name", &api::Phase::name)
        .def_rw("type", &api::Phase::type)
        .def_rw("duration", &api::Phase::duration, "[day] (Dynamic: [s])")
        .def_rw("time_steps", &api::Phase::time_steps)
        .def_rw("design_approach", &api::Phase::design_approach)
        .def_rw("seismic_wave", &api::Phase::seismic_wave)
        .def_rw("seismic_amp", &api::Phase::seismic_amp)
        .def_rw("seismic_freq", &api::Phase::seismic_freq)
        .def_rw("damping_ratio", &api::Phase::damping_ratio)
        .def_rw("rayleigh_f1", &api::Phase::rayleigh_f1)
        .def_rw("rayleigh_f2", &api::Phase::rayleigh_f2)
        .def_rw("accel_record", &api::Phase::accel_record)
        .def_rw("record_dt", &api::Phase::record_dt)
        .def_rw("seismic_free_field", &api::Phase::seismic_free_field,
                "Dynamic: Lysmer free-field lateral boundaries instead of free/roller sides")
        .def_rw("seismic_compliant_base", &api::Phase::seismic_compliant_base,
                "Dynamic: absorbing (compliant) base; input becomes the upward wave, halved")
        .def_rw("dynamic_nonlinear", &api::Phase::dynamic_nonlinear,
                "Dynamic: full Newton on f_int(u) per step (plasticity during shaking); opt-in")
        .def_rw("ec8_enabled", &api::Phase::ec8_enabled)
        .def_rw("ec8_agr", &api::Phase::ec8_agr, "EC8 reference PGA on type A ground [m/s2]")
        .def_rw("ec8_gamma", &api::Phase::ec8_gamma, "EC8 importance factor gamma_I")
        .def_rw("ec8_ground", &api::Phase::ec8_ground, "EC8 ground type: 0=A .. 4=E")
        .def_rw("ec8_type", &api::Phase::ec8_type, "EC8 spectrum: 0=Type 1, 1=Type 2")
        .def_rw("tbdy_ss", &api::Phase::tbdy_ss, "TBDY short-period map acceleration S_s")
        .def_rw("tbdy_s1", &api::Phase::tbdy_s1, "TBDY 1 s map acceleration S_1")
        .def_rw("site_class", &api::Phase::site_class, "TBDY site class: 0=ZA .. 4=ZE")
        // The activity vectors are std::vector<char> in the schema (0/1 flags); bind
        // them as int lists -- nanobind would otherwise read a char as a one-character
        // string and refuse [1, 0].
        .def_prop_rw(
            "poly_active",
            [](const api::Phase& p) {
                return std::vector<int>(p.poly_active.begin(), p.poly_active.end());
            },
            [](api::Phase& p, const std::vector<int>& v) {
                p.poly_active.assign(v.begin(), v.end());
            })
        .def_prop_rw(
            "struct_active",
            [](const api::Phase& p) {
                return std::vector<int>(p.struct_active.begin(), p.struct_active.end());
            },
            [](api::Phase& p, const std::vector<int>& v) {
                p.struct_active.assign(v.begin(), v.end());
            })
        .def_prop_rw(
            "load_active",
            [](const api::Phase& p) {
                return std::vector<int>(p.load_active.begin(), p.load_active.end());
            },
            [](api::Phase& p, const std::vector<int>& v) {
                p.load_active.assign(v.begin(), v.end());
            })
        .def_prop_rw(
            "disp_active",
            [](const api::Phase& p) {
                return std::vector<int>(p.disp_active.begin(), p.disp_active.end());
            },
            [](api::Phase& p, const std::vector<int>& v) {
                p.disp_active.assign(v.begin(), v.end());
            })
        // Water conditions for this phase (staged dewatering): the phreatic polyline that
        // replaces the project's while this phase runs.
        .def_rw("water_override", &api::Phase::water_override,
                "use this phase's own phreatic line instead of the project's")
        .def_rw("wx", &api::Phase::wx, "phase phreatic polyline, x [m]")
        .def_rw("wy", &api::Phase::wy, "phase phreatic polyline, y [m]");

    nb::class_<api::MeshSettings>(m, "MeshSettings")
        .def(nb::init<>())
        .def_rw("elem_size", &api::MeshSettings::elem_size, "target element edge [m]")
        .def_rw("order", &api::MeshSettings::order, "6 = tri6, 15 = tri15")
        .def_rw("auto_refine", &api::MeshSettings::auto_refine);

    nb::class_<api::Project>(m, "Project")
        .def(nb::init<>())
        .def_rw("name", &api::Project::name)
        .def_rw("axisymmetric", &api::Project::axisymmetric)
        .def_rw("initial_procedure", &api::Project::initial_procedure)
        .def_rw("mesh", &api::Project::mesh)
        .def_rw("x_min", &api::Project::x_min)
        .def_rw("x_max", &api::Project::x_max)
        .def_rw("y_min", &api::Project::y_min)
        .def_rw("y_max", &api::Project::y_max)
        .def_rw("has_water", &api::Project::has_water)
        .def_rw("wx", &api::Project::wx)
        .def_rw("wy", &api::Project::wy)
        .def_rw("materials", &api::Project::materials)
        .def_rw("plates", &api::Project::plates)
        .def_rw("anchors", &api::Project::anchors)
        .def_rw("geogrids", &api::Project::geogrids)
        .def_rw("embedded", &api::Project::embedded)
        .def_rw("polygons", &api::Project::polygons)
        .def_rw("structs", &api::Project::structs)
        .def_rw("loads", &api::Project::loads)
        .def_rw("disps", &api::Project::disps)
        .def_rw("initial", &api::Project::initial)
        .def_rw("phases", &api::Project::phases);

    // ------------------------------------------------------------ validation --
    nb::class_<api::Issue>(m, "Issue")
        .def_ro("severity", &api::Issue::severity)
        .def_ro("path", &api::Issue::path)
        .def_ro("message", &api::Issue::message)
        .def("__repr__", [](const api::Issue& i) {
            return std::string(i.severity == api::Severity::Error ? "error: " : "warning: ") +
                   i.path + ": " + i.message;
        });
    nb::class_<api::ValidationReport>(m, "ValidationReport")
        .def_ro("issues", &api::ValidationReport::issues)
        .def_ro("errors", &api::ValidationReport::errors)
        .def_ro("warnings", &api::ValidationReport::warnings)
        .def("ok", &api::ValidationReport::ok)
        .def("__str__", &report_text);
    m.def("validate_project", &api::validate_project, "project"_a,
          "The input contract, executable: every violation at its JSON field path.");

    // --------------------------------------------------------------- file IO --
    m.def("load_project", [](const std::string& path) {
        api::Project pr;
        std::string err;
        std::vector<api::Issue> notes;
        if (!api::load_project(path, pr, &err, &notes))
            throw nb::value_error(("cannot read " + path + ": " + err).c_str());
        return nb::make_tuple(pr, notes);
    }, "path"_a, "Load a .k2d; returns (project, reader_notes). Raises ValueError on parse failure.");
    m.def("save_project", [](const api::Project& pr, const std::string& path) {
        std::string err;
        if (!api::save_project(pr, path, &err))
            throw nb::value_error(("cannot write " + path + ": " + err).c_str());
    }, "project"_a, "path"_a);
    m.def("project_to_json", &api::project_to_json, "project"_a);
    m.def("save_results", [](const std::string& path, std::uint64_t model_hash,
                             const std::vector<api::SolveResult>& results) {
        std::string err;
        if (!api::save_results(path, model_hash, results, &err))
            throw nb::value_error(("cannot write " + path + ": " + err).c_str());
    }, "path"_a, "model_hash"_a, "results"_a,
          "Write a .res results file -- the same writer `katai solve --out` uses. "
          "model_hash = fnv1a64(project_to_json(project)), the staleness stamp every "
          "front end checks before trusting a .res against a project.");
    m.def("fnv1a64", &api::fnv1a64, "text"_a,
          "The canonical model hash: fnv1a64(project_to_json(project)) is what every "
          "front end stamps into a .res for staleness detection.");
    m.def("project_from_json", [](const std::string& text) {
        api::Project pr;
        std::string err;
        std::vector<api::Issue> notes;
        if (!api::project_from_json(text, pr, &err, &notes))
            throw nb::value_error(err.c_str());
        return nb::make_tuple(pr, notes);
    }, "text"_a);

    // ----------------------------------------------------------------- results --
    nb::class_<api::MeshResult>(m, "MeshResult")
        .def_ro("ok", &api::MeshResult::ok)
        .def_ro("message", &api::MeshResult::message)
        .def_prop_ro("node_count", [](const api::MeshResult& r) { return r.mesh.node_count; })
        .def_prop_ro("element_count", [](const api::MeshResult& r) { return r.mesh.element_count; })
        .def_prop_ro("x", [](const api::MeshResult& r) { return r.mesh.x; })
        .def_prop_ro("y", [](const api::MeshResult& r) { return r.mesh.y; });

    // What a run did that its file does not literally say. The severity is the consequence,
    // not the rarity: Warning = the answer stands but the model is not literally the drawn
    // one, Refusal = the run stopped rather than answer a model the engineer did not draw.
    nb::enum_<api::DiagnosticSeverity>(m, "DiagnosticSeverity", nb::is_arithmetic())
        .value("Note", api::DiagnosticSeverity::Note)
        .value("Warning", api::DiagnosticSeverity::Warning)
        .value("Refusal", api::DiagnosticSeverity::Refusal);
    nb::class_<api::Diagnostic>(m, "Diagnostic")
        .def_ro("severity", &api::Diagnostic::severity)
        .def_ro("code", &api::Diagnostic::code,
                "stable machine tag, e.g. 'K2D-G003' -- match on this, never on the prose")
        .def_ro("subject", &api::Diagnostic::subject, "the object as the file names it")
        .def_ro("message", &api::Diagnostic::message)
        .def("__repr__", [](const api::Diagnostic& d) {
            const char* sev = d.severity == api::DiagnosticSeverity::Refusal  ? "refusal"
                            : d.severity == api::DiagnosticSeverity::Warning  ? "warning"
                                                                              : "note";
            return "<Diagnostic " + std::string(sev) + " " + d.code +
                   (d.subject.empty() ? "" : " " + d.subject) + ": " + d.message + ">";
        });

    nb::class_<api::SolveResult>(m, "SolveResult")
        .def_ro("ok", &api::SolveResult::ok)
        .def_ro("message", &api::SolveResult::message)
        .def_ro("max_disp", &api::SolveResult::max_disp, "[m]")
        .def_ro("load_factor", &api::SolveResult::load_factor)
        .def_ro("fos", &api::SolveResult::fos, "-1 when not a Safety run")
        .def_ro("fos_is_lower_bound", &api::SolveResult::fos_lower_bound,
                "True: no mechanism up to the cap; fos is a LOWER BOUND, report '>' not '='")
        .def_prop_ro("displacement", [](const api::SolveResult& r) { return r.disp; },
                     "full DOF vector (2*node_count), [m] -- copied to numpy")
        .def_prop_ro("pore", [](const api::SolveResult& r) { return r.pore; },
                     "nodal pore pressure [kN/m2]")
        .def_prop_ro("node_x", [](const api::SolveResult& r) { return r.mesh.x; })
        .def_prop_ro("node_y", [](const api::SolveResult& r) { return r.mesh.y; })
        .def_prop_ro("stress", [](const api::SolveResult& r) {
            // nodal effective [sxx, syy, sxy] as an (n, 3) array, copied.
            const size_t n = r.stress.stress.size();
            Eigen::MatrixXd s(n, 3);
            for (size_t i = 0; i < n; ++i) s.row(i) = r.stress.stress[i].transpose();
            return s;
        }, "nodal effective stress (n, 3): sxx, syy, sxy [kN/m2] -- copied to numpy")
        .def_prop_ro("reaction", [](const api::SolveResult& r) { return r.reaction; },
                     "support reactions at fixed dofs [kN/m], (2*node_count) like displacement; "
                     "soil contribution, static phases only (empty otherwise)")
        .def_prop_ro("diagnostics", [](const api::SolveResult& r) { return r.diagnostics; },
                     "everything the run did that the file does not literally say: clipped "
                     "geometry, a fallback taken, the refusal that stopped it (list[Diagnostic])")
        .def_prop_ro("consol_time", [](const api::SolveResult& r) { return r.consol_time; })
        .def_prop_ro("consol_settlement",
                     [](const api::SolveResult& r) { return r.consol_settlement; })
        .def_prop_ro("consol_excess_pore",
                     [](const api::SolveResult& r) { return r.consol_excess_pore; });

    // -------------------------------------------------------------------- job --
    nb::class_<katai::jobs::PhaseTiming>(m, "PhaseTiming")
        .def_ro("name", &katai::jobs::PhaseTiming::name)
        .def_ro("seconds", &katai::jobs::PhaseTiming::seconds);

    nb::class_<api::Job>(m, "Job")
        .def(nb::init<api::Project>(), "project"_a,
             "One headless run: validate -> mesh from the project's own settings -> "
             "every phase through the driver. run() blocks; the progress callback "
             "fires on the calling thread.")
        .def("set_on_phase", &api::Job::set_on_phase, "callback"_a,
             "callback(index, total, name), called before each phase solve")
        // run() releases the progress callback when it returns: the Job is one-shot,
        // the callback belongs to the run -- and a C++-held reference to a Python
        // callable is an edge the garbage collector cannot see. Left in place it
        // closes an uncollectable cycle (job -> C++ -> callable -> module dict ->
        // job) and nanobind reports every instance still alive at interpreter
        // shutdown (measured: two examples printed a leak dump on exit).
        .def("run", [](api::Job& j) {
            const bool ok = j.run();
            j.set_on_phase({});
            return ok;
        })
        .def("request_cancel", &api::Job::request_cancel)
        .def("state", &api::Job::state)
        .def("message", &api::Job::message)
        .def("report", &api::Job::report, nb::rv_policy::reference_internal)
        .def("results", &api::Job::results, nb::rv_policy::reference_internal)
        .def("timings", &api::Job::timings, nb::rv_policy::reference_internal);

    // ------------------------------------------------------------- provenance --
    m.def("backend_name", &api::backend_name,
          "the linear-solver backend actually linked into this module");
    m.attr("__version__") = api::kVersion;   // the ONE identity (<katai/api/version.hpp>)
    m.attr("__version_date__") = api::kVersionDate;
    m.attr("PROJECT_FILE_VERSION") = api::kProjectFileVersion;
    m.attr("RESULTS_FILE_VERSION") = (unsigned)api::kResultsFileVersion;
}
