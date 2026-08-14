#pragma once
// Project save / load -- versioned JSON, no third-party dependency. The writer emits compact
// standard JSON ("%.17g" doubles -> IEEE round-trip exact); the reader is a small strict
// recursive-descent parser (objects / arrays / strings / numbers / bools / null). Missing keys
// keep the default-constructed value, so OLD files keep loading as the model grows (forward
// compatibility); "katai2d": <version> guards breaking changes.
//
// Round-trip equality across every model field is pinned by tests/test_project_io.cpp.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <katai/io/issue.hpp>
#include <katai/model/project.hpp>

namespace katai::model {

// v2 (2026-08): line prescribed displacements (`disps` + phases[].disp). The bump is the
// point, not a formality: an older build reading a v2 file would silently DROP the
// prescribed displacements and solve a different problem -- the version guard turns that
// silent wrong number into an honest "newer version" refusal. (The mesh/initial_procedure
// additions stayed v1 because an old build ignoring them still solves the same problem,
// just on its own defaults.)
//
// v3 (2026-08): anchor PRESTRESS (`anchors[].prestress`). Same rule, same reason: an older
// build would read the file, ignore the lock-off force, and report a wall that deflects far
// more than the one described -- an unsafe-sided silent difference, so the version refuses
// instead. Until this field existed no anchored excavation could be modelled as it is built,
// which is why the case that exposed it was an industry benchmark rather than a unit test.
//
// v4 (2026-08): PER-PHASE water conditions (`phases[].water_override` + `wx`/`wy`). The same
// benchmark needed it one line further down its specification: the pit is dewatered from -3 m
// to -17.9 m BEFORE the first excavation step, and an older build would run every phase at the
// project's level and report an excavation that was never dewatered.
// v5 (2026-08): the DILATANCY CUT-OFF (`materials[].dilatancy_cutoff`, `e_max`). An older
// build would read the file, ignore the cut-off, and let a dense sand dilate without limit --
// which over-predicts bearing capacity. Unsafe and silent, so the guard refuses the file.
//
// v6 (2026-08): the PRESCRIBED BOUNDARY FLUX (`polygons[].edge_flux`, plus `Flux` in
// `edge_flow`). An older build reads `edge_flow` as an int it does not know and would fall
// through to "closed" -- rainfall infiltration or a recharge boundary silently removed from
// the problem, and a phreatic surface reported lower than the one described. The bump makes
// that an honest refusal.
//
// v7 (2026-08): the per-phase NUMERICAL CONTROLS (`phases[].tol`, `loadsteps`, `maxiter`).
// These exist so that a published run can be reproduced by someone else, which makes dropping
// them the exact failure they were added to prevent: an older build would read the file, solve
// it at its own tolerance, and report a number that cannot be compared with the one the file
// was written to carry. The keys are written only when set, so a project that never touches
// numerics is byte-identical apart from this version digit.
//
// v8 (2026-08): the STAGED-CONSTRUCTION TARGET and the UNDRAINED SWITCH (`phases[].mstage`,
// `ignoreund`). Both change what the phase IS, and an older build would drop both silently: a
// half-applied excavation lift would be reported as a completed one (unsafe -- the wall carries
// the full stage in the file but not in the run), and a phase asked to ignore undrained
// behaviour would generate excess pore pressures the author excluded on purpose.
//
// v9 (2026-08): the PER-MATERIAL UNDRAINED STIFFNESS (`materials[].und_mode`, `nu_u`,
// `skempton_B`). An older build reads none of them and solves every undrained material at
// nu_u = 0.495 -- the value PLAXIS uses when the user says nothing, applied to a user who said
// something. A clay entered with a measured Skempton B of 0.90 would be run at 0.978, generating
// more excess pore pressure and less effective stress than the file describes, silently and in
// whichever direction the case happens to be sensitive to.
//
// v10 (2026-08): UNDRAINED (C), the total-stress drainage type (`materials[].drainage` = 4).
// An older build meets a drainage value it does not know, and what the reader does with any
// unknown enum is load it "for display as Drained". Here that would be a soil whose undrained
// parameters are read as effective ones -- a c_u of 60 kPa taken for a c' of 60 kPa, with pore
// pressures subtracted from stresses that were never meant to have any. The version refuses in
// the reader instead, which no path can skip.
//
// v11 (2026-08): WELLS AND DRAINS (`hydros`, `phases[].hydro`). An older build ignores the array
// and solves the same ground with no dewatering in it: the head stays where the boundaries put
// it, the excavation floor is reported dry when the file says it is being pumped, and every
// uplift and seepage number that follows belongs to a different problem.
//
// v12 (2026-08): CROSS PERMEABILITY of walls and interfaces (`structs[].flow_barrier`, `hyd_res`).
// An older build reads neither and computes the flow straight THROUGH a cut-off wall: the head
// difference the wall was drawn to hold does not appear, the uplift under an excavation comes out
// too low and the inflow too high. Unsafe in both, and invisible in the result.
// v13 (2026-08): the embedded beam's CONNECTION POINT (`structs[].conn`, 0 hinged / 1 free), and
// with it a corrected default. Until now the beam's top was always FREE -- coupled to the soil
// only through the skin springs -- which is not what PLAXIS does: with no structure sharing the
// point it connects the beam node HINGED to the soil node there. The consequence was not subtle
// and not visible: a point load at a pile head is carried by the nearest SOIL node, so with a
// free top the pile barely engaged (measured: a 10 m concrete pile row changed max|u| by 0.7%),
// and a pile row could not be loaded at its head at all. An older build reads no `conn` and
// solves the free version of every file, which is a different structure with the same drawing.
//
// v14 (2026-08): RESET SMALL STRAIN (`phases[].resetsmall`), PLAXIS's phase option of the same
// name (Material Models Manual sec. 7.6). The flag says the Hardening Soil small history is to
// be cleared at the start of the phase, so the soil meets it at G0. An older build reads no
// such key and runs the phase with whatever history the earlier phases accumulated -- which is
// the degraded stiffness, not the reset one. The direction is systematic: the run comes out
// SOFTER than the file asks for, with more settlement and more wall deflection, and nothing in
// the output says a stiffness reset was requested and dropped. The file's own reason for the
// reset -- typically a surcharge placed and removed to leave a preconsolidation pressure, whose
// strain history ageing would long since have erased -- is invisible to that older build, so
// the phase it runs is a different problem with the same drawing.
inline constexpr int kProjectFileVersion = 14;

// ---------------------------------------------------------------- minimal JSON value + parser --
struct Json {
    enum Type { Null, Bool, Num, Str, Arr, Obj };
    Type type = Null;
    bool b = false;
    double n = 0.0;
    std::string s;
    std::vector<Json> a;
    std::vector<std::pair<std::string, Json>> o;   // small objects: ordered, linear lookup

    const Json* find(const char* key) const {
        for (const auto& kv : o)
            if (kv.first == key) return &kv.second;
        return nullptr;
    }
    double num(const char* key, double dflt) const {
        const Json* v = find(key); return v && v->type == Num ? v->n : dflt;
    }
    bool flag(const char* key, bool dflt) const {
        const Json* v = find(key); return v && v->type == Bool ? v->b : dflt;
    }
    std::string str(const char* key, const std::string& dflt) const {
        const Json* v = find(key); return v && v->type == Str ? v->s : dflt;
    }
    // Numeric array helpers (absent key -> empty).
    std::vector<double> nums(const char* key) const {
        std::vector<double> r;
        const Json* v = find(key);
        if (v && v->type == Arr)
            for (const auto& it : v->a) r.push_back(it.type == Num ? it.n : 0.0);
        return r;
    }
    std::vector<int> ints(const char* key) const {
        std::vector<int> r;
        for (double d : nums(key)) r.push_back((int)d);
        return r;
    }
    std::vector<char> chars(const char* key) const {
        std::vector<char> r;
        for (double d : nums(key)) r.push_back(d != 0.0 ? 1 : 0);
        return r;
    }
};

namespace detail {

struct Parser {
    const char* p;
    const char* e;
    std::string* err;

    bool fail(const char* m) {
        if (err && err->empty()) *err = m;
        return false;
    }
    void ws() {
        while (p < e && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) ++p;
    }
    bool string(std::string& s) {
        ++p;   // opening quote
        while (p < e && *p != '"') {
            if (*p == '\\') {
                ++p;
                if (p >= e) return fail("bad escape");
                switch (*p) {
                    case '"': s += '"'; break;
                    case '\\': s += '\\'; break;
                    case '/': s += '/'; break;
                    case 'n': s += '\n'; break;
                    case 't': s += '\t'; break;
                    case 'r': s += '\r'; break;
                    case 'b': s += '\b'; break;
                    case 'f': s += '\f'; break;
                    case 'u': {   // we never emit \u; accept + decode the Latin-1 range
                        if (e - p < 5) return fail("bad \\u escape");
                        char hex[5] = {p[1], p[2], p[3], p[4], 0};
                        const long c = std::strtol(hex, nullptr, 16);
                        s += c < 256 ? (char)c : '?';
                        p += 4;
                        break;
                    }
                    default: return fail("unknown escape");
                }
                ++p;
            } else {
                s += *p++;
            }
        }
        if (p >= e) return fail("unterminated string");
        ++p;   // closing quote
        return true;
    }
    bool value(Json& v) {
        ws();
        if (p >= e) return fail("unexpected end of input");
        switch (*p) {
            case '{': {
                v.type = Json::Obj;
                ++p; ws();
                if (p < e && *p == '}') { ++p; return true; }
                while (true) {
                    ws();
                    if (p >= e || *p != '"') return fail("expected object key");
                    std::string key;
                    if (!string(key)) return false;
                    ws();
                    if (p >= e || *p != ':') return fail("expected ':'");
                    ++p;
                    Json item;
                    if (!value(item)) return false;
                    v.o.emplace_back(std::move(key), std::move(item));
                    ws();
                    if (p < e && *p == ',') { ++p; continue; }
                    if (p < e && *p == '}') { ++p; return true; }
                    return fail("expected ',' or '}'");
                }
            }
            case '[': {
                v.type = Json::Arr;
                ++p; ws();
                if (p < e && *p == ']') { ++p; return true; }
                while (true) {
                    Json item;
                    if (!value(item)) return false;
                    v.a.push_back(std::move(item));
                    ws();
                    if (p < e && *p == ',') { ++p; continue; }
                    if (p < e && *p == ']') { ++p; return true; }
                    return fail("expected ',' or ']'");
                }
            }
            case '"': v.type = Json::Str; return string(v.s);
            case 't':
                if (e - p >= 4 && std::strncmp(p, "true", 4) == 0) { v.type = Json::Bool; v.b = true; p += 4; return true; }
                return fail("bad token");
            case 'f':
                if (e - p >= 5 && std::strncmp(p, "false", 5) == 0) { v.type = Json::Bool; v.b = false; p += 5; return true; }
                return fail("bad token");
            case 'n':
                if (e - p >= 4 && std::strncmp(p, "null", 4) == 0) { v.type = Json::Null; p += 4; return true; }
                return fail("bad token");
            default: {
                char* q = nullptr;
                v.n = std::strtod(p, &q);
                if (q == p || q > e) return fail("bad number");
                v.type = Json::Num;
                p = q;
                return true;
            }
        }
    }
};

inline bool parse(const std::string& text, Json& out, std::string* err) {
    Parser ps{text.data(), text.data() + text.size(), err};
    if (!ps.value(out)) return false;
    ps.ws();
    if (ps.p != ps.e) { if (err && err->empty()) *err = "trailing characters"; return false; }
    return true;
}

// ------------------------------------------------------------------------------ JSON writing --
inline void wstr(std::string& o, const std::string& s) {
    o += '"';
    for (char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\t': o += "\\t"; break;
            case '\r': o += "\\r"; break;
            default:
                if ((unsigned char)c < 0x20) { char b[8]; std::snprintf(b, sizeof(b), "\\u%04x", c); o += b; }
                else o += c;
        }
    }
    o += '"';
}
inline void wnum(std::string& o, double v) {
    if (!(v == v) || v > 1e300 || v < -1e300) v = 0.0;   // NaN/inf never belong in a file
    char b[32];
    std::snprintf(b, sizeof(b), "%.17g", v);
    o += b;
}
inline void wkey(std::string& o, const char* k) { wstr(o, k); o += ':'; }
template <class T>
inline void warr(std::string& o, const char* k, const std::vector<T>& v) {
    wkey(o, k); o += '[';
    for (size_t i = 0; i < v.size(); ++i) { if (i) o += ','; wnum(o, (double)v[i]); }
    o += "],";
}
inline void wcolor(std::string& o, const float c[3]) {
    wkey(o, "color"); o += '[';
    for (int i = 0; i < 3; ++i) { if (i) o += ','; wnum(o, c[i]); }
    o += "],";
}
inline void wfield(std::string& o, const char* k, double v) { wkey(o, k); wnum(o, v); o += ','; }
inline void wfield(std::string& o, const char* k, bool v) { wkey(o, k); o += v ? "true" : "false"; o += ','; }
inline void wfield(std::string& o, const char* k, const std::string& v) { wkey(o, k); wstr(o, v); o += ','; }
inline void closeobj(std::string& o) {   // replace the trailing comma with the closing brace
    if (!o.empty() && o.back() == ',') o.back() = '}';
    else o += '}';
}

inline void rcolor(const Json& j, float c[3]) {
    const auto v = j.nums("color");
    for (size_t i = 0; i < 3 && i < v.size(); ++i) c[i] = (float)v[i];
}

inline void wphase(std::string& o, const char* key, const Phase& ph) {
    wkey(o, key); o += '{';
    wfield(o, "name", ph.name);
    wfield(o, "type", (double)(int)ph.type);
    wfield(o, "duration", ph.duration);
    wfield(o, "steps", (double)ph.time_steps);
    wfield(o, "design", (double)(int)ph.design_approach);
    wfield(o, "seiswave", (double)(int)ph.seismic_wave);
    wfield(o, "seisamp", ph.seismic_amp);
    wfield(o, "seisfreq", ph.seismic_freq);
    wfield(o, "damp", ph.damping_ratio);
    wfield(o, "rayf1", ph.rayleigh_f1);
    wfield(o, "rayf2", ph.rayleigh_f2);
    wfield(o, "seisff", ph.seismic_free_field ? 1.0 : 0.0);
    wfield(o, "dynnl", ph.dynamic_nonlinear ? 1.0 : 0.0);
    wfield(o, "seiscb", ph.seismic_compliant_base ? 1.0 : 0.0);
    wfield(o, "ec8on", ph.ec8_enabled ? 1.0 : 0.0);
    wfield(o, "ec8agr", ph.ec8_agr);
    wfield(o, "ec8gi", ph.ec8_gamma);
    wfield(o, "ec8gnd", (double)ph.ec8_ground);
    wfield(o, "ec8typ", (double)ph.ec8_type);
    if (!ph.accel_record.empty()) {
        warr(o, "rec", ph.accel_record);
        wfield(o, "recdt", ph.record_dt);
    }
    wfield(o, "tbdyss", ph.tbdy_ss);
    wfield(o, "tbdys1", ph.tbdy_s1);
    wfield(o, "siteclass", (double)ph.site_class);
    warr(o, "poly", ph.poly_active);
    warr(o, "struct", ph.struct_active);
    warr(o, "load", ph.load_active);
    warr(o, "disp", ph.disp_active);
    warr(o, "hydro", ph.hydro_active);
    wfield(o, "water_override", ph.water_override);
    warr(o, "wx", ph.wx); warr(o, "wy", ph.wy);
    // Numerical controls, written only when SET (like the accelerogram above). Zero is "let the
    // program choose", and a file full of zeros would say that in more bytes while making every
    // project that never touches numerics differ from the one it was yesterday.
    if (ph.tolerance > 0.0) wfield(o, "tol", ph.tolerance);
    if (ph.load_steps > 0) wfield(o, "loadsteps", (double)ph.load_steps);
    if (ph.max_iterations > 0) wfield(o, "maxiter", (double)ph.max_iterations);
    // Same rule for the staged-construction target and the undrained switch: the default IS the
    // ordinary case (the whole stage, undrained soil behaving undrained), so it costs no bytes.
    if (ph.sum_mstage != 1.0) wfield(o, "mstage", ph.sum_mstage);
    if (ph.ignore_undrained) wfield(o, "ignoreund", true);
    if (ph.reset_small_strain) wfield(o, "resetsmall", true);
    closeobj(o); o += ',';
}
inline Phase rphase(const Json& j) {
    Phase ph;
    ph.name = j.str("name", ph.name);
    ph.type = (PhaseType)(int)j.num("type", 0);
    ph.duration = j.num("duration", ph.duration);
    ph.time_steps = (int)j.num("steps", ph.time_steps);
    ph.design_approach = (DesignApproach)(int)j.num("design", 0);  // 0 = None (backward compatible)
    ph.seismic_wave = (SeismicWave)(int)j.num("seiswave", 0);
    ph.seismic_amp = j.num("seisamp", ph.seismic_amp);
    ph.seismic_freq = j.num("seisfreq", ph.seismic_freq);
    ph.damping_ratio = j.num("damp", ph.damping_ratio);
    ph.rayleigh_f1 = j.num("rayf1", ph.rayleigh_f1);
    ph.rayleigh_f2 = j.num("rayf2", ph.rayleigh_f2);
    ph.seismic_free_field = j.num("seisff", ph.seismic_free_field ? 1.0 : 0.0) != 0.0;
    ph.dynamic_nonlinear = j.num("dynnl", ph.dynamic_nonlinear ? 1.0 : 0.0) != 0.0;
    ph.seismic_compliant_base = j.num("seiscb", ph.seismic_compliant_base ? 1.0 : 0.0) != 0.0;
    ph.ec8_enabled = j.num("ec8on", ph.ec8_enabled ? 1.0 : 0.0) != 0.0;
    ph.ec8_agr = j.num("ec8agr", ph.ec8_agr);
    ph.ec8_gamma = j.num("ec8gi", ph.ec8_gamma);
    ph.ec8_ground = (int)j.num("ec8gnd", (double)ph.ec8_ground);
    ph.ec8_type = (int)j.num("ec8typ", (double)ph.ec8_type);
    ph.accel_record = j.nums("rec");
    ph.record_dt = j.num("recdt", ph.record_dt);
    ph.tbdy_ss = j.num("tbdyss", ph.tbdy_ss);
    ph.tbdy_s1 = j.num("tbdys1", ph.tbdy_s1);
    ph.site_class = (int)j.num("siteclass", ph.site_class);
    ph.poly_active = j.chars("poly");
    ph.struct_active = j.chars("struct");
    ph.load_active = j.chars("load");
    ph.disp_active = j.chars("disp");
    ph.hydro_active = j.chars("hydro");
    ph.water_override = j.flag("water_override", ph.water_override);
    ph.wx = j.nums("wx"); ph.wy = j.nums("wy");
    ph.tolerance = j.num("tol", ph.tolerance);
    ph.load_steps = (int)j.num("loadsteps", ph.load_steps);
    ph.max_iterations = (int)j.num("maxiter", ph.max_iterations);
    ph.sum_mstage = j.num("mstage", ph.sum_mstage);
    ph.ignore_undrained = j.flag("ignoreund", ph.ignore_undrained);
    ph.reset_small_strain = j.flag("resetsmall", ph.reset_small_strain);
    return ph;
}

}  // namespace detail

// ------------------------------------------------------------------------- Project -> JSON --
inline std::string project_to_json(const Project& p) {
    using namespace detail;
    std::string o;
    o.reserve(4096);
    o += '{';
    wfield(o, "katai2d", (double)kProjectFileVersion);
    wfield(o, "name", p.name);
    wfield(o, "axisymmetric", p.axisymmetric);
    wfield(o, "initial_procedure", (double)(int)p.initial_procedure);
    wfield(o, "x_min", p.x_min); wfield(o, "x_max", p.x_max);
    wfield(o, "y_min", p.y_min); wfield(o, "y_max", p.y_max);
    wfield(o, "has_water", p.has_water);
    warr(o, "wx", p.wx);
    warr(o, "wy", p.wy);
    wkey(o, "mesh"); o += '{';
    wfield(o, "elem_size", p.mesh.elem_size);
    wfield(o, "order", (double)p.mesh.order);
    wfield(o, "auto_refine", p.mesh.auto_refine);
    closeobj(o); o += ',';

    wkey(o, "materials"); o += '[';
    for (size_t i = 0; i < p.materials.size(); ++i) {
        const Material& m = p.materials[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", m.name);
        wfield(o, "model", (double)(int)m.model);
        wfield(o, "drainage", (double)(int)m.drainage);
        wcolor(o, m.color);
        wfield(o, "gamma_unsat", m.gamma_unsat); wfield(o, "gamma_sat", m.gamma_sat);
        wfield(o, "e_init", m.e_init);
        wfield(o, "E", m.E); wfield(o, "nu", m.nu);
        wfield(o, "c", m.c); wfield(o, "phi", m.phi); wfield(o, "psi", m.psi);
        wfield(o, "E_inc", m.E_inc); wfield(o, "c_inc", m.c_inc); wfield(o, "y_ref", m.y_ref);
        wfield(o, "tension_cutoff", m.tension_cutoff);
    wfield(o, "dilatancy_cutoff", m.dilatancy_cutoff);
    wfield(o, "e_max", m.e_max);
        wfield(o, "tensile_strength", m.tensile_strength);
        wfield(o, "E50ref", m.E50ref); wfield(o, "Eoedref", m.Eoedref); wfield(o, "Eurref", m.Eurref);
        wfield(o, "m", m.m); wfield(o, "nu_ur", m.nu_ur); wfield(o, "p_ref", m.p_ref);
        wfield(o, "Rf", m.Rf);
        wfield(o, "k0nc_auto", m.k0nc_auto); wfield(o, "k0nc", m.k0nc);
        wfield(o, "G0ref", m.G0ref); wfield(o, "gamma07", m.gamma07);
        wfield(o, "lamstar", m.lam_star); wfield(o, "kapstar", m.kap_star);
        wfield(o, "mustar", m.mu_star);
        wfield(o, "kx", m.kx); wfield(o, "ky", m.ky);
        wfield(o, "und_mode", (double)m.und_mode);
        wfield(o, "nu_u", m.nu_u); wfield(o, "skempton_B", m.skempton_B);
        wfield(o, "gw_ga", m.gw_ga); wfield(o, "gw_gn", m.gw_gn);
        wfield(o, "gw_gl", m.gw_gl); wfield(o, "gw_Sres", m.gw_Sres);
        wfield(o, "rinter_rigid", m.rinter_rigid); wfield(o, "Rinter", m.Rinter);
        wfield(o, "k0_auto", m.k0_auto); wfield(o, "k0", m.k0);
        wfield(o, "oc_mode", (double)m.oc_mode); wfield(o, "OCR", m.OCR); wfield(o, "POP", m.POP);
        closeobj(o);
    }
    o += "],";

    wkey(o, "plates"); o += '[';
    for (size_t i = 0; i < p.plates.size(); ++i) {
        const PlateMaterial& m = p.plates[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", m.name); wcolor(o, m.color);
        wfield(o, "elastoplastic", m.elastoplastic);
        wfield(o, "EA", m.EA); wfield(o, "EI", m.EI); wfield(o, "w", m.w);
        wfield(o, "nu", m.nu); wfield(o, "Mp", m.Mp); wfield(o, "Np", m.Np);
        closeobj(o);
    }
    o += "],";
    wkey(o, "anchors"); o += '[';
    for (size_t i = 0; i < p.anchors.size(); ++i) {
        const AnchorMaterial& m = p.anchors[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", m.name); wcolor(o, m.color);
        wfield(o, "elastoplastic", m.elastoplastic);
        wfield(o, "EA", m.EA); wfield(o, "Fmax_tens", m.Fmax_tens);
        wfield(o, "Fmax_comp", m.Fmax_comp); wfield(o, "Lspacing", m.Lspacing);
        wfield(o, "prestress", m.prestress);
        closeobj(o);
    }
    o += "],";
    wkey(o, "geogrids"); o += '[';
    for (size_t i = 0; i < p.geogrids.size(); ++i) {
        const GeogridMaterial& m = p.geogrids[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", m.name); wcolor(o, m.color);
        wfield(o, "elastoplastic", m.elastoplastic);
        wfield(o, "EA", m.EA); wfield(o, "Np", m.Np);
        closeobj(o);
    }
    o += "],";
    wkey(o, "embedded"); o += '[';
    for (size_t i = 0; i < p.embedded.size(); ++i) {
        const EmbeddedBeamMaterial& m = p.embedded[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", m.name); wcolor(o, m.color);
        wfield(o, "E", m.E); wfield(o, "gamma", m.gamma);
        wfield(o, "diameter", m.diameter); wfield(o, "Lspacing", m.Lspacing);
        wfield(o, "Tskin_max", m.Tskin_max); wfield(o, "Fmax_base", m.Fmax_base);
        closeobj(o);
    }
    o += "],";

    wkey(o, "polygons"); o += '[';
    for (size_t i = 0; i < p.polygons.size(); ++i) {
        const SoilPolygon& P = p.polygons[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", P.name);
        wfield(o, "material", (double)P.material);
        wfield(o, "coarseness", P.coarseness);
        warr(o, "x", P.x);
        warr(o, "y", P.y);
        warr(o, "edge_bc", P.edge_bc);
        warr(o, "edge_flow", P.edge_flow);
        warr(o, "edge_head", P.edge_head);
        warr(o, "edge_flux", P.edge_flux);
        closeobj(o);
    }
    o += "],";

    wkey(o, "structs"); o += '[';
    for (size_t i = 0; i < p.structs.size(); ++i) {
        const StructElement& s = p.structs[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "kind", (double)(int)s.kind);
        wfield(o, "name", s.name);
        wfield(o, "x1", s.x1); wfield(o, "y1", s.y1);
        wfield(o, "x2", s.x2); wfield(o, "y2", s.y2);
        wfield(o, "material", (double)s.material);
        wfield(o, "iface_pos", s.iface_pos); wfield(o, "iface_neg", s.iface_neg);
        wfield(o, "iface_material", (double)s.iface_material);
        wfield(o, "conn", (double)s.conn);
        wfield(o, "flow_barrier", (double)s.flow_barrier);
        wfield(o, "hyd_res", s.hydraulic_resistance);
        wfield(o, "coarseness", s.coarseness);
        closeobj(o);
    }
    o += "],";

    wkey(o, "loads"); o += '[';
    for (size_t i = 0; i < p.loads.size(); ++i) {
        const Load& L = p.loads[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "kind", (double)(int)L.kind);
        wfield(o, "name", L.name);
        wfield(o, "x1", L.x1); wfield(o, "y1", L.y1);
        wfield(o, "x2", L.x2); wfield(o, "y2", L.y2);
        wfield(o, "qx1", L.qx1); wfield(o, "qy1", L.qy1);
        wfield(o, "qx2", L.qx2); wfield(o, "qy2", L.qy2);
        wfield(o, "coarseness", L.coarseness);
        closeobj(o);
    }
    o += "],";

    wkey(o, "disps"); o += '[';
    for (size_t i = 0; i < p.disps.size(); ++i) {
        const PrescribedDisp& D = p.disps[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", D.name);
        wfield(o, "x1", D.x1); wfield(o, "y1", D.y1);
        wfield(o, "x2", D.x2); wfield(o, "y2", D.y2);
        wfield(o, "set_ux", D.set_ux); wfield(o, "ux", D.ux);
        wfield(o, "set_uy", D.set_uy); wfield(o, "uy", D.uy);
        wfield(o, "coarseness", D.coarseness);
        closeobj(o);
    }
    o += "],";

    wkey(o, "hydros"); o += '[';
    for (size_t i = 0; i < p.hydros.size(); ++i) {
        const HydroLine& H = p.hydros[i];
        if (i) o += ',';
        o += '{';
        wfield(o, "name", H.name);
        wfield(o, "kind", (double)(int)H.kind);
        wfield(o, "behaviour", (double)H.behaviour);
        wfield(o, "x1", H.x1); wfield(o, "y1", H.y1);
        wfield(o, "x2", H.x2); wfield(o, "y2", H.y2);
        wfield(o, "q", H.q); wfield(o, "h_min", H.h_min); wfield(o, "head", H.head);
        wfield(o, "coarseness", H.coarseness);
        closeobj(o);
    }
    o += "],";

    detail::wphase(o, "initial", p.initial);
    wkey(o, "phases"); o += '[';
    for (size_t i = 0; i < p.phases.size(); ++i) {
        if (i) o += ',';
        std::string ph;
        detail::wphase(ph, "x", p.phases[i]);   // reuse the writer, strip the key
        o += ph.substr(4 /* "x": */, ph.size() - 5 /* trailing ',' */);
    }
    o += "],";
    closeobj(o);
    return o;
}

// ------------------------------------------------------------------------- JSON -> Project --
// `notes`, when given, collects what the reader had to decide on the caller's
// behalf: a forward-version file may carry enum values this build does not
// know, and they are CLAMPED to a safe display value so the GUI can still show
// the file. Clamping loses the original number, so it is reported here rather
// than left silent -- a front end that intends to SOLVE must treat an Error
// note exactly like a validation error.
inline bool project_from_json(const std::string& text, Project& out, std::string* err = nullptr,
                              std::vector<katai::io::Issue>* notes = nullptr) {
    using namespace detail;
    Json root;
    if (!parse(text, root, err)) return false;
    if (root.type != Json::Obj) { if (err) *err = "not a JSON object"; return false; }
    const int ver = (int)root.num("katai2d", 0);
    if (ver < 1 || ver > kProjectFileVersion) {
        if (err) *err = "not a KATAI 2D project file (or a newer version)";
        return false;
    }

    Project p;   // defaults everywhere; missing keys keep them
    p.materials.clear();
    p.name = root.str("name", p.name);
    p.axisymmetric = root.flag("axisymmetric", false);
    // Raw cast on purpose (like phases[].type): an out-of-range value is kept so the
    // validator can report it at its field path instead of a silent clamp.
    p.initial_procedure = (InitialProcedure)(int)root.num("initial_procedure", 0);
    p.x_min = root.num("x_min", p.x_min); p.x_max = root.num("x_max", p.x_max);
    p.y_min = root.num("y_min", p.y_min); p.y_max = root.num("y_max", p.y_max);
    p.has_water = root.flag("has_water", false);
    p.wx = root.nums("wx");
    p.wy = root.nums("wy");
    if (const Json* mj = root.find("mesh"); mj && mj->type == Json::Obj) {
        const Json& mo = *mj;
        p.mesh.elem_size = mo.num("elem_size", p.mesh.elem_size);
        p.mesh.order = (int)mo.num("order", (double)p.mesh.order);
        p.mesh.auto_refine = mo.flag("auto_refine", p.mesh.auto_refine);
    }

    if (const Json* arr = root.find("materials"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            Material m;
            m.name = j.str("name", m.name);
            // An out-of-range enum (a forward-version file) must not silently land on the
            // wrong model: clamp and fall back to LinearElastic (the display surfaces are
            // bounds-checked through soil_model_name as well). Clamping loses the original
            // value, so it is REPORTED through `notes` -- display may proceed, a solve must not.
            const int smv = (int)j.num("model", 0);
            m.model = (smv >= 0 && smv < kSoilModelCount) ? (SoilModel)smv : SoilModel::LinearElastic;
            if (notes && !(smv >= 0 && smv < kSoilModelCount))
                notes->push_back({katai::io::Severity::Error,
                                  "materials[" + std::to_string(p.materials.size()) + "].model",
                                  "\"" + m.name + "\": unknown soil model value " +
                                      std::to_string(smv) + " (this build knows 0.." +
                                      std::to_string(kSoilModelCount - 1) +
                                      "); loaded for display as Linear elastic -- do not solve "
                                      "with this file"});
            const int drv = (int)j.num("drainage", 0);
            const bool drv_ok = drv >= 0 && drv < kDrainageCount;
            m.drainage = drv_ok ? (Drainage)drv : Drainage::Drained;
            if (notes && !drv_ok)
                notes->push_back({katai::io::Severity::Error,
                                  "materials[" + std::to_string(p.materials.size()) + "].drainage",
                                  "\"" + m.name + "\": unknown drainage value " +
                                      std::to_string(drv) + " (this build knows 0.." +
                                      std::to_string(kDrainageCount - 1) +
                                      "); loaded for display as Drained "
                                      "-- do not solve with this file"});
            rcolor(j, m.color);
            m.gamma_unsat = j.num("gamma_unsat", m.gamma_unsat);
            m.gamma_sat = j.num("gamma_sat", m.gamma_sat);
            m.e_init = j.num("e_init", m.e_init);
            m.E = j.num("E", m.E); m.nu = j.num("nu", m.nu);
            m.c = j.num("c", m.c); m.phi = j.num("phi", m.phi); m.psi = j.num("psi", m.psi);
            m.E_inc = j.num("E_inc", m.E_inc); m.c_inc = j.num("c_inc", m.c_inc);
            m.y_ref = j.num("y_ref", m.y_ref);
            m.tension_cutoff = j.flag("tension_cutoff", m.tension_cutoff);
    m.dilatancy_cutoff = j.flag("dilatancy_cutoff", m.dilatancy_cutoff);
    m.e_max = j.num("e_max", m.e_max);
            m.tensile_strength = j.num("tensile_strength", m.tensile_strength);
            m.E50ref = j.num("E50ref", m.E50ref); m.Eoedref = j.num("Eoedref", m.Eoedref);
            m.Eurref = j.num("Eurref", m.Eurref);
            m.m = j.num("m", m.m); m.nu_ur = j.num("nu_ur", m.nu_ur); m.p_ref = j.num("p_ref", m.p_ref);
            m.Rf = j.num("Rf", m.Rf);
            m.k0nc_auto = j.flag("k0nc_auto", m.k0nc_auto); m.k0nc = j.num("k0nc", m.k0nc);
            m.G0ref = j.num("G0ref", m.G0ref); m.gamma07 = j.num("gamma07", m.gamma07);
            m.lam_star = j.num("lamstar", m.lam_star); m.kap_star = j.num("kapstar", m.kap_star);
            m.mu_star = j.num("mustar", m.mu_star);
            m.kx = j.num("kx", m.kx); m.ky = j.num("ky", m.ky);
            // Undrained stiffness definition: 0 = nu_u direct, 1 = Skempton-B based. An
            // unrecognised value must not land on 0 in silence -- that is the reading that
            // solves at 0.495 whatever the file says, which is what the v9 bump exists to stop.
            const int undv = (int)j.num("und_mode", m.und_mode);
            m.und_mode = (undv == 0 || undv == 1) ? undv : 0;
            if (notes && !(undv == 0 || undv == 1))
                notes->push_back({katai::io::Severity::Error,
                                  "materials[" + std::to_string(p.materials.size()) + "].und_mode",
                                  "\"" + m.name + "\": unknown undrained stiffness definition " +
                                      std::to_string(undv) +
                                      " (this build knows 0 = nu_u direct, 1 = Skempton B); "
                                      "loaded for display as nu_u direct -- do not solve with "
                                      "this file"});
            m.nu_u = j.num("nu_u", m.nu_u);
            m.skempton_B = j.num("skempton_B", m.skempton_B);
            m.gw_ga = j.num("gw_ga", m.gw_ga); m.gw_gn = j.num("gw_gn", m.gw_gn);
            m.gw_gl = j.num("gw_gl", m.gw_gl); m.gw_Sres = j.num("gw_Sres", m.gw_Sres);
            m.rinter_rigid = j.flag("rinter_rigid", m.rinter_rigid);
            m.Rinter = j.num("Rinter", m.Rinter);
            m.k0_auto = j.flag("k0_auto", m.k0_auto); m.k0 = j.num("k0", m.k0);
            m.oc_mode = (int)j.num("oc_mode", m.oc_mode);
            m.OCR = j.num("OCR", m.OCR); m.POP = j.num("POP", m.POP);
            p.materials.push_back(std::move(m));
        }

    if (const Json* arr = root.find("plates"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            PlateMaterial m;
            m.name = j.str("name", m.name); rcolor(j, m.color);
            m.elastoplastic = j.flag("elastoplastic", m.elastoplastic);
            m.EA = j.num("EA", m.EA); m.EI = j.num("EI", m.EI); m.w = j.num("w", m.w);
            m.nu = j.num("nu", m.nu); m.Mp = j.num("Mp", m.Mp); m.Np = j.num("Np", m.Np);
            p.plates.push_back(std::move(m));
        }
    if (const Json* arr = root.find("anchors"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            AnchorMaterial m;
            m.name = j.str("name", m.name); rcolor(j, m.color);
            m.elastoplastic = j.flag("elastoplastic", m.elastoplastic);
            m.EA = j.num("EA", m.EA);
            m.Fmax_tens = j.num("Fmax_tens", m.Fmax_tens);
            m.Fmax_comp = j.num("Fmax_comp", m.Fmax_comp);
            m.Lspacing = j.num("Lspacing", m.Lspacing);
            m.prestress = j.num("prestress", m.prestress);
            p.anchors.push_back(std::move(m));
        }
    if (const Json* arr = root.find("geogrids"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            GeogridMaterial m;
            m.name = j.str("name", m.name); rcolor(j, m.color);
            m.elastoplastic = j.flag("elastoplastic", m.elastoplastic);
            m.EA = j.num("EA", m.EA); m.Np = j.num("Np", m.Np);
            p.geogrids.push_back(std::move(m));
        }
    if (const Json* arr = root.find("embedded"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            EmbeddedBeamMaterial m;
            m.name = j.str("name", m.name); rcolor(j, m.color);
            m.E = j.num("E", m.E); m.gamma = j.num("gamma", m.gamma);
            m.diameter = j.num("diameter", m.diameter); m.Lspacing = j.num("Lspacing", m.Lspacing);
            m.Tskin_max = j.num("Tskin_max", m.Tskin_max); m.Fmax_base = j.num("Fmax_base", m.Fmax_base);
            p.embedded.push_back(std::move(m));
        }

    if (const Json* arr = root.find("polygons"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            SoilPolygon P;
            P.name = j.str("name", P.name);
            P.material = (int)j.num("material", -1);
            P.coarseness = j.num("coarseness", 1.0);
            P.x = j.nums("x");
            P.y = j.nums("y");
            P.edge_bc = j.ints("edge_bc");
            P.edge_flow = j.ints("edge_flow");
            P.edge_head = j.nums("edge_head");
        P.edge_flux = j.nums("edge_flux");
            p.polygons.push_back(std::move(P));
        }
    if (const Json* arr = root.find("structs"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            StructElement s;
            s.kind = (StructKind)(int)j.num("kind", 0);
            s.name = j.str("name", s.name);
            s.x1 = j.num("x1", 0); s.y1 = j.num("y1", 0);
            s.x2 = j.num("x2", 0); s.y2 = j.num("y2", 0);
            s.material = (int)j.num("material", -1);
            s.iface_pos = j.flag("iface_pos", false);
            s.iface_neg = j.flag("iface_neg", false);
            s.iface_material = (int)j.num("iface_material", -1);
            // Connection point of an embedded beam. An unknown value must not land on the
            // default quietly: hinged and free are different structures, not different settings.
            const int cnv = (int)j.num("conn", 0);
            s.conn = (cnv >= 0 && cnv <= 1) ? cnv : 0;
            if (notes && !(cnv >= 0 && cnv <= 1))
                notes->push_back({katai::io::Severity::Error,
                                  "structs[" + std::to_string(p.structs.size()) + "].conn",
                                  "connection type must be 0 (hinged) or 1 (free); read " +
                                      std::to_string(cnv) + ", using hinged"});
            // Cross permeability: an unknown value must not land on "fully permeable" quietly --
            // that is the reading in which a cut-off wall stops being one.
            const int fbv = (int)j.num("flow_barrier", 0);
            s.flow_barrier = (fbv >= 0 && fbv <= 2) ? fbv : 0;
            if (notes && !(fbv >= 0 && fbv <= 2))
                notes->push_back({katai::io::Severity::Error,
                                  "structs[" + std::to_string(p.structs.size()) + "].flow_barrier",
                                  "\"" + s.name + "\": unknown cross permeability " +
                                      std::to_string(fbv) +
                                      " (this build knows 0 = fully permeable, 1 = impermeable, "
                                      "2 = semi-permeable); loaded for display as fully permeable "
                                      "-- do not solve with this file"});
            s.hydraulic_resistance = j.num("hyd_res", 0.0);
            s.coarseness = j.num("coarseness", 1.0);
            p.structs.push_back(std::move(s));
        }
    if (const Json* arr = root.find("loads"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            Load L;
            L.kind = (LoadKind)(int)j.num("kind", 0);
            L.name = j.str("name", L.name);
            L.x1 = j.num("x1", 0); L.y1 = j.num("y1", 0);
            L.x2 = j.num("x2", 0); L.y2 = j.num("y2", 0);
            L.qx1 = j.num("qx1", 0); L.qy1 = j.num("qy1", 0);
            L.qx2 = j.num("qx2", 0); L.qy2 = j.num("qy2", 0);
            L.coarseness = j.num("coarseness", 1.0);
            p.loads.push_back(std::move(L));
        }
    if (const Json* arr = root.find("disps"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            PrescribedDisp D;
            D.name = j.str("name", D.name);
            D.x1 = j.num("x1", 0); D.y1 = j.num("y1", 0);
            D.x2 = j.num("x2", 0); D.y2 = j.num("y2", 0);
            D.set_ux = j.flag("set_ux", false); D.ux = j.num("ux", 0);
            D.set_uy = j.flag("set_uy", true);  D.uy = j.num("uy", 0);
            D.coarseness = j.num("coarseness", 1.0);
            p.disps.push_back(std::move(D));
        }

    if (const Json* arr = root.find("hydros"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) {
            HydroLine H;
            H.name = j.str("name", H.name);
            // An unknown kind must not land on "Well" in silence: a drain read as a well would
            // pump a discharge of zero and quietly stop holding the head it was drawn to hold.
            const int kv = (int)j.num("kind", 0);
            const bool kind_ok = kv >= 0 && kv < kHydroKindCount;
            H.kind = kind_ok ? (HydroKind)kv : HydroKind::Well;
            if (notes && !kind_ok)
                notes->push_back({katai::io::Severity::Error,
                                  "hydros[" + std::to_string(p.hydros.size()) + "].kind",
                                  "\"" + H.name + "\": unknown hydraulic condition " +
                                      std::to_string(kv) + " (this build knows 0 = well, "
                                      "1 = drain); loaded for display as a well -- do not solve "
                                      "with this file"});
            H.behaviour = (int)j.num("behaviour", 0);
            H.x1 = j.num("x1", 0); H.y1 = j.num("y1", 0);
            H.x2 = j.num("x2", 0); H.y2 = j.num("y2", 0);
            H.q = j.num("q", 0); H.h_min = j.num("h_min", 0); H.head = j.num("head", 0);
            H.coarseness = j.num("coarseness", 1.0);
            p.hydros.push_back(std::move(H));
        }

    if (const Json* ph = root.find("initial"); ph && ph->type == Json::Obj)
        p.initial = detail::rphase(*ph);
    if (const Json* arr = root.find("phases"); arr && arr->type == Json::Arr)
        for (const Json& j : arr->a) p.phases.push_back(detail::rphase(j));

    out = std::move(p);
    return true;
}

// ----------------------------------------------------------------------------------- file IO --
inline bool save_project(const Project& p, const std::string& path, std::string* err = nullptr) {
    std::ofstream f(path, std::ios::binary);
    if (!f) { if (err) *err = "cannot open file for writing: " + path; return false; }
    const std::string text = project_to_json(p);
    f.write(text.data(), (std::streamsize)text.size());
    if (!f.good()) { if (err) *err = "write failed: " + path; return false; }
    return true;
}

inline bool load_project(const std::string& path, Project& out, std::string* err = nullptr,
                         std::vector<katai::io::Issue>* notes = nullptr) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { if (err) *err = "cannot open file: " + path; return false; }
    std::ostringstream ss;
    ss << f.rdbuf();
    return project_from_json(ss.str(), out, err, notes);
}

}  // namespace katai::model
