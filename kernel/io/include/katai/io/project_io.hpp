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
inline constexpr int kProjectFileVersion = 3;

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
        wfield(o, "tensile_strength", m.tensile_strength);
        wfield(o, "E50ref", m.E50ref); wfield(o, "Eoedref", m.Eoedref); wfield(o, "Eurref", m.Eurref);
        wfield(o, "m", m.m); wfield(o, "nu_ur", m.nu_ur); wfield(o, "p_ref", m.p_ref);
        wfield(o, "Rf", m.Rf);
        wfield(o, "k0nc_auto", m.k0nc_auto); wfield(o, "k0nc", m.k0nc);
        wfield(o, "G0ref", m.G0ref); wfield(o, "gamma07", m.gamma07);
        wfield(o, "lamstar", m.lam_star); wfield(o, "kapstar", m.kap_star);
        wfield(o, "mustar", m.mu_star);
        wfield(o, "kx", m.kx); wfield(o, "ky", m.ky);
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
            m.drainage = (drv >= 0 && drv < 4) ? (Drainage)drv : Drainage::Drained;
            if (notes && !(drv >= 0 && drv < 4))
                notes->push_back({katai::io::Severity::Error,
                                  "materials[" + std::to_string(p.materials.size()) + "].drainage",
                                  "\"" + m.name + "\": unknown drainage value " +
                                      std::to_string(drv) +
                                      " (this build knows 0..3); loaded for display as Drained "
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
