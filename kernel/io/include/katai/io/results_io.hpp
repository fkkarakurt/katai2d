#pragma once
// Calculation-results save / load -- compact BINARY, deliberately NOT JSON. FEM results are big
// numeric arrays (per phase: nodal displacements, stresses, pore, activity + the mesh); as JSON
// text a double costs ~20 bytes and must be re-parsed -- a 100k-node, 10-phase run would be a
// few hundred MB of text. The binary layout stores doubles raw (8 bytes, little-endian x86):
//
//   "K2DR" magic, u32 version, u64 model_hash, mesh (stored ONCE -- phases share it), u32 phases,
//   per phase: flags + scalars + message + disp/stress/pore/active + structural force diagrams.
//
// model_hash = FNV-1a of the CANONICAL project JSON (project_to_json of the loaded model), so
// results are rejected as STALE when the model no longer matches them -- even if the user
// hand-edited the project file's formatting (canonical form is formatting-independent).
// Committed Gauss states are NOT stored: restored results are for viewing/post-processing;
// continuing a staged run re-calculates (honest, and keeps the file an order of magnitude smaller).
//
// Round-trip + stale/corrupt rejection pinned in tests/test_results_io.cpp.

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <katai/analysis/results.hpp>   // SolveResult family (engine-owned since Stage B1)
#include <katai/mesh/mesh.hpp>

// The namespace is still katai::app: this file predates the io module and every
// consumer spells its API there. The rename to katai::io joins the deferred
// namespace decision (ARCHITECTURE.md carries the exception) -- moving files and
// renaming namespaces are kept as separate steps in this repository on purpose.
namespace katai::app {

// The result types are engine vocabulary (katai/analysis/results.hpp); these
// using-declarations keep this file's and its consumers' katai::app spellings
// valid without the application driver (they lawfully repeat the driver's own).
using katai::core::SolveResult;
using katai::core::Diagnostic;
using katai::core::DiagnosticSeverity;
using katai::core::StructForce;
using katai::core::InterfaceResult;
using katai::core::ForceStation;
using katai::core::InterfaceStation;

// v2 adds interface results (tau/sigma_n/slip); v3 adds StructForce/InterfaceResult::envelope (a
// Dynamic phase's forces are an envelope over the shaking, not one instant); v4 adds `superposed`
// (the stations hold the TOTAL design action = parent static state + dynamic increment) and the
// interface Coulomb demand/capacity (utilisation / max_utilisation / over_fraction); v5 adds
// InterfaceResult::slip_checked (the envelope came from the NONLINEAR Coulomb branch, so
// [SLIPPING]/[bonded] are real findings -- false = linear elastic envelope, no slip check exists).
// Older files predate each feature, so reading the fields back as false/0 is correct for them.
inline constexpr std::uint32_t kResultsFileVersion = 5;

inline std::uint64_t fnv1a64(const std::string& s) {
    std::uint64_t h = 1469598103934665603ull;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ull; }
    return h;
}

namespace detail {

struct Writer {
    std::string buf;
    void raw(const void* p, size_t n) { buf.append(static_cast<const char*>(p), n); }
    template <class T> void put(T v) { raw(&v, sizeof(T)); }
    void str(const std::string& s) { put<std::uint32_t>((std::uint32_t)s.size()); raw(s.data(), s.size()); }
    template <class T> void vec(const std::vector<T>& v) {
        put<std::uint64_t>((std::uint64_t)v.size());
        if (!v.empty()) raw(v.data(), v.size() * sizeof(T));
    }
};

struct Reader {
    const char* p;
    const char* e;
    bool ok = true;
    bool need(size_t n) { if (!ok || (size_t)(e - p) < n) ok = false; return ok; }
    template <class T> T get() {
        T v{};
        if (need(sizeof(T))) { std::memcpy(&v, p, sizeof(T)); p += sizeof(T); }
        return v;
    }
    std::string str() {
        const std::uint32_t n = get<std::uint32_t>();
        std::string s;
        if (need(n)) { s.assign(p, n); p += n; }
        return s;
    }
    template <class T> std::vector<T> vec() {
        const std::uint64_t n = get<std::uint64_t>();
        std::vector<T> v;
        if (need((size_t)n * sizeof(T))) {
            v.resize((size_t)n);
            if (n) { std::memcpy(v.data(), p, (size_t)n * sizeof(T)); p += n * sizeof(T); }
        }
        return v;
    }
};

inline void wmesh(Writer& w, const katai::mesh::Mesh& m) {
    w.put<std::int32_t>(m.node_count);
    w.put<std::int32_t>(m.element_count);
    w.put<std::int32_t>(m.nodes_per_element);
    w.vec(m.x); w.vec(m.y);
    w.vec(m.connectivity); w.vec(m.element_material);
    w.vec(m.bottom_nodes); w.vec(m.top_nodes); w.vec(m.left_nodes); w.vec(m.right_nodes);
    w.vec(m.boundary_nodes);
}
inline katai::mesh::Mesh rmesh(Reader& r) {
    katai::mesh::Mesh m;
    m.node_count = r.get<std::int32_t>();
    m.element_count = r.get<std::int32_t>();
    m.nodes_per_element = r.get<std::int32_t>();
    m.x = r.vec<double>(); m.y = r.vec<double>();
    m.connectivity = r.vec<int>(); m.element_material = r.vec<int>();
    m.bottom_nodes = r.vec<int>(); m.top_nodes = r.vec<int>();
    m.left_nodes = r.vec<int>(); m.right_nodes = r.vec<int>();
    m.boundary_nodes = r.vec<int>();
    // Internal consistency (corrupt/truncated file must not produce an out-of-range mesh).
    if ((int)m.x.size() != m.node_count || (int)m.y.size() != m.node_count ||
        (int)m.connectivity.size() != m.element_count * m.nodes_per_element ||
        (int)m.element_material.size() != m.element_count)
        r.ok = false;
    for (int c : m.connectivity)
        if (c < 0 || c >= m.node_count) { r.ok = false; break; }
    return m;
}

}  // namespace detail

inline bool save_results(const std::string& path, std::uint64_t model_hash,
                         const std::vector<SolveResult>& phases, std::string* err = nullptr) {
    if (phases.empty()) { if (err) *err = "no results to save"; return false; }
    detail::Writer w;
    w.raw("K2DR", 4);
    w.put<std::uint32_t>(kResultsFileVersion);
    w.put<std::uint64_t>(model_hash);
    detail::wmesh(w, phases.front().mesh);   // phases share the mesh; store it once
    w.put<std::uint32_t>((std::uint32_t)phases.size());
    for (const auto& R : phases) {
        w.put<std::uint8_t>(R.ok ? 1 : 0);
        w.put<std::uint8_t>(R.nil_step ? 1 : 0);
        w.put<double>(R.max_disp);
        w.put<double>(R.load_factor);
        w.put<double>(R.fos);
        w.str(R.message);
        std::vector<double> disp(R.disp.data(), R.disp.data() + R.disp.size());
        w.vec(disp);
        std::vector<double> sv;
        sv.reserve(R.stress.stress.size() * 3);
        for (const auto& s : R.stress.stress) { sv.push_back(s(0)); sv.push_back(s(1)); sv.push_back(s(2)); }
        w.vec(sv);
        w.vec(R.pore);
        w.vec(R.active);
        w.put<std::uint32_t>((std::uint32_t)R.struct_forces.size());
        for (const auto& d : R.struct_forces) {
            w.str(d.name);
            w.put<std::int32_t>(d.kind);
            w.put<std::uint8_t>(d.yielded ? 1 : 0);
            w.put<std::uint8_t>(d.envelope ? 1 : 0);   // v3: seismic envelope, not one instant
            w.put<std::uint8_t>(d.superposed ? 1 : 0); // v4: stations are the TOTAL design action
            w.put<double>(d.max_N); w.put<double>(d.max_Q); w.put<double>(d.max_M);
            w.put<std::uint32_t>((std::uint32_t)d.stations.size());
            for (const auto& st : d.stations) {
                w.put<double>(st.s); w.put<double>(st.x); w.put<double>(st.y);
                w.put<double>(st.N); w.put<double>(st.Q); w.put<double>(st.M);
            }
        }
        // Interface results (v2): tau / sigma_n / slip / gap / slipping along each Coulomb joint.
        w.put<std::uint32_t>((std::uint32_t)R.interface_forces.size());
        for (const auto& ir : R.interface_forces) {
            w.str(ir.name);
            w.put<std::uint8_t>(ir.any_slip ? 1 : 0);
            w.put<std::uint8_t>(ir.envelope ? 1 : 0);   // v3: seismic envelope (elastic branch)
            w.put<std::uint8_t>(ir.superposed ? 1 : 0); // v4: TOTAL action + a valid capacity check
            w.put<double>(ir.max_utilisation);          // v4: Coulomb demand/capacity
            w.put<double>(ir.over_fraction);
            w.put<std::uint8_t>(ir.slip_checked ? 1 : 0);  // v5: Coulomb branch -> real slip check
            w.put<double>(ir.max_abs_tau); w.put<double>(ir.max_abs_sigma_n); w.put<double>(ir.max_abs_slip);
            w.put<std::uint32_t>((std::uint32_t)ir.stations.size());
            for (const auto& st : ir.stations) {
                w.put<double>(st.s); w.put<double>(st.x); w.put<double>(st.y);
                w.put<double>(st.tau); w.put<double>(st.sigma_n);
                w.put<double>(st.slip); w.put<double>(st.gap);
                w.put<std::uint8_t>(st.slipping ? 1 : 0);
                w.put<double>(st.utilisation);          // v4
            }
        }
    }
    std::ofstream f(path, std::ios::binary);
    if (!f) { if (err) *err = "cannot open file for writing: " + path; return false; }
    f.write(w.buf.data(), (std::streamsize)w.buf.size());
    if (!f.good()) { if (err) *err = "write failed: " + path; return false; }
    return true;
}

inline bool load_results(const std::string& path, std::uint64_t model_hash,
                         std::vector<SolveResult>& out, std::string* err = nullptr) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { if (err) *err = "no results file"; return false; }
    std::ostringstream ss;
    ss << f.rdbuf();
    const std::string data = ss.str();

    detail::Reader r{data.data(), data.data() + data.size()};
    char magic[4] = {};
    if (r.need(4)) { std::memcpy(magic, r.p, 4); r.p += 4; }
    if (!r.ok || std::memcmp(magic, "K2DR", 4) != 0) { if (err) *err = "not a KATAI results file"; return false; }
    const auto ver = r.get<std::uint32_t>();
    if (ver < 1 || ver > kResultsFileVersion) { if (err) *err = "results file from a newer version"; return false; }
    const auto hash = r.get<std::uint64_t>();
    if (hash != model_hash) {
        if (err) *err = "results are stale (the model changed since they were calculated)";
        return false;
    }
    const katai::mesh::Mesh mesh = detail::rmesh(r);
    const auto nphase = r.get<std::uint32_t>();
    if (!r.ok || nphase == 0 || nphase > 4096) { if (err) *err = "corrupt results file"; return false; }

    std::vector<SolveResult> res;
    res.reserve(nphase);
    for (std::uint32_t i = 0; i < nphase && r.ok; ++i) {
        SolveResult R;
        R.ok = r.get<std::uint8_t>() != 0;
        R.nil_step = r.get<std::uint8_t>() != 0;
        R.max_disp = r.get<double>();
        R.load_factor = r.get<double>();
        R.fos = r.get<double>();
        R.message = r.str();
        const std::vector<double> disp = r.vec<double>();
        R.disp = Eigen::Map<const Eigen::VectorXd>(disp.data(), (Eigen::Index)disp.size());
        const std::vector<double> sv = r.vec<double>();
        if (sv.size() % 3 != 0) { r.ok = false; break; }
        R.stress.stress.resize(sv.size() / 3);
        for (size_t k = 0; k < R.stress.stress.size(); ++k)
            R.stress.stress[k] = Eigen::Vector3d(sv[3 * k], sv[3 * k + 1], sv[3 * k + 2]);
        R.pore = r.vec<double>();
        R.active = r.vec<char>();
        const auto nf = r.get<std::uint32_t>();
        if (nf > 100000) { r.ok = false; break; }
        for (std::uint32_t d = 0; d < nf && r.ok; ++d) {
            StructForce sf;
            sf.name = r.str();
            sf.kind = r.get<std::int32_t>();
            sf.yielded = r.get<std::uint8_t>() != 0;
            if (ver >= 3) sf.envelope = r.get<std::uint8_t>() != 0;    // absent in v1/v2 (pre-seismic)
            if (ver >= 4) sf.superposed = r.get<std::uint8_t>() != 0;  // absent in v3 (pre-superposition)
            sf.max_N = r.get<double>(); sf.max_Q = r.get<double>(); sf.max_M = r.get<double>();
            const auto ns = r.get<std::uint32_t>();
            if (ns > 10000000) { r.ok = false; break; }
            sf.stations.reserve(ns);
            for (std::uint32_t q = 0; q < ns && r.ok; ++q) {
                katai::core::ForceStation st;
                st.s = r.get<double>(); st.x = r.get<double>(); st.y = r.get<double>();
                st.N = r.get<double>(); st.Q = r.get<double>(); st.M = r.get<double>();
                sf.stations.push_back(st);
            }
            R.struct_forces.push_back(std::move(sf));
        }
        if (ver >= 2) {   // interface results (tau/sigma_n/slip); absent in v1 files
            const auto nif = r.get<std::uint32_t>();
            if (nif > 100000) { r.ok = false; break; }
            for (std::uint32_t d = 0; d < nif && r.ok; ++d) {
                InterfaceResult ir;
                ir.name = r.str();
                ir.any_slip = r.get<std::uint8_t>() != 0;
                if (ver >= 3) ir.envelope = r.get<std::uint8_t>() != 0;   // absent in v2 (pre-seismic)
                if (ver >= 4) {                                           // absent in v3
                    ir.superposed = r.get<std::uint8_t>() != 0;
                    ir.max_utilisation = r.get<double>();
                    ir.over_fraction = r.get<double>();
                }
                if (ver >= 5) ir.slip_checked = r.get<std::uint8_t>() != 0;   // absent in v4
                ir.max_abs_tau = r.get<double>(); ir.max_abs_sigma_n = r.get<double>(); ir.max_abs_slip = r.get<double>();
                const auto ns = r.get<std::uint32_t>();
                if (ns > 10000000) { r.ok = false; break; }
                ir.stations.reserve(ns);
                for (std::uint32_t q = 0; q < ns && r.ok; ++q) {
                    katai::core::InterfaceStation st;
                    st.s = r.get<double>(); st.x = r.get<double>(); st.y = r.get<double>();
                    st.tau = r.get<double>(); st.sigma_n = r.get<double>();
                    st.slip = r.get<double>(); st.gap = r.get<double>();
                    st.slipping = r.get<std::uint8_t>() != 0;
                    if (ver >= 4) st.utilisation = r.get<double>();
                    ir.stations.push_back(st);
                }
                R.interface_forces.push_back(std::move(ir));
            }
        }
        // Per-phase sanity: nodal arrays must match the stored mesh.
        if (R.disp.size() != (Eigen::Index)mesh.node_count * 2 ||
            (int)R.stress.stress.size() != mesh.node_count ||
            (int)R.pore.size() != mesh.node_count)
            r.ok = false;
        R.mesh = mesh;   // shared topology
        res.push_back(std::move(R));
    }
    if (!r.ok) { if (err) *err = "corrupt or truncated results file"; return false; }
    out = std::move(res);
    return true;
}

}  // namespace katai::app
