// test_project_edit -- deleting or inserting an object leaves the model that was meant.
//
// Every phase holds what it does to the project's objects BY INDEX (activation of regions,
// structures, loads, prescribed displacements, wells; dry regions, region water, region material;
// plate section and contraction). Until 2026-10-22 the Studio deleted a structure, a load, a
// prescribed displacement or a well by erasing it from its list and leaving the phases alone, so
// every object after it took over its neighbour's settings in every phase: the model solved was
// another model than the one on the screen, and nothing said so. Deleting a region shifted its
// activation but not its dry flag, its water or its phase material, and deleting a material left
// the layers, the interfaces and the phase changes that named it pointing at the wrong one.
//
// The oracle is a second construction of the project that never goes through the edit
// functions: every object is a row of a table that carries its own settings in every phase, and
// the builder writes the phase vectors from the rows. "Delete object i" must give, setting by
// setting, the project built from the table without row i; "split region i" the project built
// with row i written twice; "delete material m" the project whose rows name a material that no
// longer exists. Materials are named by NAME in the table and resolved by the builder, so the
// index arithmetic of the edit is checked against a lookup, not against itself. A negative
// control repeats the old Studio deletion and must be told apart, so the comparison can say no.

// verify: KV-GEO-004
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, phases[i].poly / struct / load / disp / hydro / dry / pwater / pwlevel / pmat / smat / scontr: a phase addresses the project's objects by index; regions, layers, interfaces and structures address the material lists by index)
//   locator:  the table-driven builder in this file, which writes every phase vector from per-object rows and resolves materials by name; it shares no code with katai/model/project.hpp's edit functions
//   quantity: every per-object phase setting (activation, dry flag, water mode and level, region material, plate section, contraction) and every material reference of the edited project, against the project built without the deleted row
//   expected: identical, setting by setting, for every deletion of every object kind, for a region split into two pieces, and for every material deletion of every list; the old Studio deletion (list erased, phases untouched) differs
//   band:     exact -- the settings are integers, flags and copied doubles

#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include <katai/model/project.hpp>
#include <katai/model/stratigraphy.hpp>

namespace m = katai::model;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

// ---- the table -------------------------------------------------------------------------------
constexpr int kPhases = 4;   // the initial phase + three
using PerPhase = std::vector<int>;   // one entry per phase, index 0 = initial

struct RegionRow {
    std::string name, material;
    PerPhase active, dry, water_mode, material_change;   // material_change: "" index into names below
    std::vector<double> water_level;
    std::vector<std::string> phase_material;   // "" = the drawn one
};
struct StructRow {
    std::string name;
    m::StructKind kind;
    std::string material;        // set name in the list of its kind
    std::string iface_material;  // soil name, "" = the soil beside the joint
    PerPhase active;
    std::vector<std::string> section;   // plates: phase section, "" = drawn
    std::vector<double> contraction;
};
struct SimpleRow {   // loads, prescribed displacements, wells
    std::string name;
    PerPhase active;
};
struct Table {
    std::vector<std::string> soils, plates, anchors, geogrids, embedded;
    std::vector<RegionRow> regions;
    std::vector<std::pair<std::string, std::string>> layers;   // name, soil
    std::vector<StructRow> structs;
    std::vector<SimpleRow> loads, disps, hydros;
};

static int find(const std::vector<std::string>& names, const std::string& n) {
    if (n.empty()) return -1;
    for (size_t i = 0; i < names.size(); ++i)
        if (names[i] == n) return (int)i;
    return -1;
}
static const std::vector<std::string>& set_names(const Table& t, m::StructKind k) {
    switch (k) {
        case m::StructKind::Plate:   return t.plates;
        case m::StructKind::Anchor:  return t.anchors;
        case m::StructKind::Geogrid: return t.geogrids;
        default:                     return t.embedded;
    }
}

// The builder: every phase vector written in full from the rows.
static m::Project build(const Table& t) {
    m::Project p;
    for (const auto& n : t.soils) { m::Material x; x.name = n; p.materials.push_back(x); }
    for (const auto& n : t.plates) { m::PlateMaterial x; x.name = n; p.plates.push_back(x); }
    for (const auto& n : t.anchors) { m::AnchorMaterial x; x.name = n; p.anchors.push_back(x); }
    for (const auto& n : t.geogrids) { m::GeogridMaterial x; x.name = n; p.geogrids.push_back(x); }
    for (const auto& n : t.embedded) { m::EmbeddedBeamMaterial x; x.name = n; p.embedded.push_back(x); }
    for (const auto& r : t.regions) {
        m::SoilPolygon P;
        P.name = r.name;
        P.material = find(t.soils, r.material);
        P.x = {0, 1, 1}; P.y = {0, 0, 1};
        p.polygons.push_back(P);
    }
    for (const auto& L : t.layers) p.strata.push_back({L.first, find(t.soils, L.second)});
    for (const auto& s : t.structs) {
        m::StructElement e;
        e.name = s.name; e.kind = s.kind;
        e.material = s.kind == m::StructKind::Interface ? -1 : find(set_names(t, s.kind), s.material);
        e.iface_material = find(t.soils, s.iface_material);
        p.structs.push_back(e);
    }
    for (const auto& l : t.loads) { m::Load x; x.name = l.name; p.loads.push_back(x); }
    for (const auto& d : t.disps) { m::PrescribedDisp x; x.name = d.name; p.disps.push_back(x); }
    for (const auto& h : t.hydros) { m::HydroLine x; x.name = h.name; p.hydros.push_back(x); }
    p.phases.resize(kPhases - 1);
    for (int k = 0; k < kPhases; ++k) {
        m::Phase& ph = k == 0 ? p.initial : p.phases[(size_t)k - 1];
        ph.name = "phase " + std::to_string(k);
        for (const auto& r : t.regions) {
            ph.poly_active.push_back((char)r.active[(size_t)k]);
            ph.poly_dry.push_back((char)r.dry[(size_t)k]);
            ph.poly_water.push_back(r.water_mode[(size_t)k]);
            ph.poly_water_level.push_back(r.water_level[(size_t)k]);
            ph.poly_material.push_back(find(t.soils, r.phase_material[(size_t)k]));
        }
        for (const auto& s : t.structs) {
            ph.struct_active.push_back((char)s.active[(size_t)k]);
            ph.struct_material.push_back(s.kind == m::StructKind::Plate ? find(t.plates, s.section[(size_t)k]) : -1);
            ph.struct_contraction.push_back(s.contraction[(size_t)k]);
        }
        for (const auto& l : t.loads) ph.load_active.push_back((char)l.active[(size_t)k]);
        for (const auto& d : t.disps) ph.disp_active.push_back((char)d.active[(size_t)k]);
        for (const auto& h : t.hydros) ph.hydro_active.push_back((char)h.active[(size_t)k]);
    }
    return p;
}

// A deterministic pattern per object and phase, different for neighbouring objects, so a shift
// by one index cannot go unnoticed.
static int pat(int obj, int k, int mod) { return (obj * 7 + k * 3 + obj * k) % mod; }

static Table make_table() {
    Table t;
    t.soils = {"Sand", "Clay", "Peat", "Gravel"};
    t.plates = {"Wall A", "Wall B", "Lining"};
    t.anchors = {"Strand", "Bar"};
    t.geogrids = {"Grid"};
    t.embedded = {"Pile"};
    for (int i = 0; i < 5; ++i) {
        RegionRow r;
        r.name = "R" + std::to_string(i);
        r.material = t.soils[(size_t)i % t.soils.size()];
        for (int k = 0; k < kPhases; ++k) {
            r.active.push_back(pat(i, k, 3) != 0);
            r.dry.push_back(pat(i, k, 2));
            r.water_mode.push_back(pat(i, k, 3));
            r.water_level.push_back(-1.0 - 0.5 * pat(i, k, 5));
            r.phase_material.push_back(pat(i, k, 4) == 1 ? t.soils[(size_t)(i + k) % t.soils.size()] : "");
        }
        t.regions.push_back(r);
    }
    t.layers = {{"Top", "Sand"}, {"Middle", "Peat"}, {"Bottom", "Clay"}};
    const m::StructKind kinds[] = {m::StructKind::Plate, m::StructKind::Anchor, m::StructKind::Plate,
                                   m::StructKind::Geogrid, m::StructKind::Interface, m::StructKind::EmbeddedBeam,
                                   m::StructKind::Plate};
    for (int i = 0; i < 7; ++i) {
        StructRow s;
        s.name = "S" + std::to_string(i);
        s.kind = kinds[i];
        const auto& names = set_names(t, s.kind);
        s.material = s.kind == m::StructKind::Interface ? "" : names[(size_t)i % names.size()];
        s.iface_material = (i % 3 == 0) ? t.soils[(size_t)i % t.soils.size()] : "";
        for (int k = 0; k < kPhases; ++k) {
            s.active.push_back(pat(i, k, 3) != 1);
            s.section.push_back(s.kind == m::StructKind::Plate && pat(i, k, 3) == 2 ? t.plates[(size_t)(i + k) % t.plates.size()] : "");
            s.contraction.push_back(s.kind == m::StructKind::Plate ? 0.001 * pat(i, k, 4) : 0.0);
        }
        t.structs.push_back(s);
    }
    auto simple = [](const std::string& stem, int n, int salt) {
        std::vector<SimpleRow> v;
        for (int i = 0; i < n; ++i) {
            SimpleRow r;
            r.name = stem + std::to_string(i);
            for (int k = 0; k < kPhases; ++k) r.active.push_back(pat(i + salt, k, 3) != 0);
            v.push_back(r);
        }
        return v;
    };
    t.loads = simple("L", 4, 1);
    t.disps = simple("D", 3, 2);
    t.hydros = simple("W", 3, 3);
    return t;
}

// ---- the comparison: the same model, setting by setting ----------------------------------------
static bool same_model(const m::Project& a, const m::Project& b, std::string& why) {
    auto fail = [&why](const std::string& w) { why = w; return false; };
    if (a.materials.size() != b.materials.size() || a.plates.size() != b.plates.size() ||
        a.anchors.size() != b.anchors.size() || a.geogrids.size() != b.geogrids.size() ||
        a.embedded.size() != b.embedded.size())
        return fail("material list sizes");
    for (size_t i = 0; i < a.materials.size(); ++i)
        if (a.materials[i].name != b.materials[i].name) return fail("soil " + std::to_string(i));
    for (size_t i = 0; i < a.plates.size(); ++i)
        if (a.plates[i].name != b.plates[i].name) return fail("plate set " + std::to_string(i));
    if (a.polygons.size() != b.polygons.size() || a.structs.size() != b.structs.size() ||
        a.loads.size() != b.loads.size() || a.disps.size() != b.disps.size() ||
        a.hydros.size() != b.hydros.size() || a.strata.size() != b.strata.size())
        return fail("object list sizes");
    for (size_t i = 0; i < a.polygons.size(); ++i)
        if (a.polygons[i].name != b.polygons[i].name || a.polygons[i].material != b.polygons[i].material)
            return fail("region " + std::to_string(i));
    for (size_t i = 0; i < a.strata.size(); ++i)
        if (a.strata[i].name != b.strata[i].name || a.strata[i].material != b.strata[i].material)
            return fail("layer " + std::to_string(i));
    for (size_t i = 0; i < a.structs.size(); ++i)
        if (a.structs[i].name != b.structs[i].name || a.structs[i].material != b.structs[i].material ||
            a.structs[i].iface_material != b.structs[i].iface_material)
            return fail("structure " + std::to_string(i));
    for (size_t i = 0; i < a.loads.size(); ++i)
        if (a.loads[i].name != b.loads[i].name) return fail("load " + std::to_string(i));
    for (size_t i = 0; i < a.disps.size(); ++i)
        if (a.disps[i].name != b.disps[i].name) return fail("disp " + std::to_string(i));
    for (size_t i = 0; i < a.hydros.size(); ++i)
        if (a.hydros[i].name != b.hydros[i].name) return fail("well " + std::to_string(i));
    if (a.phases.size() != b.phases.size()) return fail("phase count");
    for (int k = 0; k <= (int)a.phases.size(); ++k) {
        const m::Phase& pa = k == 0 ? a.initial : a.phases[(size_t)k - 1];
        const m::Phase& pb = k == 0 ? b.initial : b.phases[(size_t)k - 1];
        const std::string at = " in phase " + std::to_string(k);
        for (size_t i = 0; i < a.polygons.size(); ++i) {
            const std::string r = "region " + std::to_string(i) + at;
            if (pa.active_poly(i) != pb.active_poly(i)) return fail(r + ": activation");
            if (pa.dry_poly(i) != pb.dry_poly(i)) return fail(r + ": dry");
            if (pa.water_mode(i) != pb.water_mode(i)) return fail(r + ": water mode");
            if (pa.water_level(i) != pb.water_level(i)) return fail(r + ": water level");
            if (pa.material_of(i, -1) != pb.material_of(i, -1)) return fail(r + ": material");
        }
        for (size_t i = 0; i < a.structs.size(); ++i) {
            const std::string s = "structure " + std::to_string(i) + at;
            if (pa.active_struct(i) != pb.active_struct(i)) return fail(s + ": activation");
            if (pa.plate_of(i, -1) != pb.plate_of(i, -1)) return fail(s + ": section");
            if (pa.contraction_of(i) != pb.contraction_of(i)) return fail(s + ": contraction");
        }
        for (size_t i = 0; i < a.loads.size(); ++i)
            if (pa.active_load(i) != pb.active_load(i)) return fail("load " + std::to_string(i) + at);
        for (size_t i = 0; i < a.disps.size(); ++i)
            if (pa.active_disp(i) != pb.active_disp(i)) return fail("disp " + std::to_string(i) + at);
        for (size_t i = 0; i < a.hydros.size(); ++i)
            if (pa.active_hydro(i) != pb.active_hydro(i)) return fail("well " + std::to_string(i) + at);
    }
    return true;
}

template <class Row>
static std::vector<Row> without(std::vector<Row> v, size_t i) {
    v.erase(v.begin() + (long)i);
    return v;
}

int main() {
    const Table T = make_table();

    // ---- 0) the comparison can say no: the old Studio deletion -------------------------------
    {
        m::Project p = build(T);
        p.structs.erase(p.structs.begin() + 1);   // what the Studio did until 2026-10-22
        Table t = T; t.structs = without(t.structs, 1);
        std::string why;
        const bool same = same_model(p, build(t), why);
        check(!same, "negative control: erasing a structure and leaving the phases is ANOTHER model (" + why + ")");
        m::Project q = build(T);
        q.polygons.erase(q.polygons.begin() + 2);
        for (auto* ph : {&q.initial, &q.phases[0], &q.phases[1], &q.phases[2]})
            if (ph->poly_active.size() > 2) ph->poly_active.erase(ph->poly_active.begin() + 2);
        Table u = T; u.regions = without(u.regions, 2);
        why.clear();
        const bool same_q = same_model(q, build(u), why);
        check(!same_q, "negative control: shifting only the activation of a deleted region is ANOTHER model (" + why + ")");
    }

    // ---- 1) every deletion of every object kind ---------------------------------------------
    struct Kind {
        const char* name;
        size_t count;
        std::function<void(m::Project&, size_t)> erase;
        std::function<void(Table&, size_t)> drop;
    };
    const Kind kinds[] = {
        {"region", T.regions.size(), [](m::Project& p, size_t i) { m::erase_polygon(p, i); },
         [](Table& t, size_t i) { t.regions = without(t.regions, i); }},
        {"structure", T.structs.size(), [](m::Project& p, size_t i) { m::erase_struct(p, i); },
         [](Table& t, size_t i) { t.structs = without(t.structs, i); }},
        {"load", T.loads.size(), [](m::Project& p, size_t i) { m::erase_load(p, i); },
         [](Table& t, size_t i) { t.loads = without(t.loads, i); }},
        {"prescribed displacement", T.disps.size(), [](m::Project& p, size_t i) { m::erase_disp(p, i); },
         [](Table& t, size_t i) { t.disps = without(t.disps, i); }},
        {"well", T.hydros.size(), [](m::Project& p, size_t i) { m::erase_hydro(p, i); },
         [](Table& t, size_t i) { t.hydros = without(t.hydros, i); }},
    };
    for (const auto& k : kinds)
        for (size_t i = 0; i < k.count; ++i) {
            m::Project p = build(T);
            k.erase(p, i);
            Table t = T;
            k.drop(t, i);
            std::string why;
            const bool same = same_model(p, build(t), why);
            check(same, std::string("delete ") + k.name + " " + std::to_string(i) +
                            " = the project built without it" + (same ? "" : " (" + why + ")"));
        }
    // ... and when the phase vectors are SHORT (files written before an object existed): the
    // missing entries read as the defaults before and after.
    {
        m::Project p = build(T);
        for (auto* ph : {&p.initial, &p.phases[0], &p.phases[1], &p.phases[2]}) {
            ph->struct_active.resize(2); ph->struct_material.resize(1); ph->struct_contraction.clear();
            ph->poly_dry.resize(1); ph->poly_material.clear();
        }
        m::Project ref = p;
        m::erase_struct(p, 4);   // beyond every short vector
        m::erase_polygon(p, 3);
        ref.structs.erase(ref.structs.begin() + 4);
        ref.polygons.erase(ref.polygons.begin() + 3);
        for (auto* ph : {&ref.initial, &ref.phases[0], &ref.phases[1], &ref.phases[2]}) {
            ph->poly_active.erase(ph->poly_active.begin() + 3);
            ph->poly_water.erase(ph->poly_water.begin() + 3);
            ph->poly_water_level.erase(ph->poly_water_level.begin() + 3);
        }
        std::string why;
        check(same_model(p, ref, why), "short phase vectors: deleting past their end changes nothing else" +
                                           (why.empty() ? std::string() : " (" + why + ")"));
    }

    // ---- 2) a region cut into two pieces: each piece is the whole, in every phase -------------
    for (size_t i = 0; i < T.regions.size(); ++i) {
        // The split replaces the region by its pieces, as region_edit does: the pieces go in
        // after the parent (insert_polygons reads `like` BEFORE the insertion), then the parent
        // is deleted.
        m::Project q = build(T);
        m::SoilPolygon a = q.polygons[i], b = q.polygons[i];
        a.name = T.regions[i].name + "a"; b.name = T.regions[i].name + "b";
        m::insert_polygons(q, i + 1, {a, b}, (int)i);
        m::erase_polygon(q, i);
        Table t = T;
        RegionRow ra = T.regions[i], rb = T.regions[i];
        ra.name = a.name; rb.name = b.name;
        t.regions.erase(t.regions.begin() + (long)i);
        t.regions.insert(t.regions.begin() + (long)i, {ra, rb});
        std::string why;
        const bool same = same_model(q, build(t), why);
        check(same, "split region " + std::to_string(i) + ": both pieces carry its settings in every phase" +
                        (same ? "" : " (" + why + ")"));
    }

    // ---- 3) every material deletion of every list --------------------------------------------
    struct List {
        const char* name;
        m::MaterialList list;
        std::vector<std::string> Table::*names;
    };
    const List lists[] = {{"soil", m::MaterialList::Soil, &Table::soils},
                          {"plate", m::MaterialList::Plate, &Table::plates},
                          {"anchor", m::MaterialList::Anchor, &Table::anchors},
                          {"geogrid", m::MaterialList::Geogrid, &Table::geogrids},
                          {"embedded beam", m::MaterialList::EmbeddedBeam, &Table::embedded}};
    for (const auto& L : lists)
        for (size_t mi = 0; mi < (T.*(L.names)).size(); ++mi) {
            m::Project p = build(T);
            const m::MaterialUses u = m::material_uses(p, L.list, mi);
            m::erase_material(p, L.list, mi);
            Table t = T;
            (t.*(L.names)).erase((t.*(L.names)).begin() + (long)mi);   // names that no longer resolve -> none
            std::string why;
            const bool same = same_model(p, build(t), why);
            check(same, std::string("delete ") + L.name + " material " + std::to_string(mi) + " (" +
                            std::to_string(u.total()) + " use(s)) = the project whose rows name a missing one" +
                            (same ? "" : " (" + why + ")"));
        }
    // material_uses counts what the builder wrote
    {
        const m::Project p = build(T);
        int all_changes = 0;
        for (size_t mi = 0; mi < T.soils.size(); ++mi) {
            const std::string& n = T.soils[mi];
            const m::MaterialUses u = m::material_uses(p, m::MaterialList::Soil, mi);
            int regions = 0, layers = 0, ifaces = 0, changes = 0;
            for (const auto& r : T.regions) {
                regions += r.material == n;
                for (const auto& s : r.phase_material) changes += s == n;
            }
            for (const auto& L : T.layers) layers += L.second == n;
            for (const auto& s : T.structs) ifaces += s.iface_material == n;
            all_changes += changes;
            check(u.regions == regions && u.layers == layers && u.interfaces == ifaces && u.phase_changes == changes,
                  "material_uses(" + n + "): " + std::to_string(u.regions) + " region(s), " + std::to_string(u.layers) +
                      " layer(s), " + std::to_string(u.interfaces) + " interface(s), " +
                      std::to_string(u.phase_changes) + " phase change(s), as written");
        }
        check(all_changes > 0, "the table changes a region's material in some phase (the count is exercised)");
        int sections = 0;
        for (size_t mi = 0; mi < T.plates.size(); ++mi)
            sections += m::material_uses(p, m::MaterialList::Plate, mi).phase_changes;
        int want = 0;
        for (const auto& s : T.structs)
            for (const auto& sec : s.section) want += !sec.empty();
        check(sections == want && want > 0, "material_uses counts the plate section changes: " + std::to_string(sections));
    }

    // ---- 4) the vertices of a region keep their edge conditions --------------------------------
    {
        m::SoilPolygon P;
        P.x = {0, 4, 4, 0}; P.y = {0, 0, 3, 3};
        P.edge_bc = {1, 2, 3, 4};
        P.edge_flow = {0, 3, 1, 2};
        P.edge_head = {0.0, 0.0, 5.0, 0.0};
        P.edge_flux = {0.0, 0.25, 0.0, 0.0};
        m::split_edge(P, 1);   // (4,0)-(4,3) at (4,1.5)
        check(P.x.size() == 5 && P.x[2] == 4.0 && P.y[2] == 1.5, "split: the new vertex is the edge's midpoint");
        check(P.edge_bc == std::vector<int>({1, 2, 2, 3, 4}) && P.edge_flow == std::vector<int>({0, 3, 3, 1, 2}) &&
                  P.edge_flux == std::vector<double>({0.0, 0.25, 0.25, 0.0, 0.0}) &&
                  P.edge_head == std::vector<double>({0.0, 0.0, 0.0, 5.0, 0.0}),
              "split: both halves keep the edge's conditions, its prescribed flux included");
        m::erase_vertex(P, 2);
        check(P.x.size() == 4 && P.edge_bc == std::vector<int>({1, 2, 3, 4}) &&
                  P.edge_flux == std::vector<double>({0.0, 0.25, 0.0, 0.0}),
              "delete the vertex again: the polygon and its edge conditions are as drawn");
        m::SoilPolygon Q;
        Q.x = {0, 4, 4}; Q.y = {0, 0, 3};
        m::split_edge(Q, 2);   // no edge conditions stated at all: nothing to copy, nothing read
        check(Q.x.size() == 4 && Q.edge_bc.empty() && Q.edge_flow.empty() && Q.edge_flux.empty(),
              "split with no edge conditions stated: the lists stay empty (the Studio read past their end)");
        m::erase_vertex(Q, 0);
        m::erase_vertex(Q, 0);   // refused: a region keeps three vertices
        check(Q.x.size() == 3, "a region keeps at least three vertices");
    }

    // ---- 5) ground regenerated from boreholes ---------------------------------------------------
    {
        m::Project pr;
        pr.x_min = 0; pr.x_max = 40; pr.y_min = 0; pr.y_max = 20;
        pr.materials.resize(2);
        pr.strata = {{"Fill", 0}, {"Clay", 1}};
        m::Borehole b; b.x = 10; b.level = {20, 15, 0}; b.head = 18;
        pr.boreholes = {b};
        std::string msg;
        check(m::apply_stratigraphy(pr, &msg), "boreholes: the ground is generated");
        pr.phases.resize(1);
        pr.phases[0].poly_active = {0, 1};
        pr.phases[0].poly_dry = {1, 0};
        pr.boreholes[0].level = {20, 12, 0};   // the same layers, re-levelled
        msg.clear();
        check(m::apply_stratigraphy(pr, &msg) && !pr.phases[0].active_poly(0) && pr.phases[0].dry_poly(0),
              "boreholes re-levelled: the same regions, so the phase keeps its excavation");
        pr.strata.push_back({"Rock", 0});
        pr.boreholes[0].level = {20, 12, 5, 0};   // a layer added: other regions
        msg.clear();
        const bool ok = m::apply_stratigraphy(pr, &msg);
        check(ok && pr.phases[0].poly_active.empty() && pr.phases[0].poly_dry.empty() &&
                  msg.find("were reset") != std::string::npos,
              "boreholes with a new layer: the region settings are reset and the message says so");
    }

    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
