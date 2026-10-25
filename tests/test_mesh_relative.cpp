// A relative element size (`mesh.rel`): the target element edge as a factor on a twentieth of the
// model's diagonal, so a model drawn wider keeps the same density of elements and a first mesh
// needs no number in metres. It is an input convenience, not a new discretisation: what it must
// do is stand for the absolute size it means, everywhere the element size is read.
//
// verify: none -- an input convenience, held by identity with the absolute element size it stands for (the same mesh, node for node), by the file round trip, and by the validator's two refusals
#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

m::Project block(double elem_size, double rel) {
    m::Project pr;
    pr.name = "relative mesh";
    pr.x_min = 0; pr.x_max = 40; pr.y_min = 0; pr.y_max = 30;   // diagonal 50 m
    m::Material s; s.name = "Soil"; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Block"; P.material = 0;
    P.x = {0, 40, 40, 0}; P.y = {0, 0, 30, 30};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.mesh.elem_size = elem_size;
    pr.mesh.rel = rel;
    return pr;
}
bool has_issue(const katai::io::ValidationReport& r, const char* path) {
    for (const auto& i : r.issues)
        if (i.path == path && i.severity == katai::io::Severity::Error) return true;
    return false;
}
}  // namespace

int main() {
    std::printf("== relative element size (mesh.rel) ==\n");
    // A factor of 1 on a 40 x 30 m model is a twentieth of its 50 m diagonal: 2.5 m. The absolute
    // size in the file is deliberately something else, so the identity below can only hold if the
    // relative one is the one read.
    const m::Project rel = block(7.0, 1.0), absolute = block(2.5, 0.0);
    check(std::fabs(m::element_size(rel) - 2.5) < 1e-12, "factor 1 on a 50 m diagonal is 2.5 m");
    check(m::element_size(absolute) == 2.5, "with no relative size the absolute one is the size");
    check(std::fabs(m::element_size(block(7.0, 0.5)) - 1.25) < 1e-12, "factor 0.5 halves it");

    const auto Mr = katai::app::mesh_from_project(rel), Ma = katai::app::mesh_from_project(absolute);
    check(Mr.ok && Ma.ok, "both mesh");
    if (Mr.ok && Ma.ok) {
        bool same = Mr.mesh.node_count == Ma.mesh.node_count &&
                    Mr.mesh.element_count == Ma.mesh.element_count;
        for (int n = 0; same && n < Mr.mesh.node_count; ++n)
            same = Mr.mesh.x[n] == Ma.mesh.x[n] && Mr.mesh.y[n] == Ma.mesh.y[n];
        std::printf("  relative: %d nodes, %d elements; absolute: %d nodes, %d elements\n",
                    Mr.mesh.node_count, Mr.mesh.element_count, Ma.mesh.node_count,
                    Ma.mesh.element_count);
        check(same, "the relative size meshes node for node as the absolute size it stands for");
        const auto Mf = katai::app::mesh_from_project(block(7.0, 0.5));
        check(Mf.ok && Mf.mesh.element_count > 3 * Mr.mesh.element_count,
              "a finer factor gives a finer mesh");
    }

    // The file: written only when set, and read back.
    const std::string jr = m::project_to_json(rel), ja = m::project_to_json(absolute);
    check(jr.find("\"rel\":1") != std::string::npos, "a relative size is written");
    check(ja.find("\"rel\"") == std::string::npos, "no relative size, no key (older files unchanged)");
    m::Project back;
    std::string err;
    check(m::project_from_json(jr, back, &err) && back.mesh.rel == 1.0 &&
              back.mesh.elem_size == 7.0,
          "it is read back, the absolute size beside it kept");

    // The validator: a negative factor, and a factor without an extent to take it from.
    check(has_issue(katai::io::validate_project(block(2.0, -1.0)), "mesh.rel"),
          "a negative relative size is an error");
    m::Project flat = block(2.0, 1.0);
    flat.x_max = flat.x_min; flat.y_max = flat.y_min;
    check(has_issue(katai::io::validate_project(flat), "mesh.rel"),
          "a relative size with an empty model extent is an error");
    check(!has_issue(katai::io::validate_project(rel), "mesh.rel"), "the valid one passes");

    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
