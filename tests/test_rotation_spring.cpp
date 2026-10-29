// A rotation spring on a hinged plate end: a semi-rigid connection.
//
// A hinged end keeps its own rotation (KV-STR-025). A spring on it ties that rotation to the one it
// is hinged against -- the ground where the end stands on a support, the joint where another
// plate carries on -- by a moment M = k (phi_end - phi_joint).
//
// THE ORACLE IS THE DIFFERENCE. A cantilever plate in free space (its ground switched off) under a
// point load P at its tip deflects by bending, by the Timoshenko shear of its section, and by the
// turn of its base. Clamped, the base does not turn; on a spring of stiffness k it turns by the
// base moment over k, P L / k, and the tip moves by that times L. Everything else is the same
// plate under the same load, so the difference of the two tip deflections is P L^2 / k, exactly,
// whatever the element and its shear stiffness. A spring at a joint part-way up -- the upper
// plate hinged on the lower one -- turns by P L_B / k under the moment at the joint, and adds
// P L_B^2 / k at the tip.
//
// verify: KV-STR-028
//   oracle:   closed_form
//   source:   statics and the rotation of a linear rotational spring under the moment it carries: a cantilever's base moment is P L, and a rigid-body turn theta of a member of length L moves its tip by theta L (any structural analysis text, e.g. Hibbeler, R.C. Structural Analysis -- semi-rigid connections)
//   locator:  (a), (b) a vertical plate 6 m long, EA = 1e6 kN/m, EI = 1e4 kNm2/m, standing on the fixed base of a 20 x 10 m block whose soil is switched off in the loaded phase; a horizontal point load P = 10 kN/m at its top. (a) its base clamped, then hinged on a spring k = 1e4 and 1e3 kNm/rad/m, then on k = 1e12; (b) drawn as two plates meeting at 3 m, rigid, then the upper one hinged on the lower one with k = 1e4; (c) two walls with interfaces in elastic ground (KV-STR-025 (c)), the upper hinged on the lower with k = 1e3, 1e5 and 1e12
//   quantity: the horizontal displacement of the top [m]; (c) the bending moment at the joint [kNm/m]
//   expected: (a) u(k) - u(clamped) = P L^2 / k = 0.036 m at k = 1e4, 0.36 m at k = 1e3, and the clamped value at k = 1e12; (b) u(k) - u(rigid) = P L_B^2 / k = 0.009 m; (c) the moment grows with k between the hinge's and the rigid joint's, and is the rigid joint's at k = 1e12 (a bound and a limit, not a closed form: the wall's joint is wired the same way as the plate's)
//   band:     1e-6 relative on the differences (a linear structure: the solver's precision; measured to the printed 9 digits); 1e-6 of the clamped deflection at k = 1e12; (c) 1e-3 on the rigid limit
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

constexpr double kP = 10.0, kL = 6.0, kEI = 1.0e4;

m::StructElement plate(const char* name, double x1, double y1, double x2, double y2) {
    m::StructElement e;
    e.kind = m::StructKind::Plate; e.name = name;
    e.x1 = x1; e.y1 = y1; e.x2 = x2; e.y2 = y2;
    e.material = 0;
    return e;
}

// The block, its soil on in the initial phase and off in the loaded one; the plates and the load
// only in the loaded one. No boundary condition is drawn, so the standard fixities apply: the
// base is fixed, and a plate standing on it is clamped there unless its end is hinged.
double top_displacement(const std::vector<m::StructElement>& plates, const char* what) {
    m::Project pr;
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    m::Material s;
    s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
    s.E = 10000.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon g;
    g.name = "Ground"; g.material = 0; g.x = {0, 20, 20, 0}; g.y = {0, 0, 10, 10};
    g.edge_bc.assign(4, 0);
    pr.polygons.push_back(g);
    m::PlateMaterial pm;
    pm.name = "Column"; pm.EA = 1.0e6; pm.EI = kEI; pm.nu = 0.0;
    pr.plates.push_back(pm);
    pr.structs = plates;
    m::Load L;
    L.kind = m::LoadKind::Point; L.name = "P";
    L.x1 = 10.0; L.y1 = kL; L.x2 = 10.0; L.y2 = kL;
    L.qx1 = L.qx2 = kP; L.qy1 = L.qy2 = 0.0;
    pr.loads.push_back(L);
    pr.initial.struct_active.assign(plates.size(), 0);
    pr.initial.load_active = {0};
    m::Phase ph = pr.initial;
    ph.name = "Load";
    ph.poly_active = {0};
    ph.struct_active.assign(plates.size(), 1);
    ph.load_active = {1};
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { check(false, std::string(what) + ": meshes"); return 0.0; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) {
        check(false, std::string(what) + ": solves (" + (res.size() == 2 ? res[1].message : std::string("no phase")) + ")");
        return 0.0;
    }
    for (const auto& f : res[1].struct_forces)
        for (const auto& st : f.stations)
            if (std::fabs(st.x - 10.0) < 1e-9 && std::fabs(st.y - kL) < 1e-9) return st.ux;
    check(false, std::string(what) + ": reports the top");
    return 0.0;
}

}  // namespace

int main() {
    std::printf("== a rotation spring on a hinged plate end (KV-STR-028) ==\n");
    {
        const double u0 = top_displacement({plate("Column", 10, 0, 10, kL)}, "clamped");
        std::printf("      clamped: top %.8f m (bending alone P L^3 / 3 EI = %.8f)\n", u0, kP * kL * kL * kL / (3.0 * kEI));
        for (double k : {1.0e4, 1.0e3, 1.0e12}) {
            m::StructElement c = plate("Column", 10, 0, 10, kL);
            c.hinge1 = true;
            c.rot_k1 = k;
            const double u = top_displacement({c}, "on a spring");
            const double expected = kP * kL * kL / k;
            std::printf("      base spring k = %.0e: top %.8f m, difference %.8e (P L^2 / k = %.8e)\n", k, u, u - u0, expected);
            if (k < 1e9)
                check(rel(u - u0, expected) < 1e-6, "a base on a spring turns by P L / k: the top moves P L^2 / k more than clamped");
            else
                check(std::fabs(u - u0) < 1e-6 * std::fabs(u0), "...and a very stiff spring is the clamped base");
        }
    }
    {
        const double LB = kL - 3.0;
        const double u0 = top_displacement({plate("Lower", 10, 0, 10, 3), plate("Upper", 10, 3, 10, kL)}, "in two pieces, rigid");
        m::StructElement up = plate("Upper", 10, 3, 10, kL);
        up.hinge1 = true;
        up.rot_k1 = 1.0e4;
        const double u = top_displacement({plate("Lower", 10, 0, 10, 3), up}, "in two pieces, a spring at the joint");
        const double expected = kP * LB * LB / up.rot_k1;
        std::printf("      joint spring: top %.8f m against %.8f rigid, difference %.8e (P L_B^2 / k = %.8e)\n", u, u0, u - u0, expected);
        check(rel(u - u0, expected) < 1e-6, "a spring at a joint turns the upper plate by P L_B / k: its top moves P L_B^2 / k more");
    }
    {
        // (c) Two walls with interfaces in the ground (KV-STR-025 (c)), the upper hinged on the
        // lower with a spring. The joint's moment has no closed form here; what the spring must do
        // is carry a moment between the hinge's none and the rigid joint's, more the stiffer it is,
        // and the rigid joint's when it is very stiff. Unwired, it would carry none.
        const auto joint_moment = [](bool hinged, double k) {
            m::Project pr;
            pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
            pr.has_water = false;
            pr.initial_procedure = m::InitialProcedure::GravityLoading;
            pr.mesh.elem_size = 1.0;
            pr.mesh.auto_refine = false;
            m::Material s;
            s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
            s.E = 10000.0; s.nu = 0.3; s.gamma_unsat = s.gamma_sat = 0.0;
            pr.materials.push_back(s);
            m::SoilPolygon g;
            g.name = "Ground"; g.material = 0; g.x = {0, 20, 20, 0}; g.y = {0, 0, 10, 10};
            g.edge_bc.assign(4, 0);
            pr.polygons.push_back(g);
            m::PlateMaterial pm;
            pm.name = "Wall"; pm.EA = 5.0e6; pm.EI = 5.0e4; pm.nu = 0.0;
            pr.plates.push_back(pm);
            pr.structs.push_back(plate("Upper", 10, 9, 10, 6));
            pr.structs.push_back(plate("Lower", 10, 6, 10, 3));
            for (auto& st : pr.structs) st.iface_pos = st.iface_neg = true;
            pr.structs[0].hinge2 = hinged;
            pr.structs[0].rot_k2 = k;
            m::Load L;
            L.kind = m::LoadKind::Point; L.name = "P";
            L.x1 = 10; L.y1 = 9; L.qx1 = 40; L.qy1 = 0;
            pr.loads.push_back(L);
            pr.initial.load_active = {1};
            const auto M = katai::app::mesh_from_project(pr);
            if (!M.ok) return -1.0;
            const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
            if (res.empty() || !res[0].ok) return -1.0;
            for (const auto& f : res[0].struct_forces) {
                if (f.name != "Lower") continue;
                for (const auto& st : f.stations)
                    if (std::fabs(st.x - 10.0) < 1e-6 && std::fabs(st.y - 6.0) < 1e-6) return std::fabs(st.M);
            }
            return -1.0;
        };
        const double rigid = joint_moment(false, 0.0), free = joint_moment(true, 0.0);
        const double m3 = joint_moment(true, 1.0e3), m5 = joint_moment(true, 1.0e5), m12 = joint_moment(true, 1.0e12);
        std::printf("      walls with interfaces, M at the joint: hinged %.4f, k = 1e3 %.4f, 1e5 %.4f, 1e12 %.4f, rigid %.4f kNm/m\n",
                    free, m3, m5, m12, rigid);
        check(rigid > 1.0 && free >= 0.0 && m3 > free && m5 > m3 && m12 > m5 && m5 < rigid,
              "walls with interfaces: the joint's moment grows with the spring, between the hinge's and the rigid joint's");
        check(rigid > 0.0 && rel(m12, rigid) < 1e-3, "...and a very stiff spring is the rigid joint");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
