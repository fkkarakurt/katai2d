// GENERAL interface (any orientation, soil-soil / standalone) -- the engine behind PLAXIS-style
// interfaces on ANY structure line or a free slip surface. split_mesh_at_segment duplicates the mesh
// nodes along an arbitrary segment; build_soil_interface joins the two sides with a Coulomb joint
// (the duplicate nodes carry the negative side's soil). This test proves the mechanism is REAL --
// i.e. "interface present" gives DIFFERENT results from "absent" -- which is exactly the property a
// user checks when they toggle an interface on/off:
//   (1) RIGID interface (huge kn/ks, never yields) == the bonded continuum (unsplit mesh) to a tight
//       tolerance -> the split + interface assembly is correct (no spurious flexibility).
//   (2) SOFT interface (finite ks) gives MORE relative slip than bonded, and a clearly nonzero
//       tangential slip across the seam -> the interface actually transmits/relaxes shear.
//   (3) Both hold for tri6 AND tri15 (3- and 5-node interface edges).
// Refs: interface-formulation.md; staged-construction (split). Core physics: test_interface (closed form).
#include <katai/analysis/general_interface.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace iface = katai::core::iface;
using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Structures;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
Eigen::VectorXd solve_ns(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    try {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric); s->factorize(k); return s->solve(r);
    } catch (...) { return Eigen::VectorXd::Zero(r.size()); }
}

constexpr double W = 2.0, H = 2.0;
constexpr int NX = 4, NY = 4;
constexpr double E = 1.0e4, nu = 0.3, q = 100.0;   // downward traction on the RIGHT half top

Mesh make_mesh(int order) {
    RectangularDomain domain{0.0, 0.0, W, H, 0};
    return order == 15 ? katai::mesh::generate_structured_tri15(domain, NX, NY)
                       : katai::mesh::generate_structured_tri6(domain, NX, NY);
}
// Top edge chain on the RIGHT half (x >= W/2), ordered by x (corner, mid, ... -- a valid traction chain).
std::vector<int> right_top_chain(const Mesh& mesh) {
    std::vector<int> c;
    for (int n : mesh.top_nodes) if (mesh.x[n] > 0.5 * W - 1e-6) c.push_back(n);
    std::sort(c.begin(), c.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    return c;
}
int right_top_corner(const Mesh& mesh) {  // node at (W, H)
    int best = mesh.top_nodes[0];
    for (int n : mesh.top_nodes) if (mesh.x[n] > mesh.x[best]) best = n;
    return best;
}

// Solve the differential-load block. split=false -> bonded continuum. Otherwise a vertical interface at
// x=W/2 with props *ip. Returns settlement |uy| at the right-top corner; *slip = max tangential slip.
double solve_block(int order, bool split, const iface::InterfaceProps* ip, double* slip) {
    Mesh mesh = make_mesh(order);
    const int rt_corner_orig = right_top_corner(mesh);   // (W,H) -- unaffected by an interior split
    Structures structures;
    katai::core::InterfaceRange rng{order, 0, 0};
    std::vector<katai::core::SegSeam> seam;
    if (split) {
        // split the FULL seam (both endpoints included) so the chain runs corner,mid,corner,... -> the
        // interface-edge triples/quintuples are aligned (a mid-node start would scramble them).
        seam = katai::core::split_mesh_at_segment(mesh, 0.5 * W, 0.0, 0.5 * W, H, -1.0, H + 1.0);
    }
    DofMap dofs(mesh.node_count, 2);
    if (split && ip) rng = katai::core::build_soil_interface(mesh, seam, dofs, *ip, order, structures);
    // BCs by COORDINATE (after split) so BOTH the original and the duplicate boundary nodes are fixed.
    for (int n = 0; n < mesh.node_count; ++n) {
        if (mesh.y[n] < 1e-9) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }  // base
        if (mesh.x[n] < 1e-9) dofs.fix_node_component(n, 0);                                     // left lateral
    }
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, right_top_chain(mesh), 0.0, -q, f);
    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, E, nu}};
    katai::core::NewtonOptions opt; opt.load_steps = 1; opt.max_iterations = 60; opt.tolerance = 1e-9;
    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_ns, opt, {}, {}, structures);

    if (slip) {   // max |uy_orig - uy_dup| across the seam (tangential slip on the vertical interface)
        *slip = 0.0;
        for (const auto& sp : seam) {
            const double du = r.displacement[dofs.global_dof(sp.orig, 1)] - r.displacement[dofs.global_dof(sp.dup, 1)];
            *slip = std::fmax(*slip, std::fabs(du));
        }
    }
    return std::fabs(r.displacement[dofs.global_dof(rt_corner_orig, 1)]);
}

void run(int order) {
    std::printf("-- element order %d --\n", order);
    const double bonded = solve_block(order, false, nullptr, nullptr);

    // WELDED: stiff + no tension cut-off (sigma_t huge) -> mimics the continuum exactly.
    iface::InterfaceProps weld; weld.kn = 1e9; weld.ks = 1e6; weld.c_i = 1e12; weld.phi_i = 0.0; weld.sigma_t = 1e12;
    double slip_weld = 0.0;
    const double u_weld = solve_block(order, true, &weld, &slip_weld);

    // SOFT: finite shear stiffness + a real tension cut-off (sigma_t=0) -> the interface slips tangentially
    // AND opens in normal tension, so the result clearly differs from bonded (this is what a user sees when
    // they toggle an interface on a structure).
    iface::InterfaceProps soft; soft.kn = 1e6; soft.ks = 5e2; soft.c_i = 1e12; soft.phi_i = 0.0; soft.sigma_t = 0.0;
    double slip_soft = 0.0;
    const double u_soft = solve_block(order, true, &soft, &slip_soft);

    std::printf("   right-top |uy|: bonded=%.6e  welded=%.6e (%+.2f%%)  soft=%.6e (%+.1f%%)\n",
                bonded, u_weld, 100.0 * (u_weld - bonded) / bonded, u_soft, 100.0 * (u_soft - bonded) / bonded);
    std::printf("   seam slip:      welded=%.3e   soft=%.3e\n", slip_weld, slip_soft);

    check(bonded > 0.0, "bonded continuum settles (sanity)");
    check(std::fabs(u_weld - bonded) < 0.01 * bonded, "WELDED interface (stiff, no cut-off) == bonded continuum");
    check(slip_weld < 0.02 * bonded, "welded interface: slip -> 0 (sides move together)");
    check(std::fabs(u_soft - bonded) > 0.03 * bonded, "SOFT interface changes the result (present != absent)");
    check(slip_soft > 0.1 * bonded, "soft interface: clear tangential slip on the seam");
    check(slip_soft > 50.0 * slip_weld, "soft slip >> welded slip (interface relaxes shear)");
}

}  // namespace

int main() {
    std::printf("General interface (any orientation, soil-soil) -- split + Coulomb joint\n\n");
    run(6);
    run(15);
    if (g_failures == 0) {
        std::printf("\nOK: general interface verified -- rigid==bonded continuum, soft interface slips "
                    "(present != absent), tri6 + tri15\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
