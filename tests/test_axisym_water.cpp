// A water table in AXISYMMETRIC mode. Until now the driver refused it -- "use plane strain, or
// remove the water table" -- which closed every circular problem that has groundwater in it: the
// tank, the silo, the shaft, the circular footing, the pile load test. The refusal was honest but
// it was not small.
//
// What was missing turned out to be two assemblies and one term inside one of them. The K0 seed
// was already effective-stress and already set the hoop; the phreatic body force only needed its
// r weight; and the pore-pressure load needed something that is NOT in the plane-strain version at
// all. In plane strain the pore load is integral B^T (u m) dA with m = [1, 1, 0]: the out-of-plane
// direction carries no equation, so there is nothing to put there. In axisymmetry the strain has
// four components, [eps_r, eps_z, gamma_rz, eps_theta], the hoop is a real strain with a real
// equation, and pore pressure is isotropic -- so m = [1, 1, 0, 1] and the hoop term lands on the
// RADIAL degree of freedom next to dN/dr.
//
// That term is exactly what a careful-looking copy of the plane-strain function would drop, and
// dropping it produces no error message and no obviously wrong picture. So the test is built
// around the identity it breaks.
//
// THE IDENTITY, AND WHERE IT IS EXACT. A K0 state is a state of equilibrium: the internal force
// of the seeded stresses must equal the external force holding them up, which in effective-stress
// form is the phreatic body force PLUS the pore-pressure load. That is algebra rather than a
// solution -- but only where the quadrature can integrate both sides. In axisymmetry the r weight
// makes the radial integrand one degree higher than the plane-strain one, and the 3-point rule of
// the 6-noded triangle cannot take a cubic exactly; the driver already knows this and refuses to
// read a nil-step from an assembled axisymmetric imbalance for that reason. Measured here, and it
// is the measurement that tells the two apart:
//
//     tri6    residual 5.26e-03      without the hoop term 3.04e-01        58x
//     tri15   residual 8.13e-13      without the hoop term 3.07e-01     4e+11x
//
// So the identity comes back to round-off the moment the rule can integrate it, which is what
// says the tri6 residual is the quadrature and not the physics -- and the hoop term is worth
// eleven orders of magnitude, not a correction.
//
// Section (2) runs the same cylinder through the driver. The K0 phase alone cannot fail there --
// on level ground with a level table nothing is ramped, so a zero displacement is arithmetic
// rather than evidence -- so what is asked of it is its stresses and its pore field, and then a
// second phase LOWERS the water table, which is a real staged change with a closed form of its
// own and in which the only thing that changes is the water.
//
// verify: KV-CST-013
//   oracle:   closed_form
//   source:   Terzaghi effective stress with a hydrostatic water table (docs/references/effective-stress-formulation.md); the axisymmetric equilibrium of a horizontally layered state -- with sigma_r = sigma_theta the (sigma_r - sigma_theta)/r term of the radial equation vanishes, so the K0 state is in equilibrium in r-z exactly as it is in plane strain; and the one-dimensional compression of a laterally confined column under the effective-stress increase that lowering a water table produces. Axisymmetric kinematics: strain [eps_r, eps_z, gamma_rz, eps_theta], r-weighted integration, kernel/fem/elements/axisymmetric.hpp
//   locator:  sigma'_v(z) = -[gamma_unsat (H - z_wt) + (gamma_sat - gamma_w)(z_wt - z)] below the table and -gamma_unsat (H - z) above it; sigma'_r = sigma'_theta = K0 sigma'_v; u(z) = gamma_w max(0, z_wt - z). Equilibrium: integral B^T sigma' r dA = f_gravity (gamma_sat below / gamma_unsat above, r-weighted) + integral B^T (u m) r dA with m = [1, 1, 0, 1]. Dewatering from z1 to z2: every metre of lowering adds (gamma_unsat - gamma_sat + gamma_w) to sigma'_v below it, so s = (gamma_unsat - gamma_sat + gamma_w)/Eoed * [z2 (z1 - z2) + (z1 - z2)^2 / 2], Eoed = E(1-nu)/((1+nu)(1-2nu))
//   quantity: the relative residual ||f_int(sigma_K0) - (f_gravity + f_pore)|| / ||f_gravity|| of the assembled axisymmetric system, on both element orders and with the plane-strain pore load substituted [-]; and end to end, sigma'_v and the pore pressure at three depths [kPa] plus the settlement of lowering the table by 3 m [m]
//   expected: the residual is an identity where the quadrature can integrate it -- round-off on tri15, the known r-weighted cubic residue on tri6 -- and substituting the plane-strain pore load must break it by orders rather than by a correction. sigma'_v = the buoyant closed form, the pore pressure = gamma_w (z_wt - z), and the dewatering settlement = 4.173557e-03 m for this column (Eoed = 26923.08 kPa, 6.81 kPa per metre of lowering, 16.5 m^2 of depth integral)
//   band:     1e-12 on the tri15 identity, measured 8.13e-13; on tri6 the residual is asserted to LIE BETWEEN 1e-4 and 5e-2 (measured 5.26e-03) because a round-off there would mean the quadrature had changed and a large one that the physics had. The hoop term is asserted only as a ratio, > 20x, since its real margin is 58x on tri6 and 4e+11x on tri15 and the point is that it is never small. End to end: 0.5% on sigma'_v (nodal interpolation of a Gauss-point field), 1e-6 kPa on the pore pressure (it is read from the same closed form the load was built from, so it is a plumbing check), 1% on the dewatering settlement (measured -0.05%; the residue is the FE discretisation of a 1-D closed form on a 0.5 m mesh) and max |u_r| < 2% of it (measured 0.08%)
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;
using katai::core::DofMap;
using katai::core::GaussState;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kR = 5.0, kH = 10.0;        // cylinder: radius, height [m]
constexpr double kWt = 7.0;                  // water table elevation [m]
constexpr double kWtLow = 4.0;               // where a later phase lowers it to [m]
constexpr double kGunsat = 17.0, kGsat = 20.0;
constexpr double kE = 20000.0, kNu = 0.3;
constexpr double kK0 = 0.5;

double gamma_w() { return katai::app::kGammaWater; }
double pore_at(double z) { return gamma_w() * std::fmax(0.0, kWt - z); }
// Effective vertical stress (compression negative), stated here rather than taken from the K0 code.
double sigma_v_eff(double z) {
    if (z >= kWt) return -kGunsat * (kH - z);
    return -(kGunsat * (kH - kWt) + (kGsat - gamma_w()) * (kWt - z));
}

// The K0 seed from the closed form above, for either element order (compression negative;
// stress_zz is the hoop, which carries the same K0 stress as the radial one).
template <class E>
void seed_k0_impl(const Mesh& mesh, std::vector<GaussState>& states) {
    states.assign((size_t)mesh.element_count * E::kGaussCount, GaussState{});
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        typename E::NodeCoords X;
        for (int k = 0; k < E::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            X(k, 0) = mesh.x[n]; X(k, 1) = mesh.y[n];
        }
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto N = E::shape_functions(gp[g].xi, gp[g].eta);
            double z = 0.0;
            for (int k = 0; k < E::kNodeCount; ++k) z += N(k) * X(k, 1);
            const double sv = sigma_v_eff(z), sh = kK0 * sv;
            GaussState& s = states[(size_t)e * E::kGaussCount + g];
            s.stress << sh, sv, 0.0;
            s.stress_zz = sh;
        }
    }
}
void seed_k0(const Mesh& mesh, int order, std::vector<GaussState>& states) {
    if (order == 15) seed_k0_impl<katai::core::Tri15Element>(mesh, states);
    else             seed_k0_impl<katai::core::Tri6Element>(mesh, states);
}

// ---------------------------------------------------------------------------------------------
// (1) The identity, on the assembled vectors.
// ---------------------------------------------------------------------------------------------
void test_equilibrium_identity() {
    std::printf("-- (1) the K0 state is in equilibrium in r-z, and the hoop pore term is why --\n");
    for (int order : {6, 15}) {
        RectangularDomain domain{0.0, 0.0, kR, kH, 0};
        Mesh mesh = order == 15 ? katai::mesh::generate_structured_tri15(domain, 5, 10)
                                : katai::mesh::generate_structured_tri6(domain, 5, 10);
        DofMap dofs(mesh.node_count, 2);
        for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);   // the axis: u_r = 0 by symmetry
        for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
        dofs.finalize();

        // The seeded K0 state, from the closed form above (compression negative; stress_zz = hoop).
        std::vector<GaussState> states;
        seed_k0(mesh, order, states);

        const int neq = dofs.equation_count();
        const std::vector<double> gu{kGunsat}, gs{kGsat};
        const auto wt = [](double) { return kWt; };
        const auto pore = [](double, double z) { return pore_at(z); };

        Eigen::VectorXd fint = Eigen::VectorXd::Zero(neq);
        katai::core::assemble_axisym_internal_force(mesh, dofs, states, fint);
        Eigen::VectorXd fgrav = Eigen::VectorXd::Zero(neq);
        katai::core::assemble_axisym_gravity_phreatic(mesh, dofs, gu, gs, wt, fgrav);
        Eigen::VectorXd fext = fgrav;
        katai::core::assemble_axisym_pore_pressure_load(mesh, dofs, pore, fext);

        const double scale = fgrav.norm();
        const double resid = (fint - fext).norm() / scale;
        // The control: the SAME assembly with the plane-strain pore load, which has no hoop term.
        Eigen::VectorXd fext_ps = fgrav;
        katai::core::assemble_pore_pressure_load(mesh, dofs, pore, fext_ps);
        const double resid_ps = (fint - fext_ps).norm() / scale;

        std::printf("   tri%-2d  residual %.3e   without the hoop term %.3e   (%.0fx)\n",
                    order, resid, resid_ps, resid_ps / std::fmax(resid, 1e-300));
        check(scale > 1.0, "the model has weight to balance");
        if (order == 15) {
            // THE IDENTITY IS EXACT WHERE THE QUADRATURE CAN INTEGRATE IT. The r weight makes the
            // radial integrand d(rN)/dr * sigma one degree higher than the plane-strain one, and
            // the 3-point tri6 rule cannot take a cubic exactly -- the driver already knows this
            // and refuses to read a nil-step from an assembled axisymmetric imbalance for exactly
            // that reason. The 15-noded element's rule can, and there the identity comes back to
            // round-off: that is what says the residual on tri6 is the QUADRATURE and not the
            // physics.
            check(resid < 1e-12, "tri15: the state balances its own weight and water to round-off");
        } else {
            check(resid > 1e-4 && resid < 5e-2,
                  "tri6: what is left is the known r-weighted quadrature residual, not an imbalance");
        }
        check(resid_ps > 20.0 * resid,
              "and WITHOUT the hoop term the imbalance is far louder -- the term is load-bearing");
    }
}

// ---------------------------------------------------------------------------------------------
// (2) The same column through the driver, end to end.
// ---------------------------------------------------------------------------------------------
void test_through_the_driver() {
    std::printf("\n-- (2) the same cylinder through the calculation path --\n");
    // The K0 phase alone cannot fail here: on level ground with a level table nothing is ramped,
    // so its displacement is zero by construction rather than by equilibrium. What it CAN be
    // asked is whether the stresses it seeded and the pore field it reports are the right ones,
    // and then a second phase LOWERS the water table -- which is a real staged change with a
    // closed form, and the only thing that changes in it is the water.
    m::Project pr;
    pr.axisymmetric = true;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu; s.gamma_unsat = kGunsat; s.gamma_sat = kGsat;
    s.k0_auto = false; s.k0 = kK0;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kR, kR, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed};
    P.edge_head = {0, 0, 0, 0};
    pr.polygons.push_back(P);
    pr.has_water = true;
    pr.wx = {0.0, kR};
    pr.wy = {kWt, kWt};

    m::Phase dw; dw.name = "Dewatering"; dw.type = m::PhaseType::Plastic;
    dw.water_override = true;
    dw.wx = {0.0, kR};
    dw.wy = {kWtLow, kWtLow};
    pr.phases.push_back(dw);

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "the cylinder meshed");
    if (!M.ok) { std::printf("   (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2 && res[0].ok,
          "an axisymmetric model WITH a water table now solves (it used to be refused)");
    if (res.size() != 2 || !res[0].ok) {
        if (!res.empty()) std::printf("   (%s)\n", res[0].message.c_str());
        return;
    }
    const auto& K0 = res[0];

    // Effective stress: buoyant below the table, moist above. Read at nodes.
    std::printf("   z [m]   sigma'_v FE      closed form     pore FE     closed form\n");
    int checked = 0;
    for (double zq : {2.0, 5.0, 9.0}) {
        int best = -1; double bd = 1e30;
        for (int n = 0; n < K0.mesh.node_count; ++n) {
            const double d = std::hypot(K0.mesh.x[n] - 0.5 * kR, K0.mesh.y[n] - zq);
            if (d < bd) { bd = d; best = n; }
        }
        if (best < 0) continue;
        const double z = K0.mesh.y[best];
        const double svfe = K0.stress.stress[best](1), svref = sigma_v_eff(z);
        const double pfe = best < (int)K0.pore.size() ? K0.pore[best] : 0.0;
        std::printf("   %5.2f  %12.4f  %12.4f  %10.4f  %10.4f\n", z, svfe, svref, pfe, pore_at(z));
        check(std::fabs(svfe - svref) < 0.005 * std::fabs(svref) + 0.05,
              "effective vertical stress = the buoyant closed form");
        check(std::fabs(pfe - pore_at(z)) < 1e-6 + 1e-6 * pore_at(z),
              "and the reported pore pressure is the hydrostatic one");
        ++checked;
    }
    check(checked == 3, "checked three depths");

    // --- the phase that can fail: lowering the table from kWt to kWtLow --------------------
    // Every point below the OLD table loses its buoyancy over the drained interval, and every
    // point between the two tables changes from saturated to moist weight as well. The net
    // effective-stress increase is (gamma_unsat - gamma_sat + gamma_w) per metre of lowering,
    // and a laterally confined column settles by its integral over the depth:
    //     s = (gamma_unsat - gamma_sat + gamma_w)/Eoed * [z_low (z_wt - z_low) + (z_wt - z_low)^2 / 2]
    const auto& D = res[1];
    check(D.ok, "the dewatering phase solved");
    if (!D.ok) { std::printf("   (%s)\n", D.message.c_str()); return; }
    const double eoed = kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));
    const double dgam = kGunsat - kGsat + gamma_w();
    const double drop = kWt - kWtLow;
    const double s_ref = dgam / eoed * (kWtLow * drop + 0.5 * drop * drop);
    double uz = 0.0, ur = 0.0;
    for (int n = 0; n < D.mesh.node_count; ++n) {
        uz = std::fmax(uz, std::fabs(D.disp(2 * n + 1)));
        ur = std::fmax(ur, std::fabs(D.disp(2 * n)));
    }
    std::printf("   dewatering %.1f m -> %.1f m: settlement %.6e m vs closed form %.6e m (%+.2f%%)\n",
                kWt, kWtLow, uz, s_ref, 100.0 * (uz / s_ref - 1.0));
    std::printf("   max |u_r| = %.3e m, %.2f%% of the settlement (the column is 1-D)\n",
                ur, 100.0 * ur / uz);
    check(std::fabs(uz / s_ref - 1.0) < 0.01,
          "the settlement of lowering the water table = the closed form within 1%");
    check(ur < 0.02 * uz,
          "and the ground barely moves radially -- which it would not do if the pore load were "
          "unbalanced in r");
}

}  // namespace

int main() {
    std::printf("Axisymmetric analysis with a water table\n\n");
    test_equilibrium_identity();
    test_through_the_driver();
    if (g_failures == 0) {
        std::printf("\nOK: saturated ground stands still in r-z, and the hoop pore term is why\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
