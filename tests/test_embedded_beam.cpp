// Embedded beam (pile row) -- skin interaction, step 1 verification (Faz A.4). A vertical pile
// (Timoshenko beam on its own DOFs) embedded in soil, coupled by skin springs (PLAXIS Sci.Man 7.5).
// With the soil held fixed, the skin reduces to a distributed axial spring k_a per unit length, so
// an axial head load P gives the classical "rod on an elastic foundation" head displacement:
//     u_head = (P / (EA lambda)) coth(lambda L),  lambda = sqrt(k_a / EA).
// This verifies the beam stiffness + the skin coupling assembly + point location of the pile in the
// (non-conforming) soil mesh. (Soil deformation coupling via N_s is exercised in a later step.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/structural_forces.hpp>   // embedded_beam_force_diagram (Output)
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/embedded_beam.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace ebeam = katai::core::ebeam;
namespace plate = katai::core::plate;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_axial_pile_fixed_soil() {
    // Soil mesh (only used to locate the pile; its DOFs are held fixed here).
    const RectangularDomain domain{0.0, 0.0, 5.0, 10.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 5, 10);

    // Vertical pile at x = 2.3 (off the mesh columns -> non-conforming, exercises N_s location).
    const double xw = 2.3, y0 = 0.5, y1 = 9.5, L = y1 - y0;
    const int ne = 18;                       // beam elements -> 2*ne+1 nodes
    const int nbn = 2 * ne + 1;
    std::vector<double> bx(nbn, xw), by(nbn);
    for (int i = 0; i < nbn; ++i) by[i] = y0 + L * i / (nbn - 1);
    std::vector<std::array<int, 3>> belem;
    for (int e = 0; e < ne; ++e) belem.push_back({2 * e, 2 * e + 2, 2 * e + 1});  // [A,B,mid]

    // Pile/skin properties.
    const double Ep = 3.0e7, d = 0.4;
    const double EA = Ep * d, EI = Ep * d * d * d / 12.0;
    const double k_axial = 1.0e4, k_lateral = 1.0e4;   // skin stiffness per unit length
    const double P = 100.0;
    const double lambda = std::sqrt(k_axial / EA);
    const double u_exact = (P / (EA * lambda)) * (1.0 / std::tanh(lambda * L));

    const auto skin = ebeam::build_skin_points(mesh, bx, by, belem, k_axial, k_lateral);
    int located = 0; for (const auto& sp : skin) located += sp.ok ? 1 : 0;
    std::printf("  skin points: %zu, located in soil: %d\n", skin.size(), located);
    check(located == static_cast<int>(skin.size()), "all skin points located in the soil mesh");

    // Beam DOF system: 3 DOF/node (ux, uy, phi). Soil fixed -> excluded.
    const int ndof = 3 * nbn;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
    plate::PlateProps pp; pp.EA = EA; pp.EI = EI; pp.nu = 0.0;
    // Beam (Timoshenko) stiffness.
    for (const auto& el : belem) {
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = bx[el[k]]; X(k, 1) = by[el[k]]; }
        const plate::ElementMatrix Ke = plate::stiffness(X, pp);
        int g[9];
        for (int k = 0; k < 3; ++k) { g[3 * k] = 3 * el[k]; g[3 * k + 1] = 3 * el[k] + 1; g[3 * k + 2] = 3 * el[k] + 2; }
        for (int a = 0; a < 9; ++a) for (int b = 0; b < 9; ++b) K(g[a], g[b]) += Ke(a, b);
    }
    // Skin (beam-beam block; soil fixed -> soil block drops): K_bb = wJ Nb_i Nb_j T (translations).
    for (const auto& sp : skin) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const Eigen::Matrix2d blk = sp.wJ * sp.Nb(i) * sp.Nb(j) * sp.T;
                const int bi = 3 * sp.beam_node[i], bj = 3 * sp.beam_node[j];
                K(bi, bj) += blk(0, 0); K(bi, bj + 1) += blk(0, 1);
                K(bi + 1, bj) += blk(1, 0); K(bi + 1, bj + 1) += blk(1, 1);
            }
    }

    // Axial head load (vertical pile -> axial = uy) at the top node.
    Eigen::VectorXd f = Eigen::VectorXd::Zero(ndof);
    const int head = nbn - 1;       // top node (y = y1)
    f(3 * head + 1) = -P;            // downward
    const Eigen::VectorXd u = K.ldlt().solve(f);
    const double u_head = std::fabs(u(3 * head + 1));

    std::printf("  axial pile: u_head=%.6e  closed-form=%.6e  lambda L=%.3f  (%.1f%%)\n",
                u_head, u_exact, lambda * L, 100.0 * (u_head - u_exact) / u_exact);
    check(std::fabs(u_head - u_exact) < 0.02 * u_exact,
          "head displacement = (P/EA lambda) coth(lambda L) (rod on elastic foundation)");
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// Solver-integrated embedded beam in a DEFORMING elastic soil: an axial head load is shed to the
// soil through the skin. Verifies the beam<->soil coupling via N_s (load transfer + equilibrium).
void test_axial_pile_in_soil() {
    const RectangularDomain domain{0.0, 0.0, 12.0, 12.0, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 12, 12);

    const double xw = 6.3, y0 = 1.0, y1 = 11.0;   // pile off the mesh columns
    const int ne = 10, nbn = 2 * ne + 1;
    std::vector<double> bx(nbn, xw), by(nbn);
    for (int i = 0; i < nbn; ++i) by[i] = y0 + (y1 - y0) * i / (nbn - 1);

    DofMap dofs(mesh.node_count, 2);
    const double Ep = 3.0e7, d = 0.4;
    plate::PlateProps pp; pp.EA = Ep * d; pp.EI = Ep * d * d * d / 12.0; pp.nu = 0.0;
    const double k_axial = 5.0e3, k_lateral = 5.0e3;
    auto beam = ebeam::build_embedded_beam(mesh, dofs, bx, by, pp, k_axial, k_lateral);

    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const double P = 100.0;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const int head = nbn - 1;
    f(dofs.equation(beam.dof_y[head])) = -P;

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, 2.0e4, 0.3}};
    katai::core::Structures st; st.embedded_beams = {beam};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_unsym, {1, 30, 1e-8},
                                                {}, {}, st);
    check(r.converged, "embedded beam in soil: solve converged");

    // Load transfer: total vertical skin force on the soil must equal the head load P (free toe,
    // no foot -> the pile sheds all P to the soil through the skin). Recompute from the solution.
    double skin_fy = 0.0;
    for (const auto& sp : beam.skin) {
        if (!sp.ok) continue;
        const auto Ns = katai::core::tri6::shape_functions(sp.xi_s, sp.eta_s);
        Eigen::Vector2d ub(0, 0), us(0, 0);
        for (int i = 0; i < 3; ++i) {
            ub(0) += sp.Nb(i) * r.displacement[beam.dof_x[sp.beam_node[i]]];
            ub(1) += sp.Nb(i) * r.displacement[beam.dof_y[sp.beam_node[i]]];
        }
        for (int j = 0; j < 6; ++j) {
            const int sn = mesh.node_of(sp.soil_elem, j);
            us(0) += Ns(j) * r.displacement[dofs.global_dof(sn, 0)];
            us(1) += Ns(j) * r.displacement[dofs.global_dof(sn, 1)];
        }
        const Eigen::Vector2d t = sp.T * (ub - us);
        skin_fy += sp.wJ * t(1);
    }
    const double u_head = std::fabs(r.displacement[beam.dof_y[head]]);
    std::printf("  pile in soil: u_head=%.4e  total skin Fy=%.3f (vs P=%.1f, %.1f%%)\n",
                u_head, skin_fy, P, 100.0 * (std::fabs(skin_fy) - P) / P);
    check(std::fabs(std::fabs(skin_fy) - P) < 0.02 * P,
          "total skin shear = head load P (load fully transferred to soil)");
    check(u_head > 0.0, "pile settles under axial load");

    // Internal-force diagram (PLAXIS Output -> embedded beam N/Q/M): under a pure axial head load the
    // axial force at the head equals the applied P (head equilibrium) and the bending moment ~ 0.
    const auto diag = katai::core::embedded_beam_force_diagram(beam, r.displacement);
    double Nmax = 0.0, Mmax = 0.0, N_at_head = 0.0, ytop = -1e30, N_at_toe = 0.0, ybot = 1e30;
    for (const auto& st : diag) {
        Nmax = std::fmax(Nmax, std::fabs(st.N));
        Mmax = std::fmax(Mmax, std::fabs(st.M));
        if (st.y > ytop) { ytop = st.y; N_at_head = std::fabs(st.N); }
        if (st.y < ybot) { ybot = st.y; N_at_toe = std::fabs(st.N); }
    }
    std::printf("  diagram: |N| head=%.3f (vs P=%.1f), |N| toe=%.3f, max|M|=%.4f\n",
                N_at_head, P, N_at_toe, Mmax);
    check(std::fabs(N_at_head - P) < 0.05 * P, "embedded beam diagram: axial force at the head = P");
    check(N_at_toe < 0.25 * P, "embedded beam diagram: axial force decays toward the (free) toe");
    check(Mmax < 0.05 * P, "embedded beam diagram: bending ~ 0 under pure axial load");
}

// Ultimate axial pile capacity: skin friction limit T_max (per length) + foot capacity F_max.
// Limit analysis (load control): a head load above the capacity cannot be equilibrated -> the pile
// plunges; the highest equilibrated load (load_factor * P_applied) is the ultimate capacity
//   Q_ult = T_max*L + F_max   (Q_skin + Q_base, classical pile capacity). Soil held fixed.
void test_pile_capacity() {
    const RectangularDomain domain{0.0, 0.0, 6.0, 11.0, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 6, 11);

    const double xw = 3.3, y0 = 0.5, y1 = 9.5, L = y1 - y0;
    const int ne = 12, nbn = 2 * ne + 1;
    std::vector<double> bx(nbn, xw), by(nbn);
    for (int i = 0; i < nbn; ++i) by[i] = y0 + L * i / (nbn - 1);

    DofMap dofs(mesh.node_count, 2);
    const double Ep = 3.0e7, d = 0.4;
    plate::PlateProps pp; pp.EA = Ep * d; pp.EI = Ep * d * d * d / 12.0; pp.nu = 0.0;
    const double k_axial = 1.0e4, k_lateral = 1.0e4;
    const double T_max = 50.0, D_foot = 1.0e5, F_max = 200.0;   // capacity inputs
    auto beam = ebeam::build_embedded_beam(mesh, dofs, bx, by, pp, k_axial, k_lateral,
                                           T_max, D_foot, F_max);
    const double Q_ult = T_max * L + F_max;

    for (int n = 0; n < mesh.node_count; ++n) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();

    // Apply a head load above capacity -> the limit load (load_factor*P) converges to Q_ult.
    const int head = nbn - 1;
    const double P = 1.3 * Q_ult;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    f(dofs.equation(beam.dof_y[head])) = -P;   // downward

    const std::vector<MaterialModel> mm = {{MaterialType::LinearElastic, 2.0e4, 0.3}};
    katai::core::Structures st; st.embedded_beams = {beam};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_unsym, {80, 50, 1e-7}, {}, {}, st);
    const double Q_limit = r.load_factor * P;
    std::printf("  pile capacity (limit analysis): Q_limit=%.1f  Q_ult=T_max*L+F_max=%.1f (%.1f%%)  lf=%.3f\n",
                Q_limit, Q_ult, 100.0 * (Q_limit - Q_ult) / Q_ult, r.load_factor);
    check(!r.converged, "load above capacity -> pile plunges (no full convergence)");
    check(std::fabs(Q_limit - Q_ult) < 0.03 * Q_ult, "ultimate axial capacity = T_max L + F_max");
}

} // namespace

int main() {
    test_axial_pile_fixed_soil();
    test_axial_pile_in_soil();
    test_pile_capacity();
    if (g_failures == 0) {
        std::printf("OK: embedded beam skin interaction verified (axial pile, closed form)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
