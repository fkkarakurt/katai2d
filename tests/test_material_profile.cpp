// Depth-varying stiffness / cohesion: E(y) = E_ref + E_inc (y_ref - y), c(y) = c_ref + c_inc (y_ref - y)
// evaluated at each GAUSS (stress) point -- materials/material_model.hpp MaterialProfile.
//
// This is the commonest real soil profile (stiffness grows with depth), and in a seismic analysis E(y)
// IS the Vs profile. The material editor has offered E'_inc / y_ref all along, with help text spelling
// out the formula, and project_io persisted it -- but NO analysis implemented it, so every phase
// quietly solved the uniform E'_ref model. The seismic audit caught that and made it an honest refusal;
// this is the capability that replaces the refusal.
//
// Two oracles, neither of which shares code with the gradient path:
//
//   (a) UNIFORM LIMIT. E_inc = c_inc = 0 must reproduce the constant-E answer BIT-FOR-BIT. That is not
//       a formality: it proves the new per-Gauss integration collapses exactly onto the long-validated
//       element_stiffness path, so every existing verification still covers the uniform case.
//   (b) SUB-LAYER EQUIVALENCE. The same profile, approximated by N thin layers each with its own
//       CONSTANT E = E(y_mid), solved through the already-validated constant-E path. As N grows the two
//       must converge -- and the convergence itself is the evidence, because a sub-layer stack shares
//       nothing with the gradient code but the element routines. This is what a per-element (rather
//       than per-Gauss) gradient would fail on a coarse mesh.
//   (c) CLOSED FORM. A confined (oedometer) column under a surface load has eps = q/Eoed(y), so
//       u_top = q * integral(dy / Eoed(y)) -- an exact analytic settlement for a linear E(y).
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::core::MaterialModel;
using katai::core::MaterialProfile;
using katai::core::MaterialType;
using katai::geometry::RectangularDomain;
using katai::math::CsrMatrix;
using katai::math::SparseMatrixBuilder;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
Eigen::MatrixXd dense(const CsrMatrix& A) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows, A.cols);
    for (int r = 0; r < A.rows; ++r)
        for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p) D(r, A.col_indices[p]) = A.values[p];
    return D;
}
Eigen::VectorXd solve_dense(const CsrMatrix& A, const Eigen::VectorXd& b) {
    return dense(A).ldlt().solve(b);
}

constexpr double kH = 10.0, kW = 2.0;      // column: 2 m wide, 10 m tall
constexpr double kEref = 5000.0, kNu = 0.3;
constexpr double kEinc = 2000.0;           // E grows 2 MPa per metre of depth below y_ref = kH
constexpr double kQ = 100.0;               // surface pressure [kPa]

double E_at(double y) { return kEref + kEinc * (kH - y); }

// Oedometer (confined, ux = 0 on the sides): Eoed = E (1-nu) / ((1+nu)(1-2nu)).
double eoed(double E) { return E * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }

// A uniform vertical pressure q on the top of a laterally confined column: sigma_y = -q everywhere, so
// eps_y(y) = -q / Eoed(y) and the surface settlement is the exact integral
//     u_top = q * int_0^H dy / Eoed(y) = q (1+nu)(1-2nu)/((1-nu) Einc) * ln(E(0)/E(H)).
double settlement_closed_form() {
    const double k = (1.0 + kNu) * (1.0 - 2.0 * kNu) / (1.0 - kNu);
    return kQ * k / kEinc * std::log(E_at(0.0) / E_at(kH));
}

// Confined column, top loaded with q. `profile` empty -> uniform E_ref. `layers` > 1 -> the profile
// approximated by that many constant-E sub-layers (each element gets the material of its band).
struct Column { Mesh mesh; DofMap dofs; std::vector<LinearElastic> mats; std::vector<MaterialProfile> prof; };

// Build the column. mode: 0 = uniform E_ref, 1 = gradient via MaterialProfile, N>=2 = N sub-layers.
Column build(int ny, int mode) {
    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, ny);
    Column c{std::move(mesh), DofMap(0, 2), {}, {}};
    c.dofs = DofMap(c.mesh.node_count, 2);
    if (mode <= 1) {
        c.mats.push_back(LinearElastic{kEref, kNu});
        c.prof.push_back(MaterialProfile{});
        if (mode == 1) { c.prof[0].E_inc = kEinc; c.prof[0].y_ref = kH; }
    } else {
        // N constant-E sub-layers: each element takes the E of its band's mid-depth. Shares nothing
        // with the gradient code -- it is the validated constant-E path, N times.
        for (int i = 0; i < mode; ++i) {
            const double y_mid = kH * (i + 0.5) / mode;
            c.mats.push_back(LinearElastic{E_at(y_mid), kNu});
            c.prof.push_back(MaterialProfile{});
        }
        for (int e = 0; e < c.mesh.element_count; ++e) {
            double yc = 0.0;
            for (int k = 0; k < c.mesh.nodes_per_element; ++k) yc += c.mesh.y[c.mesh.node_of(e, k)];
            yc /= c.mesh.nodes_per_element;
            c.mesh.element_material[e] =
                std::min(mode - 1, std::max(0, (int)(yc / kH * mode)));
        }
    }
    for (int n : c.mesh.bottom_nodes) { c.dofs.fix_node_component(n, 0); c.dofs.fix_node_component(n, 1); }
    for (int n : c.mesh.left_nodes) c.dofs.fix_node_component(n, 0);    // confined -> oedometer
    for (int n : c.mesh.right_nodes) c.dofs.fix_node_component(n, 0);
    c.dofs.finalize();
    return c;
}

// Surface settlement under q, through assemble_stiffness (+ optional profile).
double solve_column(const Column& c) {
    SparseMatrixBuilder bK(c.dofs.equation_count());
    katai::core::assemble_stiffness(c.mesh, c.dofs, c.mats, bK, {}, c.prof);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(c.dofs.equation_count());
    katai::core::assemble_surface_traction(c.mesh, c.dofs, c.mesh.top_nodes, 0.0, -kQ, f);
    const Eigen::VectorXd u = solve_dense(bK.build(), f);
    double top = 0.0;
    for (int n : c.mesh.top_nodes) {
        const int eq = c.dofs.equation(c.dofs.global_dof(n, 1));
        if (eq >= 0) top = std::fmin(top, u[eq]);   // most-negative uy = settlement
    }
    return -top;
}

// (a) The uniform limit must be BIT-FOR-BIT the old path.
void test_uniform_limit() {
    std::printf("-- (a) uniform limit: E_inc = 0 reproduces the constant-E path exactly --\n");
    const Column c0 = build(10, 0);
    // Same model, but forced through the per-Gauss branch by a profile that is a no-op in disguise:
    // E_inc = 0 with a non-zero y_ref must still be recognised as uniform, and even a profile that
    // genuinely evaluates (c_inc only, on a linear-elastic material) must not touch the stiffness.
    Column c1 = build(10, 0);
    c1.prof[0].y_ref = 4.2;                     // irrelevant when both increments are zero
    check(c1.prof[0].uniform(), "a profile with E_inc = c_inc = 0 reports uniform() regardless of y_ref");

    SparseMatrixBuilder b0(c0.dofs.equation_count()), b1(c1.dofs.equation_count());
    katai::core::assemble_stiffness(c0.mesh, c0.dofs, c0.mats, b0);                 // no profile at all
    katai::core::assemble_stiffness(c1.mesh, c1.dofs, c1.mats, b1, {}, c1.prof);    // uniform profile
    const Eigen::MatrixXd K0 = dense(b0.build()), K1 = dense(b1.build());
    const double d = (K0 - K1).cwiseAbs().maxCoeff();
    std::printf("   max |K(no profile) - K(uniform profile)| = %.3e  (|K| = %.3e)\n",
                d, K0.cwiseAbs().maxCoeff());
    check(d == 0.0, "a uniform profile assembles the IDENTICAL stiffness matrix (bit-for-bit)");

    // And a profile that DOES evaluate per Gauss point, but with E_inc = 0, must give the same K to
    // round-off -- this exercises the new integration loop itself against the old element routine.
    Column c2 = build(10, 0);
    c2.prof[0].c_inc = 1.0;   // non-uniform() -> takes the per-Gauss branch; c does not affect elastic K
    check(!c2.prof[0].uniform(), "a c_inc-only profile is not uniform (it takes the per-Gauss branch)");
    SparseMatrixBuilder b2(c2.dofs.equation_count());
    katai::core::assemble_stiffness(c2.mesh, c2.dofs, c2.mats, b2, {}, c2.prof);
    const Eigen::MatrixXd K2 = dense(b2.build());
    const double d2 = (K0 - K2).cwiseAbs().maxCoeff();
    std::printf("   max |K(element routine) - K(per-Gauss loop, E_inc=0)| = %.3e (relative %.2e)\n",
                d2, d2 / K0.cwiseAbs().maxCoeff());
    check(d2 < 1e-9 * K0.cwiseAbs().maxCoeff(),
          "the per-Gauss integration reproduces E::element_stiffness when E is constant");
}

// (b) + (c) The gradient itself: closed form, and convergence of a sub-layer stack onto it.
void test_gradient() {
    std::printf("-- (b)/(c) gradient vs closed form and vs a constant-E sub-layer stack --\n");
    const double u_exact = settlement_closed_form();
    std::printf("   closed form: u_top = q (1+nu)(1-2nu)/((1-nu) E_inc) ln(E(0)/E(H)) = %.6f m\n", u_exact);

    const double u_grad = solve_column(build(20, 1));
    std::printf("   gradient (per-Gauss, ny=20):        u_top = %.6f m  (err %+.2f%%)\n",
                u_grad, 100 * (u_grad - u_exact) / u_exact);
    check(std::fabs(u_grad - u_exact) < 0.01 * u_exact,
          "gradient settlement = the exact integral q int dy/Eoed(y)");

    // The uniform model is FAR from it -- so the check above is not passing by accident.
    const double u_unif = solve_column(build(20, 0));
    std::printf("   uniform E_ref (what every phase used to solve): u_top = %.6f m  (err %+.0f%%)\n",
                u_unif, 100 * (u_unif - u_exact) / u_exact);
    check(u_unif > 1.5 * u_exact,
          "ignoring the gradient (the old silent behaviour) is grossly wrong -- the profile matters");

    // Sub-layer stack: an independent route to the same physics, through the validated constant-E path.
    double prev_err = 1e30;
    for (int n : {2, 4, 8, 16}) {
        const double u_n = solve_column(build(20, n));
        const double err = std::fabs(u_n - u_grad) / u_grad;
        std::printf("   %2d constant-E sub-layers:            u_top = %.6f m  (vs gradient %+.2f%%)\n",
                    n, u_n, 100 * (u_n - u_grad) / u_grad);
        if (n > 2) check(err < prev_err, "more sub-layers converge towards the gradient solution");
        prev_err = err;
    }
    check(prev_err < 0.01, "16 constant-E sub-layers agree with the per-Gauss gradient within 1%");
}

// The nonlinear solver must see the same profile -- and c(y) must reach the Mohr-Coulomb strength.
void test_nonlinear_path() {
    std::printf("-- gradient in the nonlinear solver (E and c per stress point) --\n");
    const Column c = build(20, 1);
    std::vector<MaterialModel> models(1);
    models[0].type = MaterialType::LinearElastic;
    models[0].youngs_modulus = kEref; models[0].poisson_ratio = kNu;

    Eigen::VectorXd f = Eigen::VectorXd::Zero(c.dofs.equation_count());
    katai::core::assemble_surface_traction(c.mesh, c.dofs, c.mesh.top_nodes, 0.0, -kQ, f);
    katai::core::NewtonOptions opt{1, 30, 1e-12};
    const auto nr = katai::core::solve_nonlinear(
        c.mesh, c.dofs, models, f,
        [](const CsrMatrix& A, const Eigen::VectorXd& b) { return solve_dense(A, b); },
        opt, {}, {}, {}, {}, {}, c.prof);
    check(nr.converged, "nonlinear solve with a gradient converged");
    double top = 0.0;
    for (int n : c.mesh.top_nodes) top = std::fmin(top, nr.displacement[c.dofs.global_dof(n, 1)]);
    const double u_nl = -top, u_exact = settlement_closed_form();
    std::printf("   nonlinear (LE) with gradient: u_top = %.6f m  (closed form %.6f, err %+.2f%%)\n",
                u_nl, u_exact, 100 * (u_nl - u_exact) / u_exact);
    // Same physics as the assembler route -> the two paths must agree with each other and with theory.
    check(std::fabs(u_nl - u_exact) < 0.01 * u_exact,
          "the nonlinear path applies the same E(y) as the linear assembler");

    // c(y) must reach the STRENGTH, not just the stiffness. Judge it by SETTLEMENT, not by a collapse
    // load: a laterally confined (oedometer) column CANNOT collapse in Tresca -- with eps_h = 0 the
    // stress path simply rides the yield surface (sigma_1 - sigma_3 = 2 su) while both stresses keep
    // growing, so there is no mechanism and the load factor is 1.0 either way. It does yield though
    // (in confined compression sigma_1 - sigma_3 = sigma_v (1 - nu/(1-nu)) = 0.571 sigma_v, so Tresca
    // bites at sigma_v ~ 17.5 kPa for su = 5), and the resulting plastic strain shows up as extra
    // settlement. A cohesion that grows with depth yields later and deeper -> LESS settlement.
    std::vector<MaterialModel> mc(1);
    mc[0].type = MaterialType::MohrCoulomb;
    mc[0].youngs_modulus = kEref; mc[0].poisson_ratio = kNu;
    mc[0].cohesion = 5.0; mc[0].friction_angle = 0.0; mc[0].dilatancy_angle = 0.0;   // Tresca, su = 5
    auto mc_settlement = [&](bool with_c_inc) {
        std::vector<MaterialProfile> p(1);
        p[0].E_inc = kEinc; p[0].y_ref = kH;                      // same stiffness profile in both
        if (with_c_inc) p[0].c_inc = 20.0;                        // su grows 20 kPa per metre of depth
        Eigen::VectorXd fq = Eigen::VectorXd::Zero(c.dofs.equation_count());
        katai::core::assemble_surface_traction(c.mesh, c.dofs, c.mesh.top_nodes, 0.0, -kQ, fq);
        katai::core::NewtonOptions o{10, 40, 1e-9};
        const auto r = katai::core::solve_nonlinear(
            c.mesh, c.dofs, mc, fq,
            [](const CsrMatrix& A, const Eigen::VectorXd& b) { return solve_dense(A, b); },
            o, {}, {}, {}, {}, {}, p);
        double t = 0.0;
        for (int n : c.mesh.top_nodes) t = std::fmin(t, r.displacement[c.dofs.global_dof(n, 1)]);
        return -t;
    };
    const double s_flat = mc_settlement(false), s_grad = mc_settlement(true);
    std::printf("   Tresca column under %.0f kPa: c = c_ref -> u_top = %.6f m | c(y) -> %.6f m (%+.1f%%)\n",
                kQ, s_flat, s_grad, 100 * (s_grad - s_flat) / s_flat);
    check(s_grad < s_flat, "a cohesion growing with depth yields less -> settles less (c_inc reaches the "
                           "Mohr-Coulomb strength; it is not ignored)");
    // And both must exceed the purely elastic answer -- otherwise nothing yielded and the check is void.
    check(s_flat > u_exact, "the flat-cohesion column really does yield (the comparison has teeth)");
}

}  // namespace

int main() {
    std::printf("Depth-varying stiffness / cohesion E(y), c(y) at the stress point\n\n");
    test_uniform_limit();
    std::printf("\n");
    test_gradient();
    std::printf("\n");
    test_nonlinear_path();
    if (g_failures == 0)
        std::printf("\nOK: uniform limit bit-for-bit, gradient = closed form, sub-layer stack converges,"
                    " nonlinear path agrees, c(y) reaches strength\n");
    else
        std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
