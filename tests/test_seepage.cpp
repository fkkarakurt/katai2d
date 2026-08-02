// Steady-state confined seepage (Darcy + continuity, ∇·(k∇h)=0). Two checks from
// docs/references/seepage-formulation.md §5:
//   1. 1D Darcy column: prescribed head h1,h2 at the two ends, impermeable sides
//      (natural q_n=0) -> LINEAR head profile h(x)=h1+(h2-h1)x/L (round-off) and a
//      UNIFORM specific discharge q_x = k_x(h1-h2)/L.
//   2. 2D confined patch test: a linear field h=a·x+b·y prescribed on the whole
//      boundary is reproduced exactly in the interior, for ANISOTROPIC k (∇²h=0).
// Exercises assemble_seepage (1-DOF-per-node DofMap, Dirichlet lift, SPD solve).
#include <katai/analysis/seepage.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::LinearElastic;
using katai::core::Permeability;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace tri6 = katai::core::tri6;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// Solve a confined seepage BVP: head fixed at `fixed` nodes to head_prescribed[n].
// Returns the full nodal head field (free DOFs solved, fixed overlaid).
Eigen::VectorXd solve_head(const Mesh& mesh, const std::vector<Permeability>& perm,
                           const std::vector<int>& fixed,
                           const std::vector<double>& head_prescribed) {
    DofMap dofs(mesh.node_count, 1);  // one DOF per node: the hydraulic head h
    for (int n : fixed) dofs.fix_node_component(n, 0);
    dofs.finalize();

    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_seepage(mesh, dofs, perm, head_prescribed, builder, rhs);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    Eigen::VectorXd h = katai::core::expand_to_full(dofs, solver->solve(rhs));
    for (int n = 0; n < mesh.node_count; ++n)
        if (dofs.is_fixed(dofs.global_dof(n, 0))) h[n] = head_prescribed[n];
    return h;
}

void test_darcy_column() {
    constexpr double L = 10.0, Hgt = 2.0;  // x in [0,L], y in [0,Hgt]
    constexpr double kx = 1e-4, ky = 1e-4;
    constexpr double h1 = 12.0, h2 = 4.0;  // head at left (x=0) / right (x=L)
    const RectangularDomain domain{0.0, 0.0, L, Hgt, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 10, 2);

    std::vector<int> fixed;
    std::vector<double> hp(mesh.node_count, 0.0);
    for (int n : mesh.left_nodes)  { fixed.push_back(n); hp[n] = h1; }
    for (int n : mesh.right_nodes) { fixed.push_back(n); hp[n] = h2; }

    const std::vector<Permeability> perm = {{kx, ky}};
    const Eigen::VectorXd h = solve_head(mesh, perm, fixed, hp);

    // Linear head profile h(x) = h1 + (h2 - h1) x / L (round-off).
    double max_err = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double exact = h1 + (h2 - h1) * mesh.x[n] / L;
        max_err = std::fmax(max_err, std::fabs(h[n] - exact));
    }

    // Specific discharge q = -k ∇h at Gauss points; expect uniform q_x=kx(h1-h2)/L,
    // q_y = 0. Recover ∇h via G = strain_displacement gradients.
    const double qx_exact = kx * (h1 - h2) / L;
    double max_qx_err = 0.0, max_qy = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri6::NodeCoords coords;
        Eigen::Matrix<double, 6, 1> he;
        for (int k = 0; k < 6; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            he(k) = h[n];
        }
        for (const auto& gp : tri6::gauss_points()) {
            const auto g = tri6::strain_displacement(coords, gp.xi, gp.eta);
            double dhdx = 0.0, dhdy = 0.0;
            for (int k = 0; k < 6; ++k) {
                dhdx += g.B(0, 2 * k) * he(k);
                dhdy += g.B(1, 2 * k + 1) * he(k);
            }
            max_qx_err = std::fmax(max_qx_err, std::fabs(-kx * dhdx - qx_exact));
            max_qy = std::fmax(max_qy, std::fabs(-ky * dhdy));
        }
    }
    std::printf("  darcy column: max|h err|=%.3e  max|qx err|=%.3e  max|qy|=%.3e\n",
                max_err, max_qx_err, max_qy);
    check(max_err < 1e-9 * (h1 - h2), "head profile linear h1+(h2-h1)x/L");
    check(max_qx_err < 1e-9 * std::fabs(qx_exact), "specific discharge qx uniform");
    check(max_qy < 1e-9 * std::fabs(qx_exact), "no transverse flux qy");
}

void test_confined_patch() {
    constexpr double W = 6.0, H = 4.0;
    constexpr double a = 0.7, b = -0.3;  // linear field h = a x + b y
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);

    // Prescribe h = a x + b y on the entire boundary; anisotropic k (kx != ky).
    std::vector<char> on_bdry(mesh.node_count, 0);
    for (int n : mesh.bottom_nodes) on_bdry[n] = 1;
    for (int n : mesh.top_nodes)    on_bdry[n] = 1;
    for (int n : mesh.left_nodes)   on_bdry[n] = 1;
    for (int n : mesh.right_nodes)  on_bdry[n] = 1;
    std::vector<int> fixed;
    std::vector<double> hp(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        if (on_bdry[n]) { fixed.push_back(n); hp[n] = a * mesh.x[n] + b * mesh.y[n]; }

    const std::vector<Permeability> perm = {{2.0e-4, 5.0e-4}};
    const Eigen::VectorXd h = solve_head(mesh, perm, fixed, hp);

    double max_err = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double exact = a * mesh.x[n] + b * mesh.y[n];
        max_err = std::fmax(max_err, std::fabs(h[n] - exact));
    }
    std::printf("  confined patch (anisotropic k): max|h err|=%.3e\n", max_err);
    check(max_err < 1e-9 * (std::fabs(a) * W + std::fabs(b) * H),
          "linear field reproduced in interior (patch test)");
}

// Build the all-boundary prescribed-head set for a structured mesh from an exact
// field, then solve and return the max nodal error |h - h_exact| over all nodes.
template <class Gen, class Field>
double mms_max_error(Gen generate, int nx, int ny, const Field& exact) {
    const RectangularDomain domain{0.0, 0.0, 2.0, 1.0, 0};
    const Mesh mesh = generate(domain, nx, ny);
    std::vector<char> bdry(mesh.node_count, 0);
    for (int n : mesh.bottom_nodes) bdry[n] = 1;
    for (int n : mesh.top_nodes)    bdry[n] = 1;
    for (int n : mesh.left_nodes)   bdry[n] = 1;
    for (int n : mesh.right_nodes)  bdry[n] = 1;
    std::vector<int> fixed;
    std::vector<double> hp(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        if (bdry[n]) { fixed.push_back(n); hp[n] = exact(mesh.x[n], mesh.y[n]); }

    const std::vector<Permeability> perm = {{1.0e-4, 1.0e-4}};  // isotropic: h harmonic
    const Eigen::VectorXd h = solve_head(mesh, perm, fixed, hp);
    double e = 0.0;
    for (int n = 0; n < mesh.node_count; ++n)
        e = std::fmax(e, std::fabs(h[n] - exact(mesh.x[n], mesh.y[n])));
    return e;
}

// Method of Manufactured Solutions (Roache 1998): h = e^x sin(y) is harmonic
// (∇²h = e^x sin y - e^x sin y = 0) -> an EXACT confined-seepage solution for
// isotropic k. It is TRANSCENDENTAL (not a polynomial), so NEITHER tri6 nor tri15
// reproduces it exactly: refining the mesh must drive the nodal error down at each
// element's theoretical rate (tri6/P2 ~O(h^3), tri15/P4 ~O(h^5)) -- a far stronger
// check than a polynomial field (a harmonic polynomial like x^3-3xy^2 is captured at
// the nodes to round-off by the Laplace Galerkin solver, so it does NOT exercise the
// approximation error).
void test_manufactured_harmonic() {
    auto exact = [](double x, double y) { return std::exp(x) * std::sin(y); };

    // tri6 (quadratic): nodal error must drop at ~O(h^3).
    const double a1 = mms_max_error(katai::mesh::generate_structured_tri6, 6, 3, exact);
    const double a2 = mms_max_error(katai::mesh::generate_structured_tri6, 12, 6, exact);
    const double a3 = mms_max_error(katai::mesh::generate_structured_tri6, 24, 12, exact);
    const double pt6 = std::log(a2 / a3) / std::log(2.0);
    std::printf("  MMS tri6  (e^x sin y): e(h)=%.3e e(h/2)=%.3e e(h/4)=%.3e  order=%.2f,%.2f\n",
                a1, a2, a3, std::log(a1 / a2) / std::log(2.0), pt6);
    check(a3 < a2 && a2 < a1, "tri6 error decreases monotonically under refinement");
    check(pt6 > 2.7, "tri6 nodal convergence order ~3 (P2 theory)");

    // tri15 (quartic): nodal error must drop at ~O(h^5) and be far smaller in absolute
    // terms. Coarse meshes keep the (very small) error well above round-off.
    const double b1 = mms_max_error(katai::mesh::generate_structured_tri15, 2, 1, exact);
    const double b2 = mms_max_error(katai::mesh::generate_structured_tri15, 4, 2, exact);
    const double b3 = mms_max_error(katai::mesh::generate_structured_tri15, 8, 4, exact);
    const double pt15 = std::log(b2 / b3) / std::log(2.0);
    std::printf("  MMS tri15 (e^x sin y): e(h)=%.3e e(h/2)=%.3e e(h/4)=%.3e  order=%.2f,%.2f\n",
                b1, b2, b3, std::log(b1 / b2) / std::log(2.0), pt15);
    check(b3 < b2 && b2 < b1, "tri15 error decreases monotonically under refinement");
    check(pt15 > 4.5, "tri15 nodal convergence order ~5 (P4 theory)");
    // tri15's higher order shows in absolute accuracy: 8x4 quartic (4.7e-7) beats the
    // 24x12 quadratic (1.4e-6) with far fewer elements.
    check(b3 < a3, "tri15 reaches higher accuracy with a much coarser mesh");
}

// Terzaghi upward-seepage column: the full coupling head -> pore -> effective stress.
// A saturated 1D column under steady UPWARD seepage of gradient i. The seepage solve
// gives the linear head h(y)=H+i(H-y); the head field drives the pore-pressure load;
// the deformation solve then returns the effective stress. Classic closed form (e.g.
// Das, Principles of Geotechnical Engineering; Craig's Soil Mechanics; Terzaghi 1943):
//     sigma'_v(y) = -(gamma' - i*gamma_w) (H - y),   gamma' = gamma_sat - gamma_w,
// i.e. the seepage force i*gamma_w per unit depth reduces the (buoyant) effective
// stress. At the CRITICAL gradient i_cr = gamma'/gamma_w the effective stress vanishes
// (quick / heave condition). Statically determinate -> exact (round-off) for any E,nu.
void test_terzaghi_upward_seepage(double i_grad, const char* label, double tol_scale) {
    constexpr double W = 2.0, H = 10.0;
    constexpr double gamma_sat = 20.0, gamma_w = 9.81;
    constexpr double gamma_eff = gamma_sat - gamma_w;  // buoyant unit weight gamma'
    constexpr double E = 1.0e4, nu = 0.3;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 10);

    // (1) Seepage solve -> head field h(y) = H + i (H - y): prescribe head at the two
    // ends (bottom = H(1+i), top = H), impermeable sides. Recovers the linear field.
    std::vector<int> fixed_h;
    std::vector<double> hp(mesh.node_count, 0.0);
    for (int n : mesh.bottom_nodes) { fixed_h.push_back(n); hp[n] = H + i_grad * H; }
    for (int n : mesh.top_nodes)    { fixed_h.push_back(n); hp[n] = H; }
    const std::vector<Permeability> perm = {{1e-4, 1e-4}};
    const Eigen::VectorXd head = solve_head(mesh, perm, fixed_h, hp);

    // (2) Deformation solve: drained stiffness, gravity(gamma_sat) + pore load built
    // from the head field. Confined column (bottom fixed, lateral rollers).
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    const std::vector<LinearElastic> materials = {{E, nu}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, materials, builder);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_gravity(mesh, dofs, {gamma_sat}, f);
    katai::core::assemble_pore_load_from_head(mesh, dofs, head, gamma_w, f);

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    const Eigen::VectorXd u = katai::core::expand_to_full(dofs, solver->solve(f));

    // (3) Gauss-point effective stress sigma' = D B u_e; compare sigma'_v to
    // -(gamma' - i gamma_w)(H - y).
    const double slope = gamma_eff - i_grad * gamma_w;  // effective unit weight w/ seepage
    const Eigen::Matrix3d D = materials[0].plane_strain_matrix();
    double max_err = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        tri6::NodeCoords coords;
        Eigen::Matrix<double, 12, 1> ue;
        for (int k = 0; k < 6; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            ue(2 * k) = u[dofs.global_dof(n, 0)];
            ue(2 * k + 1) = u[dofs.global_dof(n, 1)];
        }
        for (const auto& gp : tri6::gauss_points()) {
            const auto g = tri6::strain_displacement(coords, gp.xi, gp.eta);
            const Eigen::Vector3d sigma_eff = D * (g.B * ue);
            const tri6::ShapeValues N = tri6::shape_functions(gp.xi, gp.eta);
            double yq = 0.0;
            for (int k = 0; k < 6; ++k) yq += N(k) * coords(k, 1);
            max_err = std::fmax(max_err, std::fabs(sigma_eff(1) - (-slope * (H - yq))));
        }
    }
    std::printf("  terzaghi seepage (%s i=%.3f): slope=gamma'-i*gw=%.3f  max|sigma'_v err|=%.3e\n",
                label, i_grad, slope, max_err);
    check(max_err < tol_scale * gamma_eff * H, "effective sigma'_v matches -(gamma'-i*gw)(H-y)");
}

// Build a curved annular-sector mesh: generate a structured mesh on the logical
// rectangle (r in [r1,r2]) x (theta in [0,theta_max]) and remap each node to physical
// (x,y) = (r cos th, r sin th). tri6 mid-side nodes land on the arcs/radial lines
// (isoparametric quadratic edges). Logical left=r1 arc, right=r2 arc; theta-faces are
// radial streamlines (impermeable, natural).
Mesh make_annular_mesh(double r1, double r2, double theta_max, int nr, int nt) {
    // RectangularDomain is {x0, y0, WIDTH, HEIGHT}: logical r in [r1, r1+(r2-r1)].
    const RectangularDomain logical{r1, 0.0, r2 - r1, theta_max, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(logical, nr, nt);
    for (int n = 0; n < mesh.node_count; ++n) {
        const double r = mesh.x[n], th = mesh.y[n];
        mesh.x[n] = r * std::cos(th);
        mesh.y[n] = r * std::sin(th);
    }
    return mesh;
}

// Radial confined flow (Thiem): the first genuinely 2D curved-field check (prior tests
// are 1D or manufactured). Exact closed form h(r) = h1 + (h2-h1) ln(r/r1)/ln(r2/r1);
// the discharge through the sector is Q = k (h1-h2) theta / ln(r2/r1), independent of r
// (conservation). Validates the 2D solver AND the consistent nodal-flux discharge
// recovery + mass balance. Reference: Harr (1962); Bear (1972); steady well equation
// (Thiem 1906) -- the seepage analogue of the Lame cylinder (test_axisym_cylinder).
void test_radial_confined_flow() {
    const double r1 = 1.0, r2 = 3.0;
    const double theta_max = std::acos(-1.0) / 3.0;  // 60 degree sector
    const double h1 = 10.0, h2 = 4.0, k = 1.0e-4;
    const std::vector<Permeability> perm = {{k, k}};
    const double Q_exact = k * (h1 - h2) * theta_max / std::log(r2 / r1);

    auto run = [&](int nr, int nt, double& head_err, double& Q_in, double& Q_total) {
        Mesh mesh = make_annular_mesh(r1, r2, theta_max, nr, nt);
        std::vector<int> fixed;
        std::vector<double> hp(mesh.node_count, 0.0);
        for (int n : mesh.left_nodes)  { fixed.push_back(n); hp[n] = h1; }  // inner arc r1
        for (int n : mesh.right_nodes) { fixed.push_back(n); hp[n] = h2; }  // outer arc r2
        const Eigen::VectorXd head = solve_head(mesh, perm, fixed, hp);
        head_err = 0.0;
        for (int n = 0; n < mesh.node_count; ++n) {
            const double r = std::hypot(mesh.x[n], mesh.y[n]);
            const double exact = h1 + (h2 - h1) * std::log(r / r1) / std::log(r2 / r1);
            head_err = std::fmax(head_err, std::fabs(head[n] - exact));
        }
        const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, head);
        Q_in = 0.0;
        for (int n : mesh.left_nodes) Q_in += Q[n];  // discharge at the inner boundary
        Q_total = Q.sum();
    };

    double ec, qin_c, qt_c, ef, qin_f, qt_f;
    run(8, 4, ec, qin_c, qt_c);
    run(16, 8, ef, qin_f, qt_f);
    const double order = std::log(ec / ef) / std::log(2.0);
    std::printf("  radial flow: head e(coarse)=%.3e e(fine)=%.3e order=%.2f\n", ec, ef, order);
    std::printf("    Q_in=%.4e Q_exact=%.4e relerr=%.2e  sum(Q)=%.2e (mass balance)\n",
                std::fabs(qin_f), Q_exact, std::fabs(std::fabs(qin_f) - Q_exact) / Q_exact, qt_f);
    check(ef < ec && order > 1.8, "radial head converges to the log profile");
    check(ef < 5e-3 * (h1 - h2), "radial head accurate (<0.5% of head drop) on the fine mesh");
    check(std::fabs(qt_c) < 1e-9 * Q_exact && std::fabs(qt_f) < 1e-9 * Q_exact,
          "mass conservation: total nodal flux ~ 0 (round-off)");
    check(std::fabs(std::fabs(qin_f) - Q_exact) < 0.02 * Q_exact,
          "recovered discharge matches k(h1-h2)theta/ln(r2/r1) within 2%");
}

// Confined flow under a flat-bottom dam (the canonical PLAXIS confined-seepage
// problem). Along the impermeable base the head follows the EXACT arccos law for a
// homogeneous isotropic semi-infinite foundation (Harr 1962; Polubarinova-Kochina,
// Theory of Groundwater Movement; Das, Principles of Geotechnical Engineering):
//     h(x) = (H/pi) * arccos(x/b),   -b <= x <= b   (upstream H at -b, 0 at +b)
// The dam base is modelled as an impermeable segment on the TOP boundary (no internal
// slit). The exit gradient is SINGULAR at the dam edges x=+-b (the physical infinite
// exit gradient at the toe -- inherently mesh-sensitive for any FE code, PLAXIS
// included; engineering practice uses integral/averaged quantities). So validate the
// smooth central base against arccos (converging), with the edge region reported.
// Note: uplift force int_base h dx = H*b and h(0)=H/2 hold by ANTISYMMETRY (h(x)+
// h(-x)=H), so they are exact regardless of mesh -- the arccos SHAPE is the real test.
void test_dam_base_uplift() {
    const double b = 1.0, H = 10.0, X = 16.0, T = 16.0;  // domain large -> ~half-plane
    const double pi = std::acos(-1.0);
    const std::vector<Permeability> perm = {{1e-4, 1e-4}};

    // Solve on base resolution set by nx,ny; report max |h - arccos| over the central
    // base (|x| <= 0.5b) and over all interior base nodes; h(0); and the mass balance.
    auto run = [&](int nx, int ny, double& dev_central, double& dev_all, double& h0,
                   double& q_total) {
        const RectangularDomain dom{-X, -T, 2 * X, T, 0};  // top surface at y=0
        const Mesh mesh = katai::mesh::generate_structured_tri6(dom, nx, ny);
        std::vector<int> fixed;
        std::vector<double> hp(mesh.node_count, 0.0);
        for (int n : mesh.left_nodes)  { fixed.push_back(n); hp[n] = H; }
        for (int n : mesh.right_nodes) { fixed.push_back(n); hp[n] = 0.0; }
        const double eps = 1e-9;
        std::vector<int> base;
        for (int n : mesh.top_nodes) {
            const double x = mesh.x[n];
            if (x <= -b + eps)     { fixed.push_back(n); hp[n] = H; }
            else if (x >= b - eps) { fixed.push_back(n); hp[n] = 0.0; }
            else                   base.push_back(n);
        }
        const Eigen::VectorXd head = solve_head(mesh, perm, fixed, hp);
        dev_central = 0.0; dev_all = 0.0; h0 = 0.0;
        for (int n : base) {
            const double x = mesh.x[n];
            const double d = std::fabs(head[n] - (H / pi) * std::acos(x / b));
            dev_all = std::fmax(dev_all, d);
            if (std::fabs(x) <= 0.5 * b + eps) dev_central = std::fmax(dev_central, d);
            if (std::fabs(x) < eps) h0 = head[n];
        }
        q_total = katai::core::compute_nodal_flux(mesh, perm, head).sum();
    };

    double c1, a1, h0a, qt1, c2, a2, h0b, qt2;
    run(64, 32, c1, a1, h0a, qt1);    // base node spacing 0.25
    run(128, 64, c2, a2, h0b, qt2);   // base node spacing 0.125
    std::printf("  dam uplift (arccos law), central |x|<=0.5b / all-interior dev (%% of H):\n");
    std::printf("    coarse: central=%.2f%% all=%.2f%%\n", 100 * c1 / H, 100 * a1 / H);
    std::printf("    fine  : central=%.2f%% all=%.2f%%  h(0)=%.5f  sum(Q)=%.2e (mass balance)\n",
                100 * c2 / H, 100 * a2 / H, h0b, qt2);
    check(std::fabs(h0b - H / 2.0) < 1e-6 * H, "base center head = H/2 (antisymmetry, exact)");
    check(c2 < 0.01 * H, "central base head matches arccos within 1% of H");
    check(c2 < c1, "central arccos deviation decreases under mesh refinement");
    check(std::fabs(qt2) < 1e-10, "mass conservation: total nodal flux sum(Q) ~ 0 (round-off)");
}

// Neumann (prescribed flux) boundary: a 1D column with head h0 fixed at the left and a
// uniform inflow flux q_n on the right (top/bottom impermeable). Steady state -> uniform
// specific discharge, linear head h(x) = h0 + (q_n/k) x (reproduced to round-off). The
// recovered discharge at the left boundary equals the total inflow q_n * height (mass
// balance). Exercises assemble_seepage_flux.
void test_neumann_flux() {
    constexpr double L = 5.0, Hgt = 2.0, h0 = 3.0, k = 2.0e-4, q_n = 1.0e-4;
    const RectangularDomain dom{0.0, 0.0, L, Hgt, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(dom, 10, 4);

    DofMap dofs(mesh.node_count, 1);
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<double> hp(mesh.node_count, 0.0);
    for (int n : mesh.left_nodes) hp[n] = h0;

    const std::vector<Permeability> perm = {{k, k}};
    katai::math::SparseMatrixBuilder builder(dofs.equation_count());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_seepage(mesh, dofs, perm, hp, builder, rhs);
    katai::core::assemble_seepage_flux(mesh, dofs, mesh.right_nodes, q_n, rhs);  // inflow

    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(builder.build());
    Eigen::VectorXd head = katai::core::expand_to_full(dofs, solver->solve(rhs));
    for (int n : mesh.left_nodes) head[n] = h0;

    double max_err = 0.0;
    for (int n = 0; n < mesh.node_count; ++n)
        max_err = std::fmax(max_err, std::fabs(head[n] - (h0 + (q_n / k) * mesh.x[n])));

    const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, head);
    double q_left = 0.0;
    for (int n : mesh.left_nodes) q_left += Q[n];      // outflow at the fixed-head end
    const double q_total_in = q_n * Hgt;               // prescribed inflow at the right
    std::printf("  neumann flux: max|h - (h0+q/k x)|=%.3e  Q_left=%.4e  -q_n*H=%.4e  sum(Q)=%.2e\n",
                max_err, q_left, -q_total_in, Q.sum());
    check(max_err < 1e-9 * (q_n / k) * L, "head linear h0+(q_n/k)x under prescribed flux");
    check(std::fabs(q_left + q_total_in) < 1e-9 * q_total_in, "outflow balances prescribed inflow");
    check(std::fabs(Q.sum()) < 1e-12, "mass conservation with Neumann BC");
}

// Unconfined flow through a rectangular dam between two reservoirs (depths h1 > h2),
// the canonical free-surface benchmark. By CHARNY's theorem (1951; Polubarinova-
// Kochina; Bear) the discharge is EXACT regardless of the free-surface shape:
//     q = k (h1^2 - h2^2) / (2 L)   (Dupuit's formula, exact for the discharge)
// The phreatic surface drops from (0,h1) to (L,h2) (both submerged -> no seepage face).
// Exercises solve_unconfined_seepage (variable k_rel(psi) Picard iteration).
void test_unconfined_dam() {
    constexpr double L = 10.0, D = 6.0, h1 = 5.0, h2 = 1.0, k = 1.0e-4;
    const RectangularDomain dom{0.0, 0.0, L, D, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(dom, 40, 24);

    DofMap dofs(mesh.node_count, 1);
    std::vector<double> hp(mesh.node_count, 0.0);
    const double eps = 1e-9;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) { dofs.fix_node_component(n, 0); hp[n] = h1; }
    for (int n : mesh.right_nodes)
        if (mesh.y[n] <= h2 + eps) { dofs.fix_node_component(n, 0); hp[n] = h2; }
    dofs.finalize();

    const std::vector<Permeability> perm = {{k, k}};
    auto lin = [](const katai::math::CsrMatrix& A, const Eigen::VectorXd& b) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
        s->factorize(A);
        return s->solve(b);
    };
    katai::core::UnconfinedOptions opt;
    opt.transition = 0.05;  // ~ element size near the free surface; small -> accurate q
    const auto res = katai::core::solve_unconfined_seepage(mesh, dofs, perm, hp, lin, opt);

    const double q_exact = k * (h1 * h1 - h2 * h2) / (2.0 * L);
    const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, res.head, opt);
    double q_left = 0.0;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) q_left += Q[n];

    // Free-surface (phreatic) height at mid-span x=L/2: scan the node column for psi=h-y
    // sign change. Dupuit parabola y_s = sqrt(h1^2 - (h1^2-h2^2) x/L).
    const double xm = L / 2.0;
    double y_lo = 0.0, p_lo = -1e9, y_hi = D, p_hi = 1e9;
    for (int n = 0; n < mesh.node_count; ++n) {
        if (std::fabs(mesh.x[n] - xm) > 1e-6) continue;
        const double psi = res.head[n] - mesh.y[n];
        if (psi >= 0 && mesh.y[n] > y_lo) { y_lo = mesh.y[n]; p_lo = psi; }
        if (psi < 0 && mesh.y[n] < y_hi)  { y_hi = mesh.y[n]; p_hi = psi; }
    }
    const double y_surface = y_lo + (y_hi - y_lo) * p_lo / (p_lo - p_hi);  // psi=0 interp
    const double y_dupuit = std::sqrt(h1 * h1 - (h1 * h1 - h2 * h2) * xm / L);

    std::printf("  unconfined dam: iters=%d res=%.2e  q=%.4e q_charny=%.4e relerr=%.2e\n",
                res.iterations, res.residual, std::fabs(q_left), q_exact,
                std::fabs(std::fabs(q_left) - q_exact) / q_exact);
    std::printf("    free surface @x=L/2: y_fe=%.3f  y_dupuit=%.3f  (diff %.1f%%)\n",
                y_surface, y_dupuit, 100.0 * std::fabs(y_surface - y_dupuit) / y_dupuit);
    check(res.converged, "unconfined Picard iteration converged");
    check(std::fabs(std::fabs(q_left) - q_exact) < 0.005 * q_exact,
          "unconfined discharge matches Charny q=k(h1^2-h2^2)/(2L) within 0.5%");
    check(std::fabs(y_surface - y_dupuit) < 0.05 * y_dupuit,
          "phreatic surface near Dupuit parabola within 5%");
}

// Same unconfined dam, but with van Genuchten/Mualem k_rel -- the SAME unsaturated model the
// transient + fully-coupled flow uses (materials/water_retention.hpp) -- instead of the linear
// free-surface regularization. This is the steady <-> coupled CONSISTENCY merge (audit P2c), and it
// pins the suction SIGN MAP (suction = y - h): a flipped sign would drop k_rel < 1 in the SATURATED
// zone (h > y) and wreck the discharge. Charny's theorem still fixes q = k(h1^2-h2^2)/(2L) exactly,
// so recovering it certifies both the physics and the sign. The linear-model dam above and this one
// bracket the free-surface solver: retention = nullptr keeps the legacy linear ramp bit-for-bit.
void test_unconfined_dam_van_genuchten() {
    constexpr double L = 10.0, D = 6.0, h1 = 5.0, h2 = 1.0, k = 1.0e-4;
    const RectangularDomain dom{0.0, 0.0, L, D, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(dom, 40, 24);

    DofMap dofs(mesh.node_count, 1);
    std::vector<double> hp(mesh.node_count, 0.0);
    const double eps = 1e-9;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) { dofs.fix_node_component(n, 0); hp[n] = h1; }
    for (int n : mesh.right_nodes)
        if (mesh.y[n] <= h2 + eps) { dofs.fix_node_component(n, 0); hp[n] = h2; }
    dofs.finalize();

    const std::vector<Permeability> perm = {{k, k}};
    // van Genuchten "sand" (Carsel & Parrish 1988; PLAXIS default): a steep SWCC -> sharp phreatic
    // surface. {g_a, g_n, g_l, S_res, S_sat}.
    const std::vector<katai::core::WaterRetention> ret = {{14.5, 2.68, 0.5, 0.10, 1.0}};
    auto lin = [](const katai::math::CsrMatrix& A, const Eigen::VectorXd& b) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
        s->factorize(A);
        return s->solve(b);
    };
    katai::core::UnconfinedOptions opt;
    opt.retention = &ret;              // van Genuchten k_rel instead of the linear ramp
    opt.relax = 0.10; opt.max_iter = 1500; opt.tol = 1e-5;
    const auto res = katai::core::solve_unconfined_seepage(mesh, dofs, perm, hp, lin, opt);

    const double q_exact = k * (h1 * h1 - h2 * h2) / (2.0 * L);
    const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, res.head, opt);
    double q_left = 0.0;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) q_left += Q[n];

    // Sign-map sanity: the wet zone below the reservoir levels must stay saturated (h >= y). A
    // flipped suction would desaturate it -- check the deepest interior column.
    bool wet_saturated = true;
    for (int n = 0; n < mesh.node_count; ++n)
        if (mesh.y[n] < 0.5 && res.head[n] - mesh.y[n] < -1e-6) wet_saturated = false;

    const double relerr = std::fabs(std::fabs(q_left) - q_exact) / q_exact;
    std::printf("  unconfined dam (van Genuchten): iters=%d res=%.2e  q=%.4e q_charny=%.4e relerr=%.2e\n",
                res.iterations, res.residual, std::fabs(q_left), q_exact, relerr);
    check(res.converged, "unconfined Picard with van Genuchten k_rel converged");
    check(wet_saturated, "saturated zone stays wet (suction sign map correct)");
    check(relerr < 0.02, "van Genuchten discharge matches Charny within 2% (consistency merge)");
}

// Unconfined dam with a SEEPAGE FACE: upstream reservoir h1, downstream tailwater h2
// on the lower face and a free-draining seepage face above it (where water exits at
// atmospheric pressure, h=y, between the tailwater and the unknown exit point a0).
// By Charny's theorem the discharge is still EXACT: q = k(h1^2 - h2^2)/(2L). Exercises
// solve_unconfined_seepage_face (active-set on the downstream face + k_rel free surface).
void test_unconfined_seepage_face() {
    constexpr double L = 10.0, D = 6.0, h1 = 5.0, h2 = 1.0, k = 1.0e-4;
    const RectangularDomain dom{0.0, 0.0, L, D, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(dom, 40, 24);
    const double eps = 1e-9;

    std::vector<int> fixed_nodes, seepage_nodes;
    std::vector<double> fixed_values;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) { fixed_nodes.push_back(n); fixed_values.push_back(h1); }
    for (int n : mesh.right_nodes) {
        if (mesh.y[n] <= h2 + eps) { fixed_nodes.push_back(n); fixed_values.push_back(h2); }
        else seepage_nodes.push_back(n);  // free-draining face above the tailwater
    }

    const std::vector<Permeability> perm = {{k, k}};
    auto lin = [](const katai::math::CsrMatrix& A, const Eigen::VectorXd& b) {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
        s->factorize(A);
        return s->solve(b);
    };
    katai::core::UnconfinedOptions opt;
    opt.transition = 0.05; opt.relax = 0.05; opt.tol = 1e-5; opt.max_iter = 800;
    const auto res = katai::core::solve_unconfined_seepage_face(
        mesh, fixed_nodes, fixed_values, seepage_nodes, perm, lin, opt);

    const double q_exact = k * (h1 * h1 - h2 * h2) / (2.0 * L);
    const Eigen::VectorXd Q = katai::core::compute_nodal_flux(mesh, perm, res.head, opt);
    double q_in = 0.0;
    for (int n : mesh.left_nodes)
        if (mesh.y[n] <= h1 + eps) q_in += Q[n];

    // Exit point a0 = highest seepage node still discharging (psi ~ 0, h ~ y).
    double a0 = h2;
    for (int n : seepage_nodes)
        if (std::fabs(res.head[n] - mesh.y[n]) < 1e-3 && mesh.y[n] > a0) a0 = mesh.y[n];

    std::printf("  seepage face: iters=%d res=%.2e conv=%d  q=%.4e q_charny=%.4e relerr=%.2e  exit a0=%.2f\n",
                res.iterations, res.residual, (int)res.converged, std::fabs(q_in), q_exact,
                std::fabs(std::fabs(q_in) - q_exact) / q_exact, a0);
    check(res.converged, "seepage-face active-set + free-surface iteration converged");
    check(std::fabs(std::fabs(q_in) - q_exact) < 0.01 * q_exact,
          "discharge with seepage face matches Charny q=k(h1^2-h2^2)/(2L) within 1%");
    check(a0 > h2 + 0.1 && a0 < D, "exit point above tailwater and below crest (seepage face exists)");
}

} // namespace

int main() {
    test_darcy_column();
    test_confined_patch();
    test_manufactured_harmonic();
    // Coupling: hydrostatic (i=0 -> buoyant), upward seepage (reduced sigma'), and the
    // critical gradient i_cr = gamma'/gamma_w where the effective stress vanishes (heave).
    constexpr double i_cr = (20.0 - 9.81) / 9.81;  // gamma'/gamma_w
    test_terzaghi_upward_seepage(0.0, "hydrostatic", 1e-9);
    test_terzaghi_upward_seepage(0.3, "upward", 1e-9);
    test_terzaghi_upward_seepage(i_cr, "critical/heave", 1e-9);
    test_radial_confined_flow();
    test_dam_base_uplift();
    test_neumann_flux();
    test_unconfined_dam();
    test_unconfined_dam_van_genuchten();
    test_unconfined_seepage_face();
    if (g_failures == 0) {
        std::printf("OK: steady-state confined seepage verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
