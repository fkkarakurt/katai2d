#include <katai/analysis/seepage.hpp>

#include <algorithm>
#include <array>
#include <cmath>

#include <Eigen/Dense>

#include <katai/fem/assembly/assembler.hpp>  // expand_to_full
#include <katai/fem/elements/element_traits.hpp>

namespace katai::core {
namespace {

// Unsaturated-zone relative permeability k_rel(ψ): ψ≥0 saturated → 1; ψ≤−transition →
// k_min; linear transition in between (smoothing, eases Picard convergence).
double relative_perm(double psi, double k_min, double transition) {
    if (psi >= 0.0) return 1.0;
    if (psi <= -transition) return k_min;
    return 1.0 + (-psi / transition) * (k_min - 1.0);
}

// 1D Lagrange edge shape functions N_k and derivatives dN_k/ds, for n equally spaced nodes
// on [-1,1] (n=3 quadratic tri6 edge, n=5 quartic tri15 edge). A copy for the seepage edge
// flux of the same helper in the surface traction (assembler.cpp).
void edge_shape(int n, double s, double* N, double* dN) {
    for (int k = 0; k < n; ++k) {
        const double sk = -1.0 + 2.0 * k / (n - 1);
        double val = 1.0, der = 0.0;
        for (int m = 0; m < n; ++m) {
            if (m == k) continue;
            const double sm = -1.0 + 2.0 * m / (n - 1);
            const double inv = 1.0 / (sk - sm);
            der = der * (s - sm) * inv + val * inv;
            val *= (s - sm) * inv;
        }
        N[k] = val;
        dN[k] = der;
    }
}

// Element conductivity He = ∫ Gᵀ diag(kx,ky) G detJ (n×n). The G rows come from the
// existing strain_displacement B: G(0,i)=dN_i/dx=B(0,2i), G(1,i)=dN_i/dy=B(1,2i+1).
template <class E>
Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> element_conductivity(
    const typename E::NodeCoords& coords, const Permeability& k) {
    Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> he =
        Eigen::Matrix<double, E::kNodeCount, E::kNodeCount>::Zero();
    for (const auto& gp : E::gauss_points()) {
        const auto g = E::strain_displacement(coords, gp.xi, gp.eta);
        Eigen::Matrix<double, 1, E::kNodeCount> gx, gy;
        for (int i = 0; i < E::kNodeCount; ++i) {
            gx(0, i) = g.B(0, 2 * i);
            gy(0, i) = g.B(1, 2 * i + 1);
        }
        const double w = gp.weight * g.det_jacobian;
        he.noalias() += w * (k.kx * gx.transpose() * gx + k.ky * gy.transpose() * gy);
    }
    return he;
}

template <class E>
void assemble_seepage_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<Permeability>& perm,
                           const std::vector<double>& head_prescribed,
                           math::SparseMatrixBuilder& builder,
                           Eigen::VectorXd& rhs) {
    typename E::NodeCoords coords;
    std::array<int, E::kNodeCount> nodes;

    for (int e = 0; e < mesh.element_count; ++e) {
        const Permeability& k = perm[mesh.element_material[e]];
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int n = mesh.node_of(e, i);
            nodes[i] = n;
            coords(i, 0) = mesh.x[n];
            coords(i, 1) = mesh.y[n];
        }
        const Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> he =
            element_conductivity<E>(coords, k);

        // Scatter: free-free → builder; free-fixed → Dirichlet lift (rhs).
        for (int a = 0; a < E::kNodeCount; ++a) {
            const int eq_a = dofs.equation(dofs.global_dof(nodes[a], 0));
            if (eq_a < 0) continue;
            for (int b = 0; b < E::kNodeCount; ++b) {
                const int gb = dofs.global_dof(nodes[b], 0);
                const int eq_b = dofs.equation(gb);
                if (eq_b >= 0)
                    builder.add_entry(eq_a, eq_b, he(a, b));
                else
                    rhs[eq_a] -= he(a, b) * head_prescribed[gb];
            }
        }
    }
}

template <class E>
void assemble_pore_from_head_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const Eigen::VectorXd& head, double gamma_w,
                                  Eigen::VectorXd& rhs,
                                  const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;       // mechanical 2-DOF
    std::array<double, E::kNodeCount> hnode;          // nodal head
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;  // NonPorous/passive: no pore load
        for (int k = 0; k < E::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            hnode[k] = head[n];
            element_dofs[2 * k] = dofs.global_dof(n, 0);
            element_dofs[2 * k + 1] = dofs.global_dof(n, 1);
        }
        for (const auto& gp : E::gauss_points()) {
            const auto g = E::strain_displacement(coords, gp.xi, gp.eta);
            const typename E::ShapeValues N = E::shape_functions(gp.xi, gp.eta);
            double h_gp = 0.0, y_gp = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) {
                h_gp += N(i) * hnode[i];
                y_gp += N(i) * coords(i, 1);
            }
            const double u = gamma_w * std::max(0.0, h_gp - y_gp);  // suction cut-off
            if (u == 0.0) continue;
            const double w = gp.weight * g.det_jacobian * u;
            // Bᵀ m (m=[1,1,0]): node i → [dN/dx, dN/dy] = [B(0,2i), B(1,2i+1)].
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int ex = dofs.equation(element_dofs[2 * i]);
                const int ey = dofs.equation(element_dofs[2 * i + 1]);
                if (ex >= 0) rhs[ex] += w * g.B(0, 2 * i);
                if (ey >= 0) rhs[ey] += w * g.B(1, 2 * i + 1);
            }
        }
    }
}

// Head-based saturation gravity: γ = (ψ = h−y ≥ 0) ? γ_sat : γ_unsat, with h interpolated
// by the SAME shape functions as the pore-load kernel (saturation/pore pressure form a
// consistent pair). The flow-solution (phreatic surface = ψ=0 contour) version of
// assemble_gravity_phreatic.
template <class E>
void assemble_gravity_from_head_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                     const std::vector<double>& gamma_unsat,
                                     const std::vector<double>& gamma_sat,
                                     const Eigen::VectorXd& head, Eigen::VectorXd& rhs) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    std::array<double, E::kNodeCount> hnode;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int mat = mesh.element_material[e];
        const double g_unsat = gamma_unsat[mat], g_sat = gamma_sat[mat];
        if (g_unsat == 0.0 && g_sat == 0.0) continue;
        for (int k = 0; k < E::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            hnode[k] = head[n];
            element_dofs[2 * k] = dofs.global_dof(n, 0);
            element_dofs[2 * k + 1] = dofs.global_dof(n, 1);
        }
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues N = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jac = dn.transpose() * coords;
            const double w_det = gp.weight * jac.determinant();
            double h_gp = 0.0, y_gp = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) {
                h_gp += N(i) * hnode[i];
                y_gp += N(i) * coords(i, 1);
            }
            const double gamma = (h_gp - y_gp >= 0.0) ? g_sat : g_unsat;
            if (gamma == 0.0) continue;
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int eq = dofs.equation(element_dofs[2 * i + 1]);  // uy
                if (eq >= 0) rhs[eq] += -gamma * w_det * N(i);
            }
        }
    }
}

template <class E>
Eigen::VectorXd nodal_flux_impl(const mesh::Mesh& mesh,
                                const std::vector<Permeability>& perm,
                                const Eigen::VectorXd& head) {
    Eigen::VectorXd flux = Eigen::VectorXd::Zero(mesh.node_count);
    typename E::NodeCoords coords;
    std::array<int, E::kNodeCount> nodes;
    Eigen::Matrix<double, E::kNodeCount, 1> h_e;
    for (int e = 0; e < mesh.element_count; ++e) {
        const Permeability& k = perm[mesh.element_material[e]];
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int n = mesh.node_of(e, i);
            nodes[i] = n;
            coords(i, 0) = mesh.x[n];
            coords(i, 1) = mesh.y[n];
            h_e(i) = head[n];
        }
        // Q_e = He·h_e; scatter to nodes. (Sum ≈ 0 at free nodes; discharge at Dirichlet.)
        const Eigen::Matrix<double, E::kNodeCount, 1> qe =
            element_conductivity<E>(coords, k) * h_e;
        for (int i = 0; i < E::kNodeCount; ++i) flux(nodes[i]) += qe(i);
    }
    return flux;
}

// Relative-permeability-weighted element conductivity He = ∫ Gᵀ k·k_rel(ψ_gp) G detJ;
// ψ_gp = h_gp − y_gp (interpolated from the current head). The shared kernel of the
// unconfined assembly/flux. If `wr` is given, k_rel = van Genuchten/Mualem (suction =
// y_gp − h_gp; consistent with transient/coupled flow); if nullptr, the simple linear ramp
// relative_perm(k_min, transition).
template <class E>
Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> element_conductivity_krel(
    const typename E::NodeCoords& coords,
    const Eigen::Matrix<double, E::kNodeCount, 1>& h_e, const Permeability& k,
    double k_min, double transition, const WaterRetention* wr) {
    Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> he =
        Eigen::Matrix<double, E::kNodeCount, E::kNodeCount>::Zero();
    for (const auto& gp : E::gauss_points()) {
        const auto g = E::strain_displacement(coords, gp.xi, gp.eta);
        const typename E::ShapeValues N = E::shape_functions(gp.xi, gp.eta);
        double h_gp = 0.0, y_gp = 0.0;
        for (int i = 0; i < E::kNodeCount; ++i) {
            h_gp += N(i) * h_e(i);
            y_gp += N(i) * coords(i, 1);
        }
        const double kr = wr ? relative_permeability_psi(*wr, y_gp - h_gp)
                             : relative_perm(h_gp - y_gp, k_min, transition);
        Eigen::Matrix<double, 1, E::kNodeCount> gx, gy;
        for (int i = 0; i < E::kNodeCount; ++i) {
            gx(0, i) = g.B(0, 2 * i);
            gy(0, i) = g.B(1, 2 * i + 1);
        }
        const double w = gp.weight * g.det_jacobian * kr;
        he.noalias() += w * (k.kx * gx.transpose() * gx + k.ky * gy.transpose() * gy);
    }
    return he;
}

template <class E>
void assemble_krel_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                        const std::vector<Permeability>& perm,
                        const Eigen::VectorXd& head,
                        const std::vector<double>& head_prescribed, double k_min,
                        double transition, const std::vector<WaterRetention>* retention,
                        math::SparseMatrixBuilder& builder, Eigen::VectorXd& rhs) {
    typename E::NodeCoords coords;
    std::array<int, E::kNodeCount> nodes;
    Eigen::Matrix<double, E::kNodeCount, 1> h_e;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int mat = mesh.element_material[e];
        const Permeability& k = perm[mat];
        const WaterRetention* wr = retention ? &(*retention)[mat] : nullptr;
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int n = mesh.node_of(e, i);
            nodes[i] = n;
            coords(i, 0) = mesh.x[n];
            coords(i, 1) = mesh.y[n];
            h_e(i) = head[n];
        }
        const Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> he =
            element_conductivity_krel<E>(coords, h_e, k, k_min, transition, wr);
        for (int a = 0; a < E::kNodeCount; ++a) {
            const int eq_a = dofs.equation(dofs.global_dof(nodes[a], 0));
            if (eq_a < 0) continue;
            for (int b = 0; b < E::kNodeCount; ++b) {
                const int gb = dofs.global_dof(nodes[b], 0);
                const int eq_b = dofs.equation(gb);
                if (eq_b >= 0)
                    builder.add_entry(eq_a, eq_b, he(a, b));
                else
                    rhs[eq_a] -= he(a, b) * head_prescribed[gb];
            }
        }
    }
}

template <class E>
Eigen::VectorXd nodal_flux_krel_impl(const mesh::Mesh& mesh,
                                     const std::vector<Permeability>& perm,
                                     const Eigen::VectorXd& head, double k_min,
                                     double transition,
                                     const std::vector<WaterRetention>* retention) {
    Eigen::VectorXd flux = Eigen::VectorXd::Zero(mesh.node_count);
    typename E::NodeCoords coords;
    std::array<int, E::kNodeCount> nodes;
    Eigen::Matrix<double, E::kNodeCount, 1> h_e;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int mat = mesh.element_material[e];
        const Permeability& k = perm[mat];
        const WaterRetention* wr = retention ? &(*retention)[mat] : nullptr;
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int n = mesh.node_of(e, i);
            nodes[i] = n;
            coords(i, 0) = mesh.x[n];
            coords(i, 1) = mesh.y[n];
            h_e(i) = head[n];
        }
        const Eigen::Matrix<double, E::kNodeCount, 1> qe =
            element_conductivity_krel<E>(coords, h_e, k, k_min, transition, wr) * h_e;
        for (int i = 0; i < E::kNodeCount; ++i) flux(nodes[i]) += qe(i);
    }
    return flux;
}

// Expands the 1-DOF free solution to the full nodal field + overlays the Dirichlet nodes.
Eigen::VectorXd expand_overlay(const mesh::Mesh& mesh, const DofMap& dofs,
                               const Eigen::VectorXd& free_head,
                               const std::vector<double>& prescribed) {
    Eigen::VectorXd h = expand_to_full(dofs, free_head);
    for (int n = 0; n < mesh.node_count; ++n)
        if (dofs.is_fixed(dofs.global_dof(n, 0))) h[n] = prescribed[n];
    return h;
}

} // namespace

void assemble_seepage(const mesh::Mesh& mesh, const DofMap& dofs,
                      const std::vector<Permeability>& permeability,
                      const std::vector<double>& head_prescribed,
                      math::SparseMatrixBuilder& builder, Eigen::VectorXd& rhs) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_seepage_impl<Tri15Element>(mesh, dofs, permeability,
                                            head_prescribed, builder, rhs);
    else
        assemble_seepage_impl<Tri6Element>(mesh, dofs, permeability,
                                           head_prescribed, builder, rhs);
}

void assemble_pore_load_from_head(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const Eigen::VectorXd& head, double gamma_w,
                                  Eigen::VectorXd& rhs,
                                  const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_pore_from_head_impl<Tri15Element>(mesh, dofs, head, gamma_w, rhs, active_element);
    else
        assemble_pore_from_head_impl<Tri6Element>(mesh, dofs, head, gamma_w, rhs, active_element);
}

void assemble_gravity_from_head(const mesh::Mesh& mesh, const DofMap& dofs,
                                const std::vector<double>& gamma_unsat,
                                const std::vector<double>& gamma_sat,
                                const Eigen::VectorXd& head, Eigen::VectorXd& rhs) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_gravity_from_head_impl<Tri15Element>(mesh, dofs, gamma_unsat, gamma_sat, head, rhs);
    else
        assemble_gravity_from_head_impl<Tri6Element>(mesh, dofs, gamma_unsat, gamma_sat, head, rhs);
}

Eigen::VectorXd compute_nodal_flux(const mesh::Mesh& mesh,
                                   const std::vector<Permeability>& permeability,
                                   const Eigen::VectorXd& head) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return nodal_flux_impl<Tri15Element>(mesh, permeability, head);
    return nodal_flux_impl<Tri6Element>(mesh, permeability, head);
}

void assemble_seepage_flux(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<int>& ordered_boundary_nodes,
                           double q_n, Eigen::VectorXd& rhs) {
    // Nodes per edge: tri6 → 3 (quadratic), tri15 → 5 (quartic). 4-point 1D Gauss (degree 7
    // — more than exact for the quartic edge). Same rule as assemble_surface_traction.
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    if (chain < npe || (chain - 1) % (npe - 1) != 0) return;  // incompatible chain
    const int edge_count = (chain - 1) / (npe - 1);

    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        std::array<double, 5> integral_n{};  // ∫ N_i ds
        for (int q = 0; q < 4; ++q) {
            edge_shape(npe, point[q], N.data(), dN.data());
            double dxdt = 0.0, dydt = 0.0;
            for (int i = 0; i < npe; ++i) {
                const int ni = ordered_boundary_nodes[base + i];
                dxdt += dN[i] * mesh.x[ni];
                dydt += dN[i] * mesh.y[ni];
            }
            const double ds = std::sqrt(dxdt * dxdt + dydt * dydt);
            for (int i = 0; i < npe; ++i) integral_n[i] += weight[q] * ds * N[i];
        }
        for (int i = 0; i < npe; ++i) {
            const int ni = ordered_boundary_nodes[base + i];
            const int eq = dofs.equation(dofs.global_dof(ni, 0));
            if (eq >= 0) rhs[eq] += q_n * integral_n[i];  // inflow positive
        }
    }
}

UnconfinedResult solve_unconfined_seepage(
    const mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<Permeability>& permeability,
    const std::vector<double>& head_prescribed,
    const SeepageLinearSolve& linear_solve, const UnconfinedOptions& options) {
    const int neq = dofs.equation_count();
    const bool tri15 = mesh.nodes_per_element == Tri15Element::kNodeCount;

    // Assemble with k_rel(head) + solve → full (overlaid) nodal head field.
    auto step = [&](const Eigen::VectorXd& current) {
        math::SparseMatrixBuilder builder(neq);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(neq);
        if (tri15)
            assemble_krel_impl<Tri15Element>(mesh, dofs, permeability, current,
                                             head_prescribed, options.k_min,
                                             options.transition, options.retention, builder, rhs);
        else
            assemble_krel_impl<Tri6Element>(mesh, dofs, permeability, current,
                                            head_prescribed, options.k_min,
                                            options.transition, options.retention, builder, rhs);
        return expand_overlay(mesh, dofs, linear_solve(builder.build(), rhs),
                              head_prescribed);
    };

    // Start: the fully saturated (confined) solution — a good initial guess.
    Eigen::VectorXd head;
    {
        math::SparseMatrixBuilder builder(neq);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(neq);
        assemble_seepage(mesh, dofs, permeability, head_prescribed, builder, rhs);
        head = expand_overlay(mesh, dofs, linear_solve(builder.build(), rhs),
                              head_prescribed);
    }

    UnconfinedResult res;
    res.head = head;
    for (int it = 1; it <= options.max_iter; ++it) {
        const Eigen::VectorXd hn = step(head);
        const Eigen::VectorXd prev = head;
        head = (1.0 - options.relax) * head + options.relax * hn;
        for (int n = 0; n < mesh.node_count; ++n)
            if (dofs.is_fixed(dofs.global_dof(n, 0))) head[n] = head_prescribed[n];
        // Convergence = the step-to-step change of the RELAXED head (has the accumulator
        // settled), NOT the oscillation of the raw hn. ONLY in the saturated zone (ψ≥0): in
        // the unsaturated (dry) zone k_eff=k·k_min is tiny → the head is ill-conditioned/
        // oscillates meaninglessly without affecting the physical solution (phreatic surface
        // + discharge); measuring over all nodes gives spurious divergence.
        double dh = 0.0, scale = 1e-30;
        for (int n = 0; n < mesh.node_count; ++n) {
            if (head[n] - mesh.y[n] < 0.0) continue;
            dh = std::fmax(dh, std::fabs(head[n] - prev[n]));
            scale = std::fmax(scale, std::fabs(head[n]));
        }
        res.iterations = it;
        res.residual = dh / scale;
        if (dh < options.tol * scale) {
            res.converged = true;
            break;
        }
    }
    res.head = head;
    return res;
}

Eigen::VectorXd compute_nodal_flux(const mesh::Mesh& mesh,
                                   const std::vector<Permeability>& permeability,
                                   const Eigen::VectorXd& head,
                                   const UnconfinedOptions& options) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return nodal_flux_krel_impl<Tri15Element>(mesh, permeability, head, options.k_min,
                                                  options.transition, options.retention);
    return nodal_flux_krel_impl<Tri6Element>(mesh, permeability, head, options.k_min,
                                             options.transition, options.retention);
}

UnconfinedResult solve_unconfined_seepage_face(
    const mesh::Mesh& mesh, const std::vector<int>& fixed_nodes,
    const std::vector<double>& fixed_values, const std::vector<int>& seepage_nodes,
    const std::vector<Permeability>& permeability,
    const SeepageLinearSolve& linear_solve, const UnconfinedOptions& options) {
    const int nc = mesh.node_count;
    const bool tri15 = mesh.nodes_per_element == Tri15Element::kNodeCount;

    std::vector<char> base_fixed(nc, 0);
    std::vector<double> base_val(nc, 0.0);
    for (std::size_t i = 0; i < fixed_nodes.size(); ++i) {
        base_fixed[fixed_nodes[i]] = 1;
        base_val[fixed_nodes[i]] = fixed_values[i];
    }
    std::vector<char> exit(nc, 0);  // active discharge (h=y) state

    // One Picard step with the given active set: build the DofMap (base ∪ active exits),
    // k_rel assemble, solve → full (overlaid) head.
    auto solve_with_set = [&](const Eigen::VectorXd& current) {
        DofMap dofs(nc, 1);
        std::vector<double> hp = base_val;
        for (int n = 0; n < nc; ++n)
            if (base_fixed[n]) dofs.fix_node_component(n, 0);
        for (int n : seepage_nodes)
            if (exit[n]) { dofs.fix_node_component(n, 0); hp[n] = mesh.y[n]; }
        dofs.finalize();
        const int neq = dofs.equation_count();
        math::SparseMatrixBuilder builder(neq);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(neq);
        if (tri15)
            assemble_krel_impl<Tri15Element>(mesh, dofs, permeability, current, hp,
                                             options.k_min, options.transition, options.retention,
                                             builder, rhs);
        else
            assemble_krel_impl<Tri6Element>(mesh, dofs, permeability, current, hp,
                                            options.k_min, options.transition, options.retention,
                                            builder, rhs);
        return expand_overlay(mesh, dofs, linear_solve(builder.build(), rhs), hp);
    };

    Eigen::VectorXd head = Eigen::VectorXd::Zero(nc);
    for (int n = 0; n < nc; ++n) if (base_fixed[n]) head[n] = base_val[n];
    head = solve_with_set(head);  // start (all seepage nodes free)

    UnconfinedResult res;
    res.head = head;
    for (int it = 1; it <= options.max_iter; ++it) {
        const Eigen::VectorXd hn = solve_with_set(head);
        const Eigen::VectorXd prev = head;
        head = (1.0 - options.relax) * head + options.relax * hn;
        for (int n = 0; n < nc; ++n) {
            if (base_fixed[n]) head[n] = base_val[n];
            else if (exit[n]) head[n] = mesh.y[n];
        }
        // Convergence = the step-to-step change of the RELAXED head (has the accumulator
        // settled), NOT the oscillation of the raw hn: even if a free-surface node
        // limit-cycles in the raw map, the head converges. Saturated zone only (ψ≥0).
        double dh = 0.0, scale = 1e-30;
        for (int n = 0; n < nc; ++n) {
            if (head[n] - mesh.y[n] < 0.0) continue;
            dh = std::fmax(dh, std::fabs(head[n] - prev[n]));
            scale = std::fmax(scale, std::fabs(head[n]));
        }
        // Active-set update: a free face node with ψ>0 → discharging; a discharging node
        // with inward flux (Q>0) → freed (outflow Q<0 = true discharge, see the Neumann
        // sign convention).
        const Eigen::VectorXd Q = (mesh.nodes_per_element == Tri15Element::kNodeCount)
            ? nodal_flux_krel_impl<Tri15Element>(mesh, permeability, head, options.k_min,
                                                 options.transition, options.retention)
            : nodal_flux_krel_impl<Tri6Element>(mesh, permeability, head, options.k_min,
                                                options.transition, options.retention);
        // Dead band (hysteresis): free→exit only if ψ is clearly positive (>transition);
        // exit→free only if the flux is clearly inward. A node near the exit point is
        // ambiguous at ψ≈0 (mesh-resolution limit) — the dead band keeps it in its current
        // state and cuts the exit↔no-flow chatter (spurious divergence).
        bool set_changed = false;
        for (int n : seepage_nodes) {
            if (!exit[n]) {
                if (head[n] - mesh.y[n] > options.transition) { exit[n] = 1; set_changed = true; }
            } else if (Q[n] > 0.0) {  // inflow → not a true discharge
                exit[n] = 0; set_changed = true;
            }
        }
        res.iterations = it;
        res.residual = dh / scale;
        if (!set_changed && dh < options.tol * scale) {
            res.converged = true;
            break;
        }
    }
    res.head = head;
    return res;
}

} // namespace katai::core
