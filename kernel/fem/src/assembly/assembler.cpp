#include <katai/fem/assembly/assembler.hpp>

#include <array>
#include <cassert>
#include <cmath>

#include <Eigen/Dense>

#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>

namespace katai::core {
namespace {

// Node coordinates of element e + global indices of its local DOFs.
template <class E>
void gather_element(const mesh::Mesh& mesh, const DofMap& dofs, int e,
                    typename E::NodeCoords& coords,
                    std::array<int, E::kDofCount>& element_dofs) {
    for (int k = 0; k < E::kNodeCount; ++k) {
        const int n = mesh.node_of(e, k);
        coords(k, 0) = mesh.x[n];
        coords(k, 1) = mesh.y[n];
        element_dofs[2 * k] = dofs.global_dof(n, 0);
        element_dofs[2 * k + 1] = dofs.global_dof(n, 1);
    }
}

// An element's profiled (E varying with depth) stiffness: Ke = Σ_g w_g·J_g·Bᵀ·D(y_g)·B.
// The element API is UNCHANGED — strain_displacement / gauss_points / shape_functions are
// already public. Uses the SAME rule as E::element_stiffness; the only difference is that
// D is evaluated at the Gauss point (gives the exact same result for a uniform profile;
// test_material_profile pins this to machine precision).
// `gauss` non-null -> D comes straight from that per-Gauss array (stress-dependent stiffness, e.g. HS
// E_ur(sigma3) computed by the caller from the committed state); otherwise from the depth profile.
template <class E>
typename E::ElementMatrix element_stiffness_profiled(const typename E::NodeCoords& coords,
                                                     const LinearElastic& m, const MaterialProfile& p,
                                                     const LinearElastic* gauss) {
    typename E::ElementMatrix ke = E::ElementMatrix::Zero();
    int g = 0;
    for (const auto& gp : E::gauss_points()) {
        const auto sg = E::strain_displacement(coords, gp.xi, gp.eta);
        LinearElastic mg = m;
        if (gauss) {
            mg = gauss[g];
        } else {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            double y = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) y += n(i) * coords(i, 1);
            mg.youngs_modulus = profile_at(m.youngs_modulus, p.E_inc, p.y_ref, y);
        }
        ke += (gp.weight * sg.det_jacobian) * (sg.B.transpose() * mg.plane_strain_matrix() * sg.B);
        ++g;
    }
    return ke;
}

template <class E>
void assemble_stiffness_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<LinearElastic>& materials,
                             math::SparseMatrixBuilder& builder,
                             const std::vector<char>& active_element,
                             const std::vector<MaterialProfile>& profile,
                             const std::vector<LinearElastic>& gauss_elastic) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    const bool per_gauss =
        gauss_elastic.size() == (std::size_t)mesh.element_count * E::kGaussCount;

    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;   // staged: excavated / not yet placed
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        const int mat = mesh.element_material[e];
        const LinearElastic& material = materials[mat];
        const MaterialProfile prof = mat < (int)profile.size() ? profile[mat] : MaterialProfile{};
        const LinearElastic* gp_e = per_gauss ? &gauss_elastic[(std::size_t)e * E::kGaussCount] : nullptr;
        const typename E::ElementMatrix ke =
            (prof.uniform() && !gp_e) ? E::element_stiffness(coords, material.plane_strain_matrix())
                                      : element_stiffness_profiled<E>(coords, material, prof, gp_e);

        // Scatter the free-free contributions to equation indices; drop fixed DOFs.
        for (int a = 0; a < E::kDofCount; ++a) {
            const int eq_a = dofs.equation(element_dofs[a]);
            if (eq_a < 0) continue;
            for (int b = 0; b < E::kDofCount; ++b) {
                const int eq_b = dofs.equation(element_dofs[b]);
                if (eq_b < 0) continue;
                builder.add_entry(eq_a, eq_b, ke(a, b));
            }
        }
    }
}

template <class E>
void assemble_mass_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                        const std::vector<double>& density,
                        math::SparseMatrixBuilder& builder,
                        const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;

    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;   // staged: excavated / not yet placed
        const double rho = density[mesh.element_material[e]];
        if (rho == 0.0) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);

        // Scalar consistent mass m_ij = ρ ∫ N_i N_j dA (Σw_g detJ_g N_i N_j).
        Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> m =
            Eigen::Matrix<double, E::kNodeCount, E::kNodeCount>::Zero();
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jacobian = dn.transpose() * coords;
            const double w_det = gp.weight * jacobian.determinant();
            for (int i = 0; i < E::kNodeCount; ++i)
                for (int j = 0; j < E::kNodeCount; ++j)
                    m(i, j) += rho * w_det * n(i) * n(j);
        }

        // Scatter to the x and y DOF blocks (block-diagonal; no x-y coupling). Drop fixed DOFs.
        for (int i = 0; i < E::kNodeCount; ++i)
            for (int j = 0; j < E::kNodeCount; ++j)
                for (int d = 0; d < 2; ++d) {
                    const int eq_a = dofs.equation(element_dofs[2 * i + d]);
                    const int eq_b = dofs.equation(element_dofs[2 * j + d]);
                    if (eq_a >= 0 && eq_b >= 0) builder.add_entry(eq_a, eq_b, m(i, j));
                }
    }
}

template <class E>
void assemble_gravity_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<double>& unit_weight,
                           Eigen::VectorXd& rhs,
                           const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;

    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        const double gamma = unit_weight[mesh.element_material[e]];
        if (gamma == 0.0) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);

        // integral_N[i] = ∫ N_i dA = Σ_g w_g detJ_g N_i(gp)
        std::array<double, E::kNodeCount> integral_n{};
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jacobian = dn.transpose() * coords;
            const double w_det = gp.weight * jacobian.determinant();
            for (int i = 0; i < E::kNodeCount; ++i)
                integral_n[i] += w_det * n(i);
        }

        // Downward body force to the uy component only: -γ ∫N_i.
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int eq = dofs.equation(element_dofs[2 * i + 1]);  // uy
            if (eq >= 0) rhs[eq] += -gamma * integral_n[i];
        }
    }
}

// Consistent mass, water-table-aware: if the Gauss point is BELOW the water table the soil
// is SATURATED and the total mass density is ρ_sat = γ_sat/g (grains + pore water);
// ρ_unsat above. The TOTAL mass carries the seismic inertial force −M·r·a_g → using
// γ_unsat everywhere UNDERSTATES the mass (typically 17 vs 20 → 15%) = unsafe-sided. The
// SAME Gauss-point classification as assemble_gravity_phreatic.
template <class E>
void assemble_mass_phreatic_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                 const std::vector<double>& rho_unsat,
                                 const std::vector<double>& rho_sat,
                                 const std::function<double(double)>& water_table_y,
                                 math::SparseMatrixBuilder& builder,
                                 const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        const int mat = mesh.element_material[e];
        const double r_unsat = rho_unsat[mat], r_sat = rho_sat[mat];
        if (r_unsat == 0.0 && r_sat == 0.0) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);

        Eigen::Matrix<double, E::kNodeCount, E::kNodeCount> m =
            Eigen::Matrix<double, E::kNodeCount, E::kNodeCount>::Zero();
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jacobian = dn.transpose() * coords;
            const double w_det = gp.weight * jacobian.determinant();
            double x = 0.0, y = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) { x += n(i) * coords(i, 0); y += n(i) * coords(i, 1); }
            const double rho = (y <= water_table_y(x)) ? r_sat : r_unsat;
            if (rho == 0.0) continue;
            for (int i = 0; i < E::kNodeCount; ++i)
                for (int j = 0; j < E::kNodeCount; ++j)
                    m(i, j) += rho * w_det * n(i) * n(j);
        }
        for (int i = 0; i < E::kNodeCount; ++i)
            for (int j = 0; j < E::kNodeCount; ++j)
                for (int d = 0; d < 2; ++d) {
                    const int eq_a = dofs.equation(element_dofs[2 * i + d]);
                    const int eq_b = dofs.equation(element_dofs[2 * j + d]);
                    if (eq_a >= 0 && eq_b >= 0) builder.add_entry(eq_a, eq_b, m(i, j));
                }
    }
}

template <class E>
void assemble_gravity_phreatic_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                    const std::vector<double>& gamma_unsat,
                                    const std::vector<double>& gamma_sat,
                                    const std::function<double(double)>& water_table_y,
                                    Eigen::VectorXd& rhs,
                                    const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        const int mat = mesh.element_material[e];
        const double g_unsat = gamma_unsat[mat], g_sat = gamma_sat[mat];
        if (g_unsat == 0.0 && g_sat == 0.0) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);

        // γ per Gauss point: saturated below the water table, moist above. -γ ∫N_i (uy).
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jacobian = dn.transpose() * coords;
            const double w_det = gp.weight * jacobian.determinant();
            double x = 0.0, y = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) { x += n(i) * coords(i, 0); y += n(i) * coords(i, 1); }
            const double gamma = (y <= water_table_y(x)) ? g_sat : g_unsat;
            if (gamma == 0.0) continue;
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int eq = dofs.equation(element_dofs[2 * i + 1]);  // uy
                if (eq >= 0) rhs[eq] += -gamma * w_det * n(i);
            }
        }
    }
}

template <class E>
void assemble_internal_force_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const std::vector<GaussState>& states,
                                  Eigen::VectorXd& rhs,
                                  const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto& gp = E::gauss_points()[g];
            const auto sg = E::strain_displacement(coords, gp.xi, gp.eta);
            const double wJ = gp.weight * sg.det_jacobian;
            const Eigen::Vector3d& sigma = states[static_cast<size_t>(e) * E::kGaussCount + g].stress;
            // fe = wJ · Bᵀ σ  → node i: [B(0,2i)·sxx + B(2,2i)·sxy,  B(1,2i+1)·syy + B(2,2i+1)·sxy].
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int ex = dofs.equation(element_dofs[2 * i]);
                const int ey = dofs.equation(element_dofs[2 * i + 1]);
                if (ex >= 0) rhs[ex] += wJ * (sg.B(0, 2 * i) * sigma(0) + sg.B(2, 2 * i) * sigma(2));
                if (ey >= 0) rhs[ey] += wJ * (sg.B(1, 2 * i + 1) * sigma(1) + sg.B(2, 2 * i + 1) * sigma(2));
            }
        }
    }
}

template <class E>
void assemble_pore_load_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::function<double(double, double)>& pore,
                             Eigen::VectorXd& rhs,
                             const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        for (const auto& gp : E::gauss_points()) {
            const auto g = E::strain_displacement(coords, gp.xi, gp.eta);
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            double x = 0.0, y = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) {
                x += n(i) * coords(i, 0);
                y += n(i) * coords(i, 1);
            }
            const double u = pore(x, y);
            if (u == 0.0) continue;
            const double w = gp.weight * g.det_jacobian * u;
            // Bᵀ m with m = [1,1,0]: node i -> [dN/dx, dN/dy] = [B(0,2i), B(1,2i+1)].
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int ex = dofs.equation(element_dofs[2 * i]);
                const int ey = dofs.equation(element_dofs[2 * i + 1]);
                if (ex >= 0) rhs[ex] += w * g.B(0, 2 * i);
                if (ey >= 0) rhs[ey] += w * g.B(1, 2 * i + 1);
            }
        }
    }
}

template <class E>
void assemble_axisym_stiffness_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                    const std::vector<LinearElastic>& materials,
                                    math::SparseMatrixBuilder& builder) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        const LinearElastic& material = materials[mesh.element_material[e]];
        const typename E::ElementMatrix ke =
            axisym::element_stiffness<E>(coords, material.axisymmetric_matrix());
        for (int a = 0; a < E::kDofCount; ++a) {
            const int eq_a = dofs.equation(element_dofs[a]);
            if (eq_a < 0) continue;
            for (int b = 0; b < E::kDofCount; ++b) {
                const int eq_b = dofs.equation(element_dofs[b]);
                if (eq_b < 0) continue;
                builder.add_entry(eq_a, eq_b, ke(a, b));
            }
        }
    }
}

// Axisymmetric body force: r-weighted downward gravity, f_z += -gamma * integral(N_i r dA).
template <class E>
void assemble_axisym_gravity_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const std::vector<double>& unit_weight,
                                  Eigen::VectorXd& rhs,
                                  const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        const double gamma = unit_weight[mesh.element_material[e]];
        if (gamma == 0.0) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        std::array<double, E::kNodeCount> integral_n{};
        for (const auto& gp : E::gauss_points()) {
            const typename E::ShapeValues n = E::shape_functions(gp.xi, gp.eta);
            const auto dn = E::shape_derivatives_natural(gp.xi, gp.eta);
            const Eigen::Matrix2d jacobian = dn.transpose() * coords;
            double r = 0.0;
            for (int i = 0; i < E::kNodeCount; ++i) r += n(i) * coords(i, 0);  // x = radius
            const double w_det_r = gp.weight * jacobian.determinant() * r;
            for (int i = 0; i < E::kNodeCount; ++i) integral_n[i] += w_det_r * n(i);
        }
        for (int i = 0; i < E::kNodeCount; ++i) {
            const int eq = dofs.equation(element_dofs[2 * i + 1]);  // u_z
            if (eq >= 0) rhs[eq] += -gamma * integral_n[i];
        }
    }
}

// Axisymmetric consistent nodal internal force F = integral B^T sigma * (w detJ r); sigma carries the
// hoop component sigma_theta = GaussState::stress_zz (B's 4th row couples it to the radial DOF).
template <class E>
void assemble_axisym_internal_force_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                         const std::vector<GaussState>& states,
                                         Eigen::VectorXd& rhs,
                                         const std::vector<char>& active_element) {
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> element_dofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        gather_element<E>(mesh, dofs, e, coords, element_dofs);
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto& gp = E::gauss_points()[g];
            const auto sg = axisym::strain_displacement<E>(coords, gp.xi, gp.eta);
            const GaussState& st = states[static_cast<size_t>(e) * E::kGaussCount + g];
            Eigen::Vector4d sig(st.stress(0), st.stress(1), st.stress(2), st.stress_zz);
            const double scale = gp.weight * sg.det_jacobian * sg.radius;
            const Eigen::Matrix<double, E::kDofCount, 1> fe = scale * (sg.B.transpose() * sig);
            for (int i = 0; i < E::kNodeCount; ++i) {
                const int ex = dofs.equation(element_dofs[2 * i]);
                const int ey = dofs.equation(element_dofs[2 * i + 1]);
                if (ex >= 0) rhs[ex] += fe(2 * i);
                if (ey >= 0) rhs[ey] += fe(2 * i + 1);
            }
        }
    }
}

} // namespace

void assemble_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                        const std::vector<LinearElastic>& materials,
                        math::SparseMatrixBuilder& builder,
                        const std::vector<char>& active_element,
                        const std::vector<MaterialProfile>& profile,
                        const std::vector<LinearElastic>& gauss_elastic) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_stiffness_impl<Tri15Element>(mesh, dofs, materials, builder, active_element, profile,
                                              gauss_elastic);
    else
        assemble_stiffness_impl<Tri6Element>(mesh, dofs, materials, builder, active_element, profile,
                                             gauss_elastic);
}

void assemble_mass(const mesh::Mesh& mesh, const DofMap& dofs,
                   const std::vector<double>& density,
                   math::SparseMatrixBuilder& builder,
                   const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_mass_impl<Tri15Element>(mesh, dofs, density, builder, active_element);
    else
        assemble_mass_impl<Tri6Element>(mesh, dofs, density, builder, active_element);
}

void assemble_mass_phreatic(const mesh::Mesh& mesh, const DofMap& dofs,
                            const std::vector<double>& rho_unsat,
                            const std::vector<double>& rho_sat,
                            const std::function<double(double)>& water_table_y,
                            math::SparseMatrixBuilder& builder,
                            const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_mass_phreatic_impl<Tri15Element>(mesh, dofs, rho_unsat, rho_sat, water_table_y,
                                                  builder, active_element);
    else
        assemble_mass_phreatic_impl<Tri6Element>(mesh, dofs, rho_unsat, rho_sat, water_table_y,
                                                 builder, active_element);
}

void assemble_gravity(const mesh::Mesh& mesh, const DofMap& dofs,
                      const std::vector<double>& unit_weight,
                      Eigen::VectorXd& rhs,
                      const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_gravity_impl<Tri15Element>(mesh, dofs, unit_weight, rhs,
                                            active_element);
    else
        assemble_gravity_impl<Tri6Element>(mesh, dofs, unit_weight, rhs,
                                           active_element);
}

void assemble_gravity_phreatic(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<double>& gamma_unsat,
                               const std::vector<double>& gamma_sat,
                               const std::function<double(double)>& water_table_y,
                               Eigen::VectorXd& rhs,
                               const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_gravity_phreatic_impl<Tri15Element>(mesh, dofs, gamma_unsat, gamma_sat,
                                                     water_table_y, rhs, active_element);
    else
        assemble_gravity_phreatic_impl<Tri6Element>(mesh, dofs, gamma_unsat, gamma_sat,
                                                    water_table_y, rhs, active_element);
}

void assemble_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<GaussState>& states,
                             Eigen::VectorXd& rhs,
                             const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_internal_force_impl<Tri15Element>(mesh, dofs, states, rhs, active_element);
    else
        assemble_internal_force_impl<Tri6Element>(mesh, dofs, states, rhs, active_element);
}

void assemble_axisym_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<LinearElastic>& materials,
                               math::SparseMatrixBuilder& builder) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_axisym_stiffness_impl<Tri15Element>(mesh, dofs, materials, builder);
    else
        assemble_axisym_stiffness_impl<Tri6Element>(mesh, dofs, materials, builder);
}

void assemble_pore_pressure_load(const mesh::Mesh& mesh, const DofMap& dofs,
                                 const std::function<double(double, double)>& pore,
                                 Eigen::VectorXd& rhs,
                                 const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_pore_load_impl<Tri15Element>(mesh, dofs, pore, rhs, active_element);
    else
        assemble_pore_load_impl<Tri6Element>(mesh, dofs, pore, rhs, active_element);
}

namespace {

// 1D Lagrange edge shape functions N_k and derivatives dN_k/ds, for n equally spaced
// nodes on [-1,1] (n=3 quadratic tri6, n=5 quartic tri15 edge).
void edge_shape(int n, double s, double* N, double* dN) {
    for (int k = 0; k < n; ++k) {
        const double sk = -1.0 + 2.0 * k / (n - 1);
        double val = 1.0, der = 0.0;
        for (int m = 0; m < n; ++m) {
            if (m == k) continue;
            const double sm = -1.0 + 2.0 * m / (n - 1);
            const double inv = 1.0 / (sk - sm);
            der = der * (s - sm) * inv + val * inv;  // product rule
            val *= (s - sm) * inv;
        }
        N[k] = val;
        dN[k] = der;
    }
}

} // namespace

void assemble_surface_traction(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<int>& ordered_boundary_nodes,
                               double traction_x, double traction_y,
                               Eigen::VectorXd& rhs) {
    // Nodes per edge: tri6 → 3 (quadratic), tri15 → 5 (quartic).
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    assert(chain >= npe && (chain - 1) % (npe - 1) == 0 && "kenar zinciri uyumsuz");
    const int edge_count = (chain - 1) / (npe - 1);

    // 4-point 1D Gauss-Legendre on [-1,1] (degree 7 — more than exact for the quartic
    // edge, ready for curved/isoparametric edges too).
    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        // ∫ N_i ds for the edge's npe nodes.
        std::array<double, 5> integral_n{};
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
            const int eq_x = dofs.equation(dofs.global_dof(ni, 0));
            const int eq_y = dofs.equation(dofs.global_dof(ni, 1));
            if (eq_x >= 0) rhs[eq_x] += traction_x * integral_n[i];
            if (eq_y >= 0) rhs[eq_y] += traction_y * integral_n[i];
        }
    }
}

void assemble_boundary_dashpot(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<int>& ordered_boundary_nodes,
                               double c_n, double c_t,
                               math::SparseMatrixBuilder& builder) {
    // Lysmer-Kuhlemeyer viscous dashpot: C_b = ∫ Nᵀ D_c N ds, D_c = c_n (n⊗n) + c_t (t⊗t).
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    if (chain < npe || (chain - 1) % (npe - 1) != 0) return;
    const int edge_count = (chain - 1) / (npe - 1);

    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        // Cb[2i+a][2j+b] = ∫ N_i N_j D_c[a][b] ds  (normal/tangent vary along a curved edge -> per gp).
        double Cb[10][10] = {};
        for (int q = 0; q < 4; ++q) {
            edge_shape(npe, point[q], N.data(), dN.data());
            double dxdt = 0.0, dydt = 0.0;
            for (int i = 0; i < npe; ++i) {
                const int ni = ordered_boundary_nodes[base + i];
                dxdt += dN[i] * mesh.x[ni];
                dydt += dN[i] * mesh.y[ni];
            }
            const double ds = std::sqrt(dxdt * dxdt + dydt * dydt);
            if (ds < 1e-30) continue;
            const double tx = dxdt / ds, ty = dydt / ds;   // unit tangent
            const double nx = ty, ny = -tx;                // unit normal (sign irrelevant for D_c)
            const double D00 = c_n * nx * nx + c_t * tx * tx;   // D_c (2x2 symmetric)
            const double D01 = c_n * nx * ny + c_t * tx * ty;
            const double D11 = c_n * ny * ny + c_t * ty * ty;
            const double w_ds = weight[q] * ds;
            for (int i = 0; i < npe; ++i)
                for (int j = 0; j < npe; ++j) {
                    const double NiNj = N[i] * N[j] * w_ds;
                    Cb[2 * i][2 * j] += NiNj * D00;
                    Cb[2 * i][2 * j + 1] += NiNj * D01;
                    Cb[2 * i + 1][2 * j] += NiNj * D01;
                    Cb[2 * i + 1][2 * j + 1] += NiNj * D11;
                }
        }
        for (int i = 0; i < npe; ++i) {
            const int ni = ordered_boundary_nodes[base + i];
            const int eqi[2] = {dofs.equation(dofs.global_dof(ni, 0)),
                                dofs.equation(dofs.global_dof(ni, 1))};
            for (int j = 0; j < npe; ++j) {
                const int nj = ordered_boundary_nodes[base + j];
                const int eqj[2] = {dofs.equation(dofs.global_dof(nj, 0)),
                                    dofs.equation(dofs.global_dof(nj, 1))};
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b)
                        if (eqi[a] >= 0 && eqj[b] >= 0)
                            builder.add_entry(eqi[a], eqj[b], Cb[2 * i + a][2 * j + b]);
            }
        }
    }
}

void assemble_surface_traction_varying(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const std::vector<int>& ordered_boundary_nodes,
                                       const std::vector<double>& tx,
                                       const std::vector<double>& ty,
                                       Eigen::VectorXd& rhs) {
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    assert(chain >= npe && (chain - 1) % (npe - 1) == 0 && "kenar zinciri uyumsuz");
    assert(static_cast<int>(tx.size()) == chain && static_cast<int>(ty.size()) == chain);
    const int edge_count = (chain - 1) / (npe - 1);

    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        // Consistent load f_i += ∫ N_i (Σ_j N_j t_j) ds: the traction is interpolated by the SAME
        // edge shape functions as the geometry, so a linearly-varying q (q1->q2) is reproduced exactly.
        std::array<double, 5> fx{}, fy{};
        for (int q = 0; q < 4; ++q) {
            edge_shape(npe, point[q], N.data(), dN.data());
            double dxdt = 0.0, dydt = 0.0, tqx = 0.0, tqy = 0.0;
            for (int i = 0; i < npe; ++i) {
                const int ni = ordered_boundary_nodes[base + i];
                dxdt += dN[i] * mesh.x[ni];
                dydt += dN[i] * mesh.y[ni];
                tqx += N[i] * tx[base + i];   // interpolated traction at the Gauss point
                tqy += N[i] * ty[base + i];
            }
            const double ds = std::sqrt(dxdt * dxdt + dydt * dydt);
            for (int i = 0; i < npe; ++i) {
                fx[i] += weight[q] * ds * N[i] * tqx;
                fy[i] += weight[q] * ds * N[i] * tqy;
            }
        }
        for (int i = 0; i < npe; ++i) {
            const int ni = ordered_boundary_nodes[base + i];
            const int eq_x = dofs.equation(dofs.global_dof(ni, 0));
            const int eq_y = dofs.equation(dofs.global_dof(ni, 1));
            if (eq_x >= 0) rhs[eq_x] += fx[i];
            if (eq_y >= 0) rhs[eq_y] += fy[i];
        }
    }
}

void assemble_axisym_traction(const mesh::Mesh& mesh, const DofMap& dofs,
                              const std::vector<int>& ordered_boundary_nodes,
                              double traction_r, double traction_z,
                              Eigen::VectorXd& rhs) {
    // r weight on the surface-traction integral: ∫ N_i (t_r, t_z) r ds (per radian).
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    assert(chain >= npe && (chain - 1) % (npe - 1) == 0 && "kenar zinciri uyumsuz");
    const int edge_count = (chain - 1) / (npe - 1);

    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        std::array<double, 5> integral_n{};
        for (int q = 0; q < 4; ++q) {
            edge_shape(npe, point[q], N.data(), dN.data());
            double dxdt = 0.0, dydt = 0.0, r = 0.0;
            for (int i = 0; i < npe; ++i) {
                const int ni = ordered_boundary_nodes[base + i];
                dxdt += dN[i] * mesh.x[ni];
                dydt += dN[i] * mesh.y[ni];
                r += N[i] * mesh.x[ni];  // radius (x = r)
            }
            const double ds = std::sqrt(dxdt * dxdt + dydt * dydt);
            for (int i = 0; i < npe; ++i) integral_n[i] += weight[q] * ds * r * N[i];
        }
        for (int i = 0; i < npe; ++i) {
            const int ni = ordered_boundary_nodes[base + i];
            const int eq_r = dofs.equation(dofs.global_dof(ni, 0));
            const int eq_z = dofs.equation(dofs.global_dof(ni, 1));
            if (eq_r >= 0) rhs[eq_r] += traction_r * integral_n[i];
            if (eq_z >= 0) rhs[eq_z] += traction_z * integral_n[i];
        }
    }
}

void assemble_axisym_gravity(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<double>& unit_weight, Eigen::VectorXd& rhs,
                             const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_axisym_gravity_impl<Tri15Element>(mesh, dofs, unit_weight, rhs, active_element);
    else
        assemble_axisym_gravity_impl<Tri6Element>(mesh, dofs, unit_weight, rhs, active_element);
}

void assemble_axisym_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                                    const std::vector<GaussState>& states, Eigen::VectorXd& rhs,
                                    const std::vector<char>& active_element) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        assemble_axisym_internal_force_impl<Tri15Element>(mesh, dofs, states, rhs, active_element);
    else
        assemble_axisym_internal_force_impl<Tri6Element>(mesh, dofs, states, rhs, active_element);
}

// Axisymmetric varying edge traction: f += integral N_i (Sum_j N_j t_j) r ds (r-weighted version of
// assemble_surface_traction_varying; the traction is interpolated by the edge shape functions).
void assemble_axisym_traction_varying(const mesh::Mesh& mesh, const DofMap& dofs,
                                      const std::vector<int>& ordered_boundary_nodes,
                                      const std::vector<double>& tx, const std::vector<double>& ty,
                                      Eigen::VectorXd& rhs) {
    const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
    const int chain = static_cast<int>(ordered_boundary_nodes.size());
    assert(chain >= npe && (chain - 1) % (npe - 1) == 0 && "kenar zinciri uyumsuz");
    assert(static_cast<int>(tx.size()) == chain && static_cast<int>(ty.size()) == chain);
    const int edge_count = (chain - 1) / (npe - 1);

    constexpr double g1 = 0.3399810435848563, g2 = 0.8611363115940526;
    constexpr double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
    const std::array<double, 4> point = {-g2, -g1, g1, g2};
    const std::array<double, 4> weight = {w2, w1, w1, w2};

    std::array<double, 5> N{}, dN{};
    for (int edge = 0; edge < edge_count; ++edge) {
        const int base = edge * (npe - 1);
        std::array<double, 5> fx{}, fy{};
        for (int q = 0; q < 4; ++q) {
            edge_shape(npe, point[q], N.data(), dN.data());
            double dxdt = 0.0, dydt = 0.0, r = 0.0, tqx = 0.0, tqy = 0.0;
            for (int i = 0; i < npe; ++i) {
                const int ni = ordered_boundary_nodes[base + i];
                dxdt += dN[i] * mesh.x[ni];
                dydt += dN[i] * mesh.y[ni];
                r += N[i] * mesh.x[ni];       // radius (x = r)
                tqx += N[i] * tx[base + i];
                tqy += N[i] * ty[base + i];
            }
            const double ds = std::sqrt(dxdt * dxdt + dydt * dydt);
            for (int i = 0; i < npe; ++i) {
                fx[i] += weight[q] * ds * r * N[i] * tqx;
                fy[i] += weight[q] * ds * r * N[i] * tqy;
            }
        }
        for (int i = 0; i < npe; ++i) {
            const int ni = ordered_boundary_nodes[base + i];
            const int eq_x = dofs.equation(dofs.global_dof(ni, 0));
            const int eq_y = dofs.equation(dofs.global_dof(ni, 1));
            if (eq_x >= 0) rhs[eq_x] += fx[i];
            if (eq_y >= 0) rhs[eq_y] += fy[i];
        }
    }
}

Eigen::VectorXd expand_to_full(const DofMap& dofs,
                               const Eigen::VectorXd& free_solution) {
    Eigen::VectorXd full = Eigen::VectorXd::Zero(dofs.total_dofs());
    for (int d = 0; d < dofs.total_dofs(); ++d) {
        const int eq = dofs.equation(d);
        if (eq >= 0) full(d) = free_solution(eq);
    }
    return full;
}

} // namespace katai::core
