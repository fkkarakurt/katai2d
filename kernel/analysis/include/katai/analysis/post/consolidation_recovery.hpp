#pragma once
// Consolidation stress recovery (Stage B8: extracted from the application
// driver; the former ad-hoc consol_detail namespace resolves to the module's
// detail convention here). After the time-dependent solve, the committed Gauss
// stresses are recovered: the linear-elastic skeleton means the effective-
// stress increment is just D * (B v_final), added to the phase's initial
// committed stress (`out` starts as a copy of `init`). Plane strain: the
// out-of-plane increment is sigma_zz += nu (sigma_xx + sigma_yy), and the
// volumetric strain accumulates eps_xx + eps_yy.

#include <vector>

#include <Eigen/Core>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/tri15.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {
namespace detail {

template <class E>
inline void recover_stress_impl(const katai::mesh::Mesh& mesh, const DofMap& dofs,
                                const std::vector<MaterialModel>& models,
                                const Eigen::VectorXd& vfull, const std::vector<char>& act,
                                std::vector<GaussState>& out) {
    constexpr int N = E::kNodeCount;
    const auto gp = E::gauss_points();
    const int ng = E::kGaussCount;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!act.empty() && !act[e]) continue;
        typename E::NodeCoords X;
        Eigen::Matrix<double, 2 * N, 1> ve;
        for (int k = 0; k < N; ++k) {
            const int nd = mesh.node_of(e, k);
            X(k, 0) = mesh.x[nd]; X(k, 1) = mesh.y[nd];
            ve(2 * k) = vfull[dofs.global_dof(nd, 0)];
            ve(2 * k + 1) = vfull[dofs.global_dof(nd, 1)];
        }
        const auto& mm = models[mesh.element_material[e]];
        const Eigen::Matrix3d D = mm.elastic_plane_strain();
        const double nu = mm.poisson_ratio;
        for (int g = 0; g < ng; ++g) {
            const auto sd = E::strain_displacement(X, gp[g].xi, gp[g].eta);
            const Eigen::Vector3d deps = sd.B * ve;
            const Eigen::Vector3d dsig = D * deps;
            auto& gs = out[(size_t)e * ng + g];
            gs.stress += dsig;
            gs.stress_zz += nu * (dsig(0) + dsig(1));
            gs.eps_vol += deps(0) + deps(1);
        }
    }
}

} // namespace detail

inline void recover_consolidation_stress(const katai::mesh::Mesh& mesh, const DofMap& dofs,
                                         const std::vector<MaterialModel>& models,
                                         const Eigen::VectorXd& vfull, const std::vector<char>& act,
                                         std::vector<GaussState>& out) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        detail::recover_stress_impl<Tri15Element>(mesh, dofs, models, vfull, act, out);
    else
        detail::recover_stress_impl<Tri6Element>(mesh, dofs, models, vfull, act, out);
}

} // namespace katai::core
