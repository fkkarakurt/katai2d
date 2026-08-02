#include <katai/analysis/post/stress_recovery.hpp>

#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/fem/elements/tri15.hpp>

namespace katai::core {
namespace {

// Extrapolation matrix from the 3 Gauss points to the 6 nodes (6x3), built once.
// Linear field: s(ξ,η)=a+bξ+cη. The coefficients are found by inverting the [1,ξ,η]
// matrix at the Gauss points, then evaluated at the nodal natural coordinates.
Eigen::Matrix<double, 6, 3> build_extrapolation_matrix() {
    const auto gauss = tri6::gauss_points();
    Eigen::Matrix3d at_gauss;
    for (int g = 0; g < 3; ++g)
        at_gauss.row(g) << 1.0, gauss[g].xi, gauss[g].eta;

    // Nodal natural coordinates (corners, then mid-edges).
    Eigen::Matrix<double, 6, 2> nat;
    nat << 0.0, 0.0,  1.0, 0.0,  0.0, 1.0,
           0.5, 0.0,  0.5, 0.5,  0.0, 0.5;
    Eigen::Matrix<double, 6, 3> at_nodes;
    for (int n = 0; n < 6; ++n)
        at_nodes.row(n) << 1.0, nat(n, 0), nat(n, 1);

    return at_nodes * at_gauss.inverse();
}

const Eigen::Matrix<double, 6, 3>& extrapolation_matrix() {
    static const Eigen::Matrix<double, 6, 3> matrix = build_extrapolation_matrix();
    return matrix;
}

} // namespace

std::array<Eigen::Vector3d, tri6::kGaussCount> gauss_point_stresses(
    const tri6::NodeCoords& coords, const Eigen::Matrix3d& constitutive,
    const ElementDisplacement& element_displacement) {
    std::array<Eigen::Vector3d, tri6::kGaussCount> stresses;
    int g = 0;
    for (const auto& gp : tri6::gauss_points()) {
        const auto grad = tri6::strain_displacement(coords, gp.xi, gp.eta);
        const Eigen::Vector3d strain = grad.B * element_displacement;
        stresses[g++] = constitutive * strain;
    }
    return stresses;
}

std::array<Eigen::Vector3d, tri6::kNodeCount> extrapolate_gauss_to_nodes(
    const std::array<Eigen::Vector3d, tri6::kGaussCount>& gauss_values) {
    Eigen::Matrix3d at_gauss;  // row = Gauss point, column = stress component
    for (int g = 0; g < 3; ++g) at_gauss.row(g) = gauss_values[g].transpose();

    const Eigen::Matrix<double, 6, 3> at_nodes =
        extrapolation_matrix() * at_gauss;

    std::array<Eigen::Vector3d, tri6::kNodeCount> nodal;
    for (int n = 0; n < tri6::kNodeCount; ++n)
        nodal[n] = at_nodes.row(n).transpose();
    return nodal;
}

namespace {

// Row matrix of degree-d polynomial monomials (ξ,η) at a point set (ξ^a η^b per column, a+b≤d).
Eigen::MatrixXd poly_basis(const std::vector<std::array<double, 2>>& pts, int degree) {
    std::vector<std::array<int, 2>> terms;
    for (int d = 0; d <= degree; ++d)
        for (int a = d; a >= 0; --a) terms.push_back({a, d - a});
    Eigen::MatrixXd P(static_cast<int>(pts.size()), static_cast<int>(terms.size()));
    for (size_t i = 0; i < pts.size(); ++i)
        for (size_t j = 0; j < terms.size(); ++j)
            P(static_cast<int>(i), static_cast<int>(j)) =
                std::pow(pts[i][0], terms[j][0]) * std::pow(pts[i][1], terms[j][1]);
    return P;
}

// Gauss→node extrapolation matrix (nNode×nGauss): nodal = M·gauss. Fits a degree-d
// polynomial to the Gauss values by least squares (nGauss≥nTerms), evaluates it at the
// nodal natural coordinates.
// M = P_node · (P_gaussᵀ P_gauss)⁻¹ P_gaussᵀ.
Eigen::MatrixXd build_extrapolation(const std::vector<std::array<double, 2>>& gauss,
                                    const std::vector<std::array<double, 2>>& nodes, int degree) {
    const Eigen::MatrixXd Pg = poly_basis(gauss, degree);
    const Eigen::MatrixXd Pn = poly_basis(nodes, degree);
    const Eigen::MatrixXd pinv =
        (Pg.transpose() * Pg).ldlt().solve(Pg.transpose());   // nTerms×nGauss
    return Pn * pinv;                                          // nNode×nGauss
}

const Eigen::MatrixXd& tri6_gauss_to_node() {
    static const Eigen::MatrixXd M = [] {
        std::vector<std::array<double, 2>> g;
        for (const auto& gp : tri6::gauss_points()) g.push_back({gp.xi, gp.eta});
        const std::vector<std::array<double, 2>> n = {
            {0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5}};
        return build_extrapolation(g, n, 1);
    }();
    return M;
}

const Eigen::MatrixXd& tri15_gauss_to_node() {
    static const Eigen::MatrixXd M = [] {
        std::vector<std::array<double, 2>> g;
        for (const auto& gp : tri15::gauss_points()) g.push_back({gp.xi, gp.eta});
        // Nodal natural coordinates (tri15.cpp kBary order: (b/4, c/4)).
        const std::vector<std::array<double, 2>> n = {
            {0, 0}, {1, 0}, {0, 1},
            {0.25, 0}, {0.5, 0}, {0.75, 0},
            {0.75, 0.25}, {0.5, 0.5}, {0.25, 0.75},
            {0, 0.75}, {0, 0.5}, {0, 0.25},
            {0.25, 0.25}, {0.5, 0.25}, {0.25, 0.5}};
        return build_extrapolation(g, n, 3);   // stress ~ cubic (quartic displacement)
    }();
    return M;
}

} // namespace

NodalStressField recover_nodal_stresses_from_gauss(
    const mesh::Mesh& mesh, const std::vector<GaussState>& gauss_states,
    const std::vector<char>& active_element) {
    NodalStressField field;
    field.stress.assign(mesh.node_count, Eigen::Vector3d::Zero());
    std::vector<int> contributions(mesh.node_count, 0);

    const bool t15 = mesh.nodes_per_element == tri15::kNodeCount;
    const int nN = t15 ? tri15::kNodeCount : tri6::kNodeCount;
    const int nG = t15 ? tri15::kGaussCount : tri6::kGaussCount;
    const Eigen::MatrixXd& M = t15 ? tri15_gauss_to_node() : tri6_gauss_to_node();

    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;   // excavated: no stress shown
        for (int k = 0; k < nN; ++k) {
            Eigen::Vector3d s = Eigen::Vector3d::Zero();
            for (int g = 0; g < nG; ++g)
                s += M(k, g) * gauss_states[static_cast<size_t>(e) * nG + g].stress;
            const int n = mesh.node_of(e, k);
            field.stress[n] += s;
            ++contributions[n];
        }
    }
    for (int n = 0; n < mesh.node_count; ++n)
        if (contributions[n] > 0) field.stress[n] /= static_cast<double>(contributions[n]);
    return field;
}

namespace {

// F += B^T sigma over the active elements -- the discrete internal force, exact for the
// committed states (no recovery smoothing involved). Plane strain and axisymmetric
// accumulators; the axisymmetric one is r-weighted per radian and carries the hoop term.
template <class E>
void accumulate_f_int_ps(const mesh::Mesh& mesh, const std::vector<GaussState>& gs,
                         const std::vector<char>& act, Eigen::VectorXd& F) {
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!act.empty() && !act[e]) continue;
        typename E::NodeCoords X;
        for (int k = 0; k < E::kNodeCount; ++k) {
            X(k, 0) = mesh.x[mesh.node_of(e, k)];
            X(k, 1) = mesh.y[mesh.node_of(e, k)];
        }
        Eigen::Matrix<double, E::kDofCount, 1> fe = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sg = E::strain_displacement(X, gp[g].xi, gp[g].eta);
            fe.noalias() += (gp[g].weight * sg.det_jacobian) * sg.B.transpose() *
                            gs[static_cast<size_t>(e) * E::kGaussCount + g].stress;
        }
        for (int k = 0; k < E::kNodeCount; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
}

template <class E>
void accumulate_f_int_axi(const mesh::Mesh& mesh, const std::vector<GaussState>& gs,
                          const std::vector<char>& act, Eigen::VectorXd& F) {
    const auto gp = E::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!act.empty() && !act[e]) continue;
        typename E::NodeCoords X;
        for (int k = 0; k < E::kNodeCount; ++k) {
            X(k, 0) = mesh.x[mesh.node_of(e, k)];
            X(k, 1) = mesh.y[mesh.node_of(e, k)];
        }
        Eigen::Matrix<double, E::kDofCount, 1> fe = Eigen::Matrix<double, E::kDofCount, 1>::Zero();
        for (int g = 0; g < E::kGaussCount; ++g) {
            const auto sg = axisym::strain_displacement<E>(X, gp[g].xi, gp[g].eta);
            const auto& s = gs[static_cast<size_t>(e) * E::kGaussCount + g];
            Eigen::Vector4d sig;
            sig << s.stress(0), s.stress(1), s.stress(2), s.stress_zz;
            fe.noalias() += (gp[g].weight * sg.det_jacobian * sg.radius) * sg.B.transpose() * sig;
        }
        for (int k = 0; k < E::kNodeCount; ++k) {
            F(2 * mesh.node_of(e, k) + 0) += fe(2 * k + 0);
            F(2 * mesh.node_of(e, k) + 1) += fe(2 * k + 1);
        }
    }
}

}  // namespace

Eigen::VectorXd nodal_internal_force_from_gauss(const mesh::Mesh& mesh,
                                                const std::vector<GaussState>& gauss_states,
                                                bool axisymmetric,
                                                const std::vector<char>& active_element) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * mesh.node_count);
    const bool t15 = mesh.nodes_per_element == tri15::kNodeCount;
    if (axisymmetric) {
        if (t15) accumulate_f_int_axi<Tri15Element>(mesh, gauss_states, active_element, F);
        else     accumulate_f_int_axi<Tri6Element>(mesh, gauss_states, active_element, F);
    } else {
        if (t15) accumulate_f_int_ps<Tri15Element>(mesh, gauss_states, active_element, F);
        else     accumulate_f_int_ps<Tri6Element>(mesh, gauss_states, active_element, F);
    }
    return F;
}

NodalStressField recover_nodal_stresses(
    const mesh::Mesh& mesh, const std::vector<LinearElastic>& materials,
    const Eigen::VectorXd& full_displacement) {
    NodalStressField field;
    field.stress.assign(mesh.node_count, Eigen::Vector3d::Zero());
    std::vector<int> contributions(mesh.node_count, 0);

    tri6::NodeCoords coords;
    ElementDisplacement element_displacement;

    for (int e = 0; e < mesh.element_count; ++e) {
        for (int k = 0; k < tri6::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            coords(k, 0) = mesh.x[n];
            coords(k, 1) = mesh.y[n];
            element_displacement(2 * k) = full_displacement(2 * n);
            element_displacement(2 * k + 1) = full_displacement(2 * n + 1);
        }
        const Eigen::Matrix3d constitutive =
            materials[mesh.element_material[e]].plane_strain_matrix();

        const auto gauss = gauss_point_stresses(coords, constitutive,
                                                element_displacement);
        const auto nodal = extrapolate_gauss_to_nodes(gauss);

        for (int k = 0; k < tri6::kNodeCount; ++k) {
            const int n = mesh.node_of(e, k);
            field.stress[n] += nodal[k];
            ++contributions[n];
        }
    }

    for (int n = 0; n < mesh.node_count; ++n)
        if (contributions[n] > 0)
            field.stress[n] /= static_cast<double>(contributions[n]);

    return field;
}

} // namespace katai::core
