#pragma once
// Element traits — gathers the tri6/tri15 common interface under one type for
// compile-time (template) dispatch. Decision (DOD, see ARCHITECTURE): NO vtable on
// the hot path; assembly/solver functions are templated over the Element trait,
// selected by a single runtime branch on `mesh.nodes_per_element` at the top level,
// and the inner loop stays fully monomorphic/inlined. Since both element namespaces
// already expose the same names (kNodeCount, gauss_points, strain_displacement, ...)
// the trait merely forwards.

#include <katai/fem/elements/tri15.hpp>
#include <katai/fem/elements/tri6.hpp>

namespace katai::core {

struct Tri6Element {
    static constexpr int kNodeCount = tri6::kNodeCount;  // 6
    static constexpr int kDofCount = tri6::kDofCount;    // 12
    static constexpr int kGaussCount = tri6::kGaussCount;  // 3
    using NodeCoords = tri6::NodeCoords;
    using ElementMatrix = tri6::ElementMatrix;
    using ShapeValues = tri6::ShapeValues;
    using ShapeGradients = tri6::ShapeGradients;
    static auto gauss_points() { return tri6::gauss_points(); }
    static ShapeValues shape_functions(double xi, double eta) {
        return tri6::shape_functions(xi, eta);
    }
    static tri6::ShapeDerivs shape_derivatives_natural(double xi, double eta) {
        return tri6::shape_derivatives_natural(xi, eta);
    }
    static ShapeGradients strain_displacement(const NodeCoords& c, double xi,
                                              double eta) {
        return tri6::strain_displacement(c, xi, eta);
    }
    static ElementMatrix element_stiffness(const NodeCoords& c,
                                           const Eigen::Matrix3d& d) {
        return tri6::element_stiffness(c, d);
    }
};

struct Tri15Element {
    static constexpr int kNodeCount = tri15::kNodeCount;  // 15
    static constexpr int kDofCount = tri15::kDofCount;    // 30
    static constexpr int kGaussCount = tri15::kGaussCount;  // 12
    using NodeCoords = tri15::NodeCoords;
    using ElementMatrix = tri15::ElementMatrix;
    using ShapeValues = tri15::ShapeValues;
    using ShapeGradients = tri15::ShapeGradients;
    static auto gauss_points() { return tri15::gauss_points(); }
    static ShapeValues shape_functions(double xi, double eta) {
        return tri15::shape_functions(xi, eta);
    }
    static tri15::ShapeDerivs shape_derivatives_natural(double xi, double eta) {
        return tri15::shape_derivatives_natural(xi, eta);
    }
    static ShapeGradients strain_displacement(const NodeCoords& c, double xi,
                                              double eta) {
        return tri15::strain_displacement(c, xi, eta);
    }
    static ElementMatrix element_stiffness(const NodeCoords& c,
                                           const Eigen::Matrix3d& d) {
        return tri15::element_stiffness(c, d);
    }
};

} // namespace katai::core
