#pragma once
// DOF (degree-of-freedom) manager — node → global equation index.
//
// Phase 0: 2 DOFs per node (ux, uy). Some DOFs are fixed with zero-Dirichlet
// (base clamped, side rollers). Fixed DOFs are ELIMINATED from the equation
// system; free DOFs get compacted equation indices → the reduced system.
// This approach beats the penalty method: the matrix is both smaller (faster)
// and, its conditioning intact, the result stays exactly correct.

#include <vector>

namespace katai::core {

class DofMap {
public:
    explicit DofMap(int node_count, int dofs_per_node = 2);

    int dofs_per_node() const { return dofs_per_node_; }
    int node_count() const { return node_count_; }
    int total_dofs() const { return total_dofs_; }

    int global_dof(int node, int component) const {
        return node * dofs_per_node_ + component;
    }

    // Add an extra (structural) DOF — like the plate rotational DOF (φ). Appended AFTER
    // the base translational DOFs (node*dofs_per_node); enters the global equation system.
    // The returned global DOF index is stored by the caller (e.g. the plate node →
    // rotation-DOF mapping). Translational DOFs are SHARED with the soil (same
    // global_dof), the rotation is plate-specific. BEFORE finalize.
    int add_extra_dof();

    // Mark the DOF as zero-Dirichlet (called before finalize).
    void fix(int global_dof);
    void fix_node_component(int node, int component) {
        fix(global_dof(node, component));
    }

    // Serbest DOF'lara 0..equation_count-1 denklem indeksleri atar.
    void finalize();

    int equation_count() const { return equation_count_; }

    // Equation index for a free DOF (>=0); -1 for a fixed DOF.
    int equation(int global_dof) const { return equation_[global_dof]; }
    bool is_fixed(int global_dof) const { return equation_[global_dof] < 0; }

private:
    int dofs_per_node_;
    int node_count_;
    int total_dofs_;
    std::vector<char> fixed_;     // total_dofs_: 1 = sabit
    std::vector<int> equation_;   // total_dofs_: denklem indeksi ya da -1
    int equation_count_ = 0;
    bool finalized_ = false;
};

} // namespace katai::core
