#include <katai/fem/assembly/dof_map.hpp>

#include <cassert>

namespace katai::core {

DofMap::DofMap(int node_count, int dofs_per_node)
    : dofs_per_node_(dofs_per_node),
      node_count_(node_count),
      total_dofs_(node_count * dofs_per_node) {
    assert(node_count >= 0 && dofs_per_node > 0);
    fixed_.assign(total_dofs_, 0);
    equation_.assign(total_dofs_, -1);
}

int DofMap::add_extra_dof() {
    assert(!finalized_ && "add_extra_dof() must be called before finalize");
    const int g = total_dofs_++;
    fixed_.push_back(0);
    equation_.push_back(-1);
    return g;
}

void DofMap::fix(int global_dof) {
    assert(!finalized_ && "fix() must be called before finalize");
    assert(global_dof >= 0 && global_dof < total_dofs_);
    fixed_[global_dof] = 1;
}

void DofMap::finalize() {
    equation_count_ = 0;
    for (int d = 0; d < total_dofs_; ++d)
        equation_[d] = fixed_[d] ? -1 : equation_count_++;
    finalized_ = true;
}

} // namespace katai::core
