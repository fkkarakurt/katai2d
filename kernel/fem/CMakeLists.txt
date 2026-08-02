# katai_fem -- the finite-element system: element kernels plus DOF management
# and global assembly. Elements and assembly are one module on purpose: the
# mutual reference embedded_beam <-> dof_map is a legitimate detail inside the
# module, where between two modules it would be a cycle.
add_library(katai_fem STATIC
  src/elements/tri6.cpp
  src/elements/tri15.cpp
  src/assembly/dof_map.cpp
  src/assembly/assembler.cpp)
target_include_directories(katai_fem PUBLIC include)
target_link_libraries(katai_fem PUBLIC katai::materials katai::math katai::mesh)
add_library(katai::fem ALIAS katai_fem)
katai_use_pch(katai_fem)
