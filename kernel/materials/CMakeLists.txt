# katai_materials -- the constitutive library. Material models are pure
# functions over material state and plain arrays, so this module depends on
# nothing else in the tree (vendored Eigen only). That independence is load-
# bearing: it is what lets the registry resolve a model by name through one
# seam, and what will let a user material enter through the same seam from
# outside the tree (versioned ABI, roadmap section 8.1).
add_library(katai_materials STATIC
  src/mohr_coulomb.cpp
  src/registry.cpp)
target_include_directories(katai_materials PUBLIC include)
target_link_libraries(katai_materials PUBLIC katai::eigen)
add_library(katai::materials ALIAS katai_materials)
katai_use_pch(katai_materials)
