# katai_mesh — mesher (minimal). P0.5: rectangle → structured tri6.
# Full Delaunay + Ruppert refinement arrives in Phase 1.
add_library(katai_mesh STATIC
  src/mesh.cpp
  src/geom_predicates.cpp
  src/delaunay.cpp)
target_include_directories(katai_mesh PUBLIC include)
target_link_libraries(katai_mesh PUBLIC katai::geometry)
add_library(katai::mesh ALIAS katai_mesh)
katai_use_pch(katai_mesh)
