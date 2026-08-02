# katai_linsolve — sparse direct solvers behind one interface.
#
# The backend is chosen here, by compiling exactly one selection translation unit.
# No source file in the tree contains a preprocessor conditional on the backend:
# an optional capability is a separate target, not an #ifdef.
#
# The Eigen implementation is always compiled, even when PARDISO is selected, so a
# test can compare the two backends on the same system rather than trusting either
# on its own.

set(KATAI_LINSOLVE_SOURCES src/eigen_direct_solver.cpp)

if(KATAI_WITH_MKL AND MKL_FOUND)
  list(APPEND KATAI_LINSOLVE_SOURCES src/select_pardiso.cpp)
  set(KATAI_LINSOLVE_BACKEND "pardiso (Eigen also available for comparison)")
else()
  list(APPEND KATAI_LINSOLVE_SOURCES src/select_eigen.cpp)
  set(KATAI_LINSOLVE_BACKEND "eigen")
endif()

add_library(katai_linsolve STATIC ${KATAI_LINSOLVE_SOURCES})
target_include_directories(katai_linsolve PUBLIC include)
target_link_libraries(katai_linsolve PUBLIC katai::math)
add_library(katai::linsolve ALIAS katai_linsolve)
katai_use_pch(katai_linsolve)

message(STATUS "katai_linsolve backend: ${KATAI_LINSOLVE_BACKEND}")
