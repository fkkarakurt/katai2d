# katai_math -- math kernel. Sparse CSR matrix + the PARDISO adapter.
set(KATAI_MATH_SOURCES
  src/katai_math.cpp
  src/sparse_matrix.cpp
  src/thread_pool.cpp)

# The PARDISO adapter is compiled only when MKL is present (decision D6). Nothing
# above this line asks whether it was: the backend is reached through
# katai::linsolve, whose selection translation unit is chosen by CMake. That is why
# no KATAI_WITH_MKL compile definition is exported any more -- a capability that
# shared code can test for is a capability that multiplies the number of programs
# that exist while only one of them gets tested.
if(KATAI_WITH_MKL AND MKL_FOUND)
  list(APPEND KATAI_MATH_SOURCES src/pardiso_solver.cpp)
endif()

add_library(katai_math STATIC ${KATAI_MATH_SOURCES})
target_include_directories(katai_math PUBLIC include)
target_link_libraries(katai_math PUBLIC katai::eigen)

# PUBLIC so the adapter and the two PARDISO adapter tests get MKL's headers and
# symbols transitively.
if(KATAI_WITH_MKL AND MKL_FOUND)
  target_link_libraries(katai_math PUBLIC MKL::MKL)
endif()

add_library(katai::math ALIAS katai_math)

# The shared precompiled header lives here because every target in the tree links
# katai::math, so every one of them can reuse this one. See the root CMakeLists for
# the measurement, and for why nothing of ours is allowed in it.
#
# Reuse is what makes this pay. A target with a single source file that builds its own
# precompiled header has done the same parsing twice and gained nothing; the saving
# comes entirely from ~130 targets sharing one.
if(KATAI_PCH)
  target_precompile_headers(katai_math PRIVATE
    <Eigen/Core>
    <Eigen/Dense>
    <algorithm>
    <array>
    <cmath>
    <cstdio>
    <functional>
    <memory>
    <stdexcept>
    <string>
    <vector>)
endif()

# Give a target the shared precompiled header. Kept as a function so the reason lives
# in one place and call sites stay a single line.
function(katai_use_pch target)
  if(KATAI_PCH)
    target_precompile_headers(${target} REUSE_FROM katai_math)
  endif()
endfunction()
