# third_party -- vendored (drop-in) dependencies (decision T-ENV-4).
# Eigen (MPL-2.0): headers plus the upstream license files, nothing else.

# --- Eigen: header-only -> INTERFACE target --------------------------------
add_library(eigen INTERFACE)
target_include_directories(eigen SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
add_library(katai::eigen ALIAS eigen)
