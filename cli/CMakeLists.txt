# katai -- the command-line front end (layer 4, a composition root). It programs
# against the published facade only (enforced by check_architecture's `frontend`
# check) and, being a whole program, selects the solver backend by linking
# katai::linsolve -- main.cpp names no backend.
# Binary identity (Windows). The version fields come from <katai/api/version.hpp>, the one
# source of truth, so the executable's Properties tab cannot drift from what `katai info`
# prints. An unsigned binary that also carries NO identity is the easiest possible target for
# a generic reputation heuristic; this is the part of that problem which costs nothing to fix.
set(KATAI_CLI_SOURCES main.cpp)
if(WIN32 AND MSVC)
  file(READ "${PROJECT_SOURCE_DIR}/kernel/api/include/katai/api/version.hpp" _katai_version_hpp)
  string(REGEX MATCH "kVersion = \"([^\"]+)\"" _katai_ver_match "${_katai_version_hpp}")
  set(KATAI_VERSION "${CMAKE_MATCH_1}")
  if(NOT KATAI_VERSION)
    message(FATAL_ERROR "cli: kVersion not found in version.hpp")
  endif()
  # "0.7.0" -> 0;7;0. A pre-release suffix ("-dev") stays out of the numeric fields, which
  # VERSIONINFO defines as integers, while the string fields keep the version verbatim.
  string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" _katai_ver_numeric "${KATAI_VERSION}")
  string(REPLACE "." ";" _katai_ver_parts "${_katai_ver_numeric}")
  list(GET _katai_ver_parts 0 KATAI_VER_MAJOR)
  list(GET _katai_ver_parts 1 KATAI_VER_MINOR)
  list(GET _katai_ver_parts 2 KATAI_VER_PATCH)
  configure_file(katai.rc.in "${CMAKE_CURRENT_BINARY_DIR}/katai.rc" @ONLY)
  list(APPEND KATAI_CLI_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/katai.rc")
endif()

add_executable(katai_cli ${KATAI_CLI_SOURCES})
set_target_properties(katai_cli PROPERTIES
  OUTPUT_NAME katai
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
target_link_libraries(katai_cli PRIVATE katai::api katai::linsolve)
katai_use_pch(katai_cli)
