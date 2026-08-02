# katai -- the command-line front end (layer 4, a composition root). It programs
# against the published facade only (enforced by check_architecture's `frontend`
# check) and, being a whole program, selects the solver backend by linking
# katai::linsolve -- main.cpp names no backend.
add_executable(katai_cli main.cpp)
set_target_properties(katai_cli PROPERTIES
  OUTPUT_NAME katai
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
target_link_libraries(katai_cli PRIVATE katai::api katai::linsolve)
katai_use_pch(katai_cli)
