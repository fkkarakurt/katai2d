# katai._core -- the nanobind extension binding the published facade (layer 4, a front
# end: core_bindings.cpp includes <katai/api/katai.hpp> and nothing else of the engine;
# check_architecture's `frontend` check enforces it).
#
# Optional module (composability rule 3): without a Python development environment or
# nanobind the module is SKIPPED with a message, never imitated. nanobind is a
# build-time dependency installed into the build Python (pip install nanobind), not
# vendored -- it ships its own CMake config.
#
# Distribution decision (2026-08-02): pip wheel -- the extension plus the pure-Python
# `katai` package run under the user's own Python. For the build tree, the package is
# staged under ${CMAKE_BINARY_DIR}/python/katai (extension + copied .py files) and the
# smoke test imports it from there.

# The tree uses ONE FindPython module everywhere (the NEW one, which nanobind
# requires); tests/ already found the Interpreter component, this call adds the
# Development.Module artifacts to the same cached interpreter. A stable-ABI
# build additionally needs the SABIModule artifacts -- without that component
# nanobind silently falls back to a version-specific build (measured: the
# first wheel came out cp314-tagged with STABLE_ABI=ON).
# The wheel preset builds against the CPython stable ABI (3.12+): one wheel
# serves every current Python instead of one per minor version. The dev builds
# keep the full (faster, version-specific) ABI.
option(KATAI_PYTHON_STABLE_ABI
       "Build _core against the CPython stable ABI (one abi3 wheel, Python >= 3.12)" OFF)

if(KATAI_PYTHON_STABLE_ABI)
  find_package(Python 3.12 COMPONENTS Interpreter Development.Module Development.SABIModule QUIET)
else()
  find_package(Python 3.10 COMPONENTS Interpreter Development.Module QUIET)
endif()
if(NOT Python_FOUND)
  message(STATUS "katai python: no Python development environment; skipping the bindings")
  return()
endif()

execute_process(
  COMMAND ${Python_EXECUTABLE} -m nanobind --cmake_dir
  OUTPUT_VARIABLE KATAI_NANOBIND_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT KATAI_NANOBIND_DIR)
  message(STATUS "katai python: nanobind is not installed in ${Python_EXECUTABLE}; skipping")
  return()
endif()
find_package(nanobind CONFIG REQUIRED PATHS ${KATAI_NANOBIND_DIR} NO_DEFAULT_PATH)

if(KATAI_PYTHON_STABLE_ABI)
  nanobind_add_module(katai_python_core NB_STATIC STABLE_ABI src/core_bindings.cpp)
else()
  nanobind_add_module(katai_python_core NB_STATIC src/core_bindings.cpp)
endif()
set_target_properties(katai_python_core PROPERTIES
  OUTPUT_NAME "_core"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/katai
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/katai)
target_link_libraries(katai_python_core PRIVATE katai::api katai::linsolve)

# Stage the pure-Python half of the package beside the extension (never build INTO the
# source tree). configure_file, not file(COPY): COPY runs once at configure time and a
# later edit leaves the staged copy STALE (measured -- the wheel script read an old
# __init__.py); configure_file registers the source as a configure dependency, so
# editing a .py re-stages it on the next build.
file(GLOB KATAI_PY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/katai/*.py)
foreach(katai_py ${KATAI_PY_SOURCES})
  get_filename_component(katai_py_name ${katai_py} NAME)
  configure_file(${katai_py} ${CMAKE_BINARY_DIR}/python/katai/${katai_py_name} COPYONLY)
endforeach()

# The extension is opaque to IDEs without a stub; nanobind's stubgen derives
# _core.pyi from the built module (signatures, IntEnum members, docstrings), so
# completion works in any editor -- part of the easy-to-learn rule. Generated
# beside the extension, never committed; test_python_stubs pins completeness.
add_custom_command(TARGET katai_python_core POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E env
          "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
          "KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}"
          "PYTHONDONTWRITEBYTECODE=1"
          ${Python_EXECUTABLE} -m nanobind.stubgen -q -m katai._core
          -o ${CMAKE_BINARY_DIR}/python/katai/_core.pyi
  COMMENT "Generating katai/_core.pyi type stub")

# The smoke test: import the package, express a corpus case, run it, be refused
# honestly. Needs the staged package on PYTHONPATH and the MKL runtime on PATH.
add_test(NAME test_python_smoke
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.py
                 ${PROJECT_SOURCE_DIR}/tests/corpus)
set_tests_properties(test_python_smoke PROPERTIES
  LABELS "python;api"
  ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}")

# The engineer-facing DSL's acceptance test: every corpus case built through the
# easy surface is BYTE-identical to its checked-in .k2d, and the slope case runs
# end to end with the same FoS as the file-loaded run.
add_test(NAME test_python_dsl
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/dsl_corpus.py
                 ${PROJECT_SOURCE_DIR}/tests/corpus)
set_tests_properties(test_python_dsl PROPERTIES
  LABELS "python;api"
  ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}")

# Structural elements from the script -- the layer that used to stop exactly where
# real engineering starts -- plus the summary layer that makes a result readable.
# Checked by SOLVING against the checked-in KV-STR-005 rather than by byte-identity;
# the file's own builder writes redundant all-active vectors, so the two describe the
# same run and only the ANSWER can be compared bit for bit. The test says why.
add_test(NAME test_python_structures
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/structures.py
                 ${PROJECT_SOURCE_DIR}/tests/corpus)
set_tests_properties(test_python_structures PROPERTIES
  LABELS "python;api;model"
  ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}")

# The schema-coverage gate: every key the .k2d writer emits must be reachable,
# faithful and writable through the bindings. The key universe is extracted by
# importing check_k2d_spec.py itself, so this gate and the format gate can never
# disagree about what the format contains.
add_test(NAME test_python_schema_coverage
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/schema_coverage.py
                 ${PROJECT_SOURCE_DIR})
# PYTHONDONTWRITEBYTECODE: this test imports check_k2d_spec.py from scripts/;
# without it, Python drops a __pycache__ into the SOURCE tree.
set_tests_properties(test_python_schema_coverage PROPERTIES
  LABELS "python;api;io"
  ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR};PYTHONDONTWRITEBYTECODE=1")

# The product-text gate: no other FE program is named in anything a user reads --
# the command line's output and the package's docstrings here, the Studio's labels
# and reports on its own side. Source COMMENTS and the validation record are out of
# scope on purpose: naming the program a number was compared against is the point of
# a comparison. A rule nobody was there for is a rule that comes back under the next
# tooltip, so it is a test and not an edit.
add_test(NAME test_product_text
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/product_text.py
                 ${PROJECT_SOURCE_DIR}/cli ${CMAKE_CURRENT_SOURCE_DIR}/katai)
set_tests_properties(test_product_text PROPERTIES LABELS "python;docs")

# The stub-completeness pin: every public _core name must appear in the
# generated _core.pyi (see the stubgen POST_BUILD step above).
add_test(NAME test_python_stubs
         COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/stub_check.py)
set_tests_properties(test_python_stubs PROPERTIES
  LABELS "python;api"
  ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}")

# The examples are executable documentation and run as tests: each solves its
# problem and refuses (non-zero exit) if the result leaves the corpus band the
# same case is pinned to. An example that drifts from the API or the physics
# fails the suite, not the first reader who tries it.
foreach(example slope_stability consolidation_terzaghi staged_excavation anchored_excavation)
  add_test(NAME test_python_example_${example}
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples/${example}.py)
  set_tests_properties(test_python_example_${example} PROPERTIES
    LABELS "python;api"
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;KATAI_DLL_PATH=${KATAI_MKL_RUNTIME_DIR}")
endforeach()
