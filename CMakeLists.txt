# KATAI 2D — root CMake
# Phase 0 / P0.1: the build scaffold. Build speed is a KPI (see docs/DECISIONS D15):
#   Ninja + ccache + MSVC embedded debug info (/Z7, cacheable) + modular targets.
cmake_minimum_required(VERSION 3.25)

project(KATAI2D
  VERSION 0.0.1
  DESCRIPTION "KATAI 2D — Geotechnical 2D Finite Element Analysis"
  LANGUAGES C CXX)

# --- Build KPI: ccache launcher -------------------------------------------
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "ccache found: ${CCACHE_PROGRAM}")
  set(CMAKE_C_COMPILER_LAUNCHER   "${CCACHE_PROGRAM}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# --- Build KPI: parallel job pools (memory-bound compilation) --------------
# If too many heavy TUs (build_problem.hpp, ~38 template-heavy tests) compile at
# once, RAM runs out → MSVC C1060 "out of heap". Pin concurrent COMPILATION to
# the core count (a safe ceiling on a 4c/16GB box), keep LINK in its own pool.
# That way `cmake --build` (Ninja default -j) runs safely+fast without a manual
# `-j2` (2→4 concurrent compiles = ~2× speedup). KPI: docs/DECISIONS D15.
if(CMAKE_GENERATOR MATCHES "Ninja")
  cmake_host_system_information(RESULT _ncore QUERY NUMBER_OF_PHYSICAL_CORES)
  if(NOT _ncore OR _ncore LESS 2)
    set(_ncore 4)
  endif()
  set_property(GLOBAL PROPERTY JOB_POOLS compile_pool=${_ncore} link_pool=3)
  set(CMAKE_JOB_POOL_COMPILE compile_pool)
  set(CMAKE_JOB_POOL_LINK    link_pool)
  message(STATUS "Ninja job pools: compile=${_ncore}, link=3 (C1060 guard + parallelism)")
endif()

# --- Source encoding -------------------------------------------------------
# Every source file in this tree is UTF-8. Without /utf-8 MSVC reads them in the
# system code page instead, which silently mangles any non-ASCII text: several files
# once carried comments stored back mis-decoded (the mojibake later fixed in the
# language sweep). Setting it once for the whole build is both the fix for that and the
# reason the shared precompiled header works -- /utf-8 sets the source and execution
# character sets, and MSVC refuses a precompiled header built with different ones
# (C2855), so it cannot be a per-target flag.
if(MSVC)
  add_compile_options(/utf-8)
endif()

# --- Build KPI: precompiled header ----------------------------------------
# Measured 2026-07-28. A test translation unit parses 179,118 lines of Eigen from 422
# files, against 255 to 13,029 lines of this project's own headers -- Eigen is 93% to
# 99.9% of the parse volume, and it is byte-identical in every unit. A full build
# therefore re-parsed roughly 21 million lines of it, 120 times over, to produce the
# same result every time.
#
# This also corrects a claim recorded earlier, that the hour-long build is caused by
# core being header-only and by build_problem.hpp being re-parsed by its consumers.
# Those are real architectural defects, but they are about 7% of the parse volume.
# The build cost is Eigen.
#
# Only third-party and standard headers may go in the precompiled header. Nothing
# from this project does, and that restriction is the whole design: editing our own
# headers -- which is what development consists of -- must never invalidate it.
option(KATAI_PCH "Precompile the third-party headers every translation unit shares" ON)

# --- Language standard -----------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- MSVC: ccache-friendly embedded debug info (/Z7), instead of /Zi (CMP0141) ---
# /Zi writes a separate PDB and does not ccache; /Z7 embeds into the .obj → cache hit.
if(POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
      "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")
endif()

# Default build type on single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

# Collect outputs in one directory (DLLs/exes side by side → easy to run)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Intel oneMKL (PARDISO + BLAS/LAPACK) -- decision D6 -------------------
# Selects the PARDISO backend for katai::linsolve. It is genuinely optional: with
# this OFF the tree builds and the whole suite runs on the vendored Eigen backend,
# which is what makes the verification record reproducible without a proprietary
# component.
#
# We deliberately do not depend on the top-level setvars.bat -- it is broken for
# installation paths containing spaces ("Program Files (x86)") because it does not
# quote its vars.bat calls. build.ps1 derives MKL_DIR and passes it with -D; the
# ONEAPI_ROOT fallback below is for using a preset or IDE directly.
option(KATAI_WITH_MKL "Intel oneMKL (PARDISO) solver backend" ON)
if(KATAI_WITH_MKL)
  if(NOT MKL_DIR AND DEFINED ENV{ONEAPI_ROOT})
    set(MKL_DIR "$ENV{ONEAPI_ROOT}/mkl/latest/lib/cmake/mkl")
  endif()
  # MKL_THREADING=tbb_thread would parallelize PARDISO later.
  #
  # MKL_LINK defaults to dynamic, and that is a measurement rather than a taste.
  # Static linking copies MKL into every executable that solves: measured 2026-07-28,
  # test_solve.exe was 130.7 MB and took 31.4 s to link, against 4.6 MB and 5.0 s
  # dynamically -- 6x the link time and 28x the size, repeated across the ~76 tests
  # that use the solver. Linking was 57 of the ~175 minutes a full build costs.
  #
  # The price is that the MKL runtime has to be found at load time, and it cannot be
  # solved by copying a DLL or two: the executable imports only mkl_sequential, which
  # then loads mkl_core and a CPU-specific kernel (mkl_avx512, mkl_avx2, ...) at run
  # time, each of them 45-72 MB. So the directory goes on PATH instead -- for the test
  # suite through CTest's own environment, so it works from any shell.
  #
  # Set MKL_LINK=static for a build that has to run somewhere MKL is not installed.
  set(MKL_INTERFACE lp64       CACHE STRING "MKL integer interface (lp64/ilp64)")
  set(MKL_THREADING sequential CACHE STRING "MKL threading model")
  set(MKL_LINK      dynamic    CACHE STRING "MKL link mode (static/dynamic)")
  find_package(MKL CONFIG)
  if(MKL_FOUND)
    message(STATUS "oneMKL found: ${MKL_ROOT} (link: ${MKL_LINK})")
    if(MKL_LINK STREQUAL "dynamic")
      set(KATAI_MKL_RUNTIME_DIR "${MKL_ROOT}/bin")
    endif()
  else()
    message(WARNING "oneMKL not found although KATAI_WITH_MKL=ON -- falling back to "
                    "the Eigen backend. Pass MKL_DIR via build.ps1 if that is not intended.")
    set(KATAI_WITH_MKL OFF)
  endif()
endif()

# --- Composition -----------------------------------------------------------
# What the program is made of is decided here and nowhere else: shared code
# never tests a KATAI_WITH_* macro (architecture gate, check 'config'). The
# 'engine' preset turns this off to prove the physics needs no project schema.
option(KATAI_WITH_MODEL "Build the project schema module and its consumers" ON)

add_subdirectory(third_party)

# kernel/ -- the engine: everything needed to run and reproduce an analysis.
add_subdirectory(kernel/math)
add_subdirectory(kernel/linsolve)
add_subdirectory(kernel/geometry)
add_subdirectory(kernel/mesh)
add_subdirectory(kernel/materials)
add_subdirectory(kernel/fem)
add_subdirectory(kernel/analysis)
if(KATAI_WITH_MODEL)
  add_subdirectory(kernel/model)
  add_subdirectory(kernel/io)   # file formats: needs the schema (and the engine's result types)
  add_subdirectory(kernel/jobs) # layer 2: the analysis/mesh/flow drivers over the schema (Stage D)
  add_subdirectory(kernel/api)  # layer 3: the published facade front ends program against (Stage D)
endif()

# cli/ -- the command-line front end (needs the schema + facade).
option(KATAI_WITH_CLI "Build the katai command-line front end" ON)
if(KATAI_WITH_CLI AND KATAI_WITH_MODEL)
  add_subdirectory(cli)
endif()

# Tests are part of every development build. The option exists for consumers:
# a tree that pulls this project in via add_subdirectory or FetchContent wants
# the libraries and front ends, not a hundred test targets -- and in a static-MKL
# composition every test executable links the whole solver (~130 MB each), so
# building them is far from free.
option(KATAI_WITH_TESTS "Build the test suite" ON)
enable_testing()
if(KATAI_WITH_TESTS)
  add_subdirectory(tests)
endif()

# python/ -- the scripting front end (needs the schema + facade). Skips itself,
# with a message, when the Python development environment or nanobind is absent:
# an optional module is refused, never imitated (composability rule 3).
option(KATAI_WITH_PYTHON "Build the katai Python bindings" ON)
if(KATAI_WITH_PYTHON AND KATAI_WITH_MODEL)
  add_subdirectory(python)
endif()
