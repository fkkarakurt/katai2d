# katai_add_test -- how a test becomes an executable and one or more CTest entries.
#
# This lives outside tests/ for one reason: a composition may want the helper without wanting
# the core's own hundred-odd test executables. The studio superbuild is exactly that -- it
# consumes this core as a subproject and adds its own tests, but the core suite has already run
# in the core's own build tree minutes earlier, and compiling it a second time was measured at
# 7587 s of the studio build's 11308 s. So the root includes this file unconditionally and adds
# tests/ only when KATAI_WITH_CORE_TESTS is on.
#
# Adding a test: katai_add_test(<name> [link targets...]). The default link is katai::math.
# katai::linsolve is added to every test unconditionally, because a test executable is a
# composition root: it is a whole program, and a program that solves has to contain a backend.
# Which backend that is comes from the configuration, never from the test -- see
# kernel/linsolve/CMakeLists.txt.
#
# The named modules do two more jobs. They label the test, so `ctest -L materials` runs exactly
# the tests whose subject is that module. And they scope it: in a composition that does not build
# one of the named modules the test is skipped, never linked against a substitute -- a missing
# module is refused, not imitated (composability rule 3).
function(katai_add_test name)
  # A composition may want only part of the core suite -- see KATAI_CORE_TESTS_MATCHING in the
  # root. The filter is a directory-scoped variable set by tests/CMakeLists.txt, so it narrows
  # the core's own tests and never a superbuild's.
  if(KATAI_TEST_FILTER AND NOT name MATCHES "${KATAI_TEST_FILTER}")
    set_property(GLOBAL APPEND PROPERTY KATAI_FILTERED_TESTS ${name})
    return()
  endif()
  set(libs ${ARGN})
  if(NOT libs)
    set(libs katai::math)
  endif()
  set(labels "")
  foreach(lib ${libs})
    if(lib MATCHES "^katai::(.+)$")
      if(NOT TARGET ${lib})
        set_property(GLOBAL APPEND PROPERTY KATAI_SKIPPED_TESTS ${name})
        return()
      endif()
      list(APPEND labels ${CMAKE_MATCH_1})
    endif()
  endforeach()
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ${libs} katai::linsolve)
  katai_use_pch(${name})
  # ONE EXECUTABLE, OPTIONALLY SEVERAL CTEST ENTRIES. Set KATAI_TEST_CASES before the call and
  # each case becomes `<name>_<case>` running `<name> <case>`. This exists for wall clock, not
  # for tidiness: with `ctest -j N` the suite cannot finish faster than its single longest test,
  # so a file holding several independent verification cases sets the floor for everything else.
  # Splitting it lets the cases run side by side. The executable still runs all of them when
  # invoked with no argument, so a developer keeps the whole story in one command.
  if(KATAI_TEST_CASES)
    foreach(case ${KATAI_TEST_CASES})
      add_test(NAME ${name}_${case} COMMAND ${name} ${case})
      list(APPEND _katai_entries ${name}_${case})
    endforeach()
  else()
    add_test(NAME ${name} COMMAND ${name})
    set(_katai_entries ${name})
  endif()
  set_tests_properties(${_katai_entries} PROPERTIES LABELS "${labels}")
  # With MKL linked dynamically its runtime has to be locatable at load time. Putting
  # it in the test's own environment rather than expecting the caller's PATH means
  # `ctest` works from any shell, including a bare one and CI.
  if(KATAI_MKL_RUNTIME_DIR)
    set_tests_properties(${_katai_entries} PROPERTIES
      ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${KATAI_MKL_RUNTIME_DIR}")
  endif()
endfunction()
