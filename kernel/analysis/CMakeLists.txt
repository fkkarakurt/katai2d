# katai_analysis -- solution procedures and post-processing over the FEM
# system: nonlinear and dynamic solvers, strength reduction, seepage, stress
# recovery. The linear solver enters as a callback (LinearSolve), so the
# backend is chosen by the composition root that builds the callback -- this
# module neither names nor links one.
add_library(katai_analysis STATIC
  src/consolidation.cpp
  src/coupled_flow_deformation.cpp
  src/structural_dynamics.cpp
  src/transient_flow.cpp
  src/nonlinear_solver.cpp
  src/dynamics_nonlinear.cpp
  src/strength_reduction.cpp
  src/seepage.cpp
  src/post/stress_recovery.cpp)
target_include_directories(katai_analysis PUBLIC include)
target_link_libraries(katai_analysis PUBLIC
  katai::fem katai::materials katai::math katai::mesh)
add_library(katai::analysis ALIAS katai_analysis)
katai_use_pch(katai_analysis)
