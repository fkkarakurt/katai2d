# katai_jobs -- layer 2, the execution seam: the analysis driver (schema -> neutral
# resolution -> engine phase strategies), the model-to-mesh builder, the groundwater-flow
# driver and the Job object. Every front end runs a project through this module, so a
# .k2d means the same run everywhere (Stage D).
#
# A COMPILED library on purpose (section 5.2): the driver bodies used to be
# header-inline, and every one of the ~26 driver-including test TUs paid ~220 s of
# codegen recompiling them. As library symbols the suite pays for them exactly once;
# the headers keep the types, the small helpers and the declarations.
add_library(katai_jobs STATIC
  src/driver.cpp
  src/flow_driver.cpp
  src/mesh_builder.cpp)
target_include_directories(katai_jobs PUBLIC include)
target_link_libraries(katai_jobs PUBLIC
  katai::analysis
  katai::fem
  katai::io
  katai::linsolve
  katai::materials
  katai::mesh
  katai::model)
katai_use_pch(katai_jobs)
add_library(katai::jobs ALIAS katai_jobs)
