# katai_model — GUI project data model (materials, soil layers, borehole). Header-only (INTERFACE).
add_library(katai_model INTERFACE)
target_include_directories(katai_model INTERFACE include)
add_library(katai::model ALIAS katai_model)
