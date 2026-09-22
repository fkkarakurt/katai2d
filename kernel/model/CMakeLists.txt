# katai_model — GUI project data model (materials, soil layers, borehole). Header-only (INTERFACE).
add_library(katai_model INTERFACE)
target_include_directories(katai_model INTERFACE include)
target_link_libraries(katai_model INTERFACE katai::geometry)   # region_edit.hpp
add_library(katai::model ALIAS katai_model)
