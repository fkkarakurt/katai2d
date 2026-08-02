# katai_io -- the project (.k2d, versioned JSON) and results (.res, versioned
# binary) file formats, unified in one module so the two version disciplines
# live side by side (Stage C). Header-only (INTERFACE).
add_library(katai_io INTERFACE)
target_include_directories(katai_io INTERFACE include)
target_link_libraries(katai_io INTERFACE katai::model katai::analysis katai::mesh)
add_library(katai::io ALIAS katai_io)
