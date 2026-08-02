# katai_geometry — 2D CAD (minimal). P0.4: a rectangular domain with a material.
# Header-only for now → INTERFACE library.
add_library(katai_geometry INTERFACE)
target_include_directories(katai_geometry INTERFACE include)
add_library(katai::geometry ALIAS katai_geometry)
