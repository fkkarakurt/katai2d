# katai_api -- layer 3, the published facade: the one surface front ends (CLI, GUI,
# Python) program against (Stage D; the section 7.3 boundary). Header-only (INTERFACE).
add_library(katai_api INTERFACE)
target_include_directories(katai_api INTERFACE include)
target_link_libraries(katai_api INTERFACE katai::io katai::jobs katai::model)
add_library(katai::api ALIAS katai_api)
