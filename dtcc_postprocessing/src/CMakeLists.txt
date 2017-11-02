include_directories(${PROJECT_SOURCE_DIR}/dtcc_postprocessing/inc/)
include_directories(${PROJECT_SOURCE_DIR}/iniparser/inc/)

add_executable(dtcc_post_analysis dynamic_analysis.c)
target_link_libraries(dtcc_post_analysis iniparser)