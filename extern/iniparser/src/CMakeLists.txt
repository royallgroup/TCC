include_directories(${PROJECT_SOURCE_DIR}/extern/iniparser/inc/)
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/extern/iniparser/src/*.c")

add_library(iniparser ${SOURCES})