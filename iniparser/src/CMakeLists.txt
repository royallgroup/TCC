include_directories(${PROJECT_SOURCE_DIR}/iniparser/inc/)
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/iniparser/src/*.c")

add_library(iniparser ${SOURCES})

