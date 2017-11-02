include_directories(${PROJECT_SOURCE_DIR}/tcc/inc/)
include_directories(${PROJECT_SOURCE_DIR}/iniparser/inc/)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/tcc/src/*.c")

add_executable(tcc ${SOURCES})
target_link_libraries(tcc iniparser)