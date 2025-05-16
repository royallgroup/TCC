include_directories(${PROJECT_SOURCE_DIR}/inc/)
include_directories(${PROJECT_SOURCE_DIR}/inc/clusters/)
include_directories(${PROJECT_SOURCE_DIR}/extern/iniparser/inc/)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.c")

add_executable(tcc ${SOURCES})
target_link_libraries(tcc iniparser)
target_link_libraries(tcc m)

install(TARGETS tcc DESTINATION ${CMAKE_INSTALL_PREFIX})