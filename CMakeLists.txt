cmake_minimum_required(VERSION 3.6)
project(TCC)
set(CMAKE_C_STANDARD 90)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(iniparser/src)
add_subdirectory(tcc/src)
add_subdirectory(dtcc_postprocessing/src)