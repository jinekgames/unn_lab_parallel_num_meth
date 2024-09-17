cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 14)

project(PARALLEL_NUM_METH
        VERSION 0.1
        DESCRIPTION "Parallel Numerical Methods"
        LANGUAGES CXX)


option(BUILD_LAB_1 ON)
option(BUILD_LAB_2 OFF)

message("Build 1 Lab: ${BUILD_LAB_1}")
message("Build 2 Lab: ${BUILD_LAB_2}")


if(${BUILD_LAB_1})
    add_subdirectory(1_lu_decomposition)
endif()

if(${BUILD_LAB_2})
    add_subdirectory(2_biconjugate_gradient)
endif()
