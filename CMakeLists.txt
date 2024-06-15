cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

project(PARALLEL_NUM_METH
        VERSION 0.1
        DESCRIPTION "Parallel Numerical Methods"
        LANGUAGES CXX)


set(BUILD_LAB_1 ON)


if(${BUILD_LAB_1})
    add_subdirectory(1_lu_decomposition)
endif()
