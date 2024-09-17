project(PARALLEL_NUM_METH_1_LU_DECOMPOSE
        VERSION 0.1
        DESCRIPTION "Parallel Numerical Methods. Lab 1. LU-decomposition"
        LANGUAGES CXX)

add_library(1_implementation STATIC ./implementation.cpp)
add_executable(1_test test.cpp)
target_link_libraries(1_test PUBLIC 1_implementation)
