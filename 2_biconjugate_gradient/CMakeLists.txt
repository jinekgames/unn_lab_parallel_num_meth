project(PARALLEL_NUM_METH_1_LU_DECOMPOSE
        VERSION 0.1
        DESCRIPTION "Parallel Numerical Methods. Lab 1. LU-decomposition"
        LANGUAGES CXX)

add_library(2_implementation SHARED ./implementation.cpp)
add_executable(2_test test.cpp)
target_link_libraries(2_test PUBLIC 2_implementation)
