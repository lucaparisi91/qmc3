cmake_minimum_required(VERSION 3.10)

# set the project name
project(qmc VERSION 3.0)

add_subdirectory(dmc)

enable_testing()
include(GoogleTest)
add_subdirectory(tst)
add_subdirectory(pimc)
add_subdirectory(external/particleKernels)
