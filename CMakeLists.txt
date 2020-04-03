cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set the project name
project(qmc VERSION 3.0)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/googletest)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/eigen)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCES  src/geometry.cpp src/wavefunction/jastrows/jastrowSquareWell.cpp src/wavefunction/wavefunction.cpp src/wavefunction/productWavefunction.cpp src/tableDistances.cpp src/potential.cpp src/walkers.cpp src/initializer.cpp src/energy.cpp)


add_executable(main src/main.cpp ${SOURCES})

target_include_directories(main PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           ${PROJECT_BINARY_DIR}/../external/eigen)
target_link_libraries (main Eigen3::Eigen)

set(TEST_SOURCES  tst/parameters-test.cpp tst/jastrow-test.cpp tst/wavefunction-test.cpp)
add_executable(tst tst/main.cpp ${TEST_SOURCES} ${SOURCES})
target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           )
target_link_libraries(tst PUBLIC ${PROJECT_BINARY_DIR}/../external gtest Eigen3::Eigen)
install(TARGETS main DESTINATION /home/luca/source/qmc3)
install(TARGETS tst DESTINATION /home/luca/source/qmc3)

