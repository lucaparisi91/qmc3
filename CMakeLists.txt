cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set the project name
project(qmc VERSION 3.0)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/googletest)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/eigen)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCES  src/geometry.cpp src/wavefunction/jastrows/jastrowSquareWell.cpp src/wavefunction/wavefunction.cpp src/wavefunction/productWavefunction.cpp src/tableDistances.cpp src/potential.cpp src/walkers.cpp src/initializer.cpp src/energy.cpp src/estimators.cpp src/driver.cpp src/estimatorCollection.cpp src/tools.cpp src/moves/vmcMoves.cpp src/moves/dmcMoves.cpp src/dmcDriver.cpp src/metropolis.cpp)

add_library(main_lib STATIC ${SOURCES})

target_include_directories(main_lib PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
			   )
target_include_directories(main_lib SYSTEM PUBLIC ${PROJECT_BINARY_DIR}/../external/eigen
		           )
			 
target_link_libraries (main_lib Eigen3::Eigen)
target_link_libraries (main_lib Eigen3::Eigen)

add_executable(main src/main.cpp)
target_link_libraries(main main_lib)


set(TEST_SOURCES  tst/parameters-test.cpp tst/jastrow-test.cpp tst/wavefunction-test.cpp)
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           )
target_link_libraries(tst PUBLIC ${PROJECT_BINARY_DIR}/../external gtest Eigen3::Eigen main_lib)			 
