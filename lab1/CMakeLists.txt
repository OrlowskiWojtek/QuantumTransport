cmake_minimum_required(VERSION 3.10)
project(QuantumDotSpectr VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(QuantumDotSpectr main.cpp utils/src/solver.cpp 
               utils/src/basis_funcs.cpp
               utils/src/params.cpp)
target_link_libraries (QuantumDotSpectr Eigen3::Eigen)

target_include_directories(QuantumDotSpectr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
