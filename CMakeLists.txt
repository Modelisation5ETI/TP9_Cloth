#CMakeLists.txt to build project using cmake
#run $cmake path/to/CMakeLists.txt

cmake_minimum_required( VERSION 2.6 )

project( Cloth )

set( CMAKE_BUILD_TYPE Debug)

# Part 1 : Springs
add_subdirectory( "${CMAKE_SOURCE_DIR}/Part1_Spring" 
  "${CMAKE_BINARY_DIR}/Part1_Spring-${CMAKE_BUILD_TYPE}" )


