# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project("3pg")

# Include sub-projects.
include_directories("include")
add_subdirectory("apps")

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_subdirectory("tests")
endif()
