# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.8)
set (CMAKE_CXX_STANDARD 17)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project ("3pg2")

# Include sub-projects.
include_directories("include")
add_subdirectory ("apps")