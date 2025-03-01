# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

option(BUILD_SHARED_LIBS "build with dynamic libraries except for the standalone executable" ON)
option(BUILD_GOOGLE_TESTS "build 3pg with tests included" OFF)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project("3pg")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add include directories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(NOT BUILD_SHARED_LIBS)
	include_directories($ENV{GDAL_STATIC_INCLUDE_DIR})
endif()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add link directories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(NOT BUILD_SHARED_LIBS)
	link_directories($ENV{GDAL_STATIC_LIB_DIR})
endif()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add subdirectories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_subdirectory("apps")
if(BUILD_GOOGLE_TESTS)
	add_subdirectory("tests")
endif()
