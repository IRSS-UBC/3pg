# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.11)
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

option(BUILD_SHARED_LIBS "build with dynamic libraries except for the standalone executable" ON)
option(BUILD_GOOGLE_TESTS "build 3pg with tests included" OFF)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project("3pg")

#if we're building a release version, add the speediest compiler optimizations
#if(NOT BUILD_SHARED_LIBS AND MSVC)
#	add_compile_options("-O2")
#endif()

if(BUILD_SHARED_LIBS)
	set(VCPKG_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg)
	
	set(GDAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/gdal_x64-windows/share/gdal)
	find_package(GDAL CONFIG REQUIRED)
endif()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add include directories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows/include/boost)
if(NOT BUILD_SHARED_LIBS)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/include)
endif()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add link directories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(NOT BUILD_SHARED_LIBS)
	#link_directories($ENV{GDAL_STATIC_LIB_DIR})
	link_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/lib)
endif()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add subdirectories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_subdirectory("apps")
if(BUILD_GOOGLE_TESTS)
	add_subdirectory("tests")
endif()
