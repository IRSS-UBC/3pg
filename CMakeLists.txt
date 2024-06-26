# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required(VERSION 3.14)
set (CMAKE_CXX_STANDARD 17)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project (
	"3pg2"
	VERSION 0.1
	DESCRIPTION "A modernized 3PG model"
	LANGUAGES CXX)


include(FetchContent)
FetchContent_Declare(
  boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG        boost-1.85.0
)
FetchContent_MakeAvailable(boost)
# Adds Boost::boost

FetchContent_Declare(
  gdal
  GIT_REPOSITORY https://github.com/OSGeo/gdal.git
  GIT_TAG        v3.9.1
)
FetchContent_MakeAvailable(gdal)
# Adds GDAL::GDAL

# The executable code is here
add_subdirectory ("apps")