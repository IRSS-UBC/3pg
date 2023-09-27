# CMakeLists.txt
add_definitions(-std=c++20)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
cmake_minimum_required(VERSION 3.26)

# set(Boost_USE_STATIC_LIBS ON)
# set(DVCPKG_TARGET_TRIPLET x64-windows-static)
project(3pg)

find_package(GDAL CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options) # algorithm and program_options are header-only, doesn't need to include a COMPONENTS section
if ( Boost_FOUND )
    include_directories( ${Boost_INCLUDE_DIRS} )
else()
    message( FATAL_ERROR "Required Boost packages not found. Perhaps add -DBOOST_ROOT?" )
endif()


include_directories(include)

set(SOURCE_FILES
	3pg_sc/3pg_src/3pg.cpp
	3pg_sc/3pg_src/util.cpp
	3pg_sc/3pg_src/The_3PG_model.cpp
	3pg_sc/3pg_src/Data_io.cpp
	3pg_sc/3pg_src/MYDate.cpp
    3pg_sc/3pg_src/GDALRasterImage.cpp
)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(3pg PRIVATE GDAL::GDAL Boost::boost Boost::program_options)