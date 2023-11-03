# CMakeLists.txt
add_definitions(-std=c++20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
cmake_minimum_required(VERSION 3.26)

set(Boost_USE_STATIC_LIBS ON)

project(3pg)

# Make sure to set the PROJ_LIB environment variable to the location of the proj.db file!
# see: https://github.com/OSGeo/gdal/issues/2248
# If PROJ_LIB isn't set properly before building you'll run into an error like: "proj_create_from_database: cannot find proj.db"
find_package(GDAL CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options) # algorithm and program_options are header-only, doesn't need to include a COMPONENTS section

# include_directories(${Boost_INCLUDE_DIRS})
# include_directories(${GDAL_INCLUDE_DIR})
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