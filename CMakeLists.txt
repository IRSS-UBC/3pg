# CMakeLists.txt
add_definitions(-std=c++20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.5)
project(3pg)

find_package(GDAL CONFIG REQUIRED)

include_directories(include)

set(SOURCE_FILES
	3pg_sc/3pg_src/3pg.cpp
	3pg_sc/3pg_src/util.cpp
	3pg_sc/3pg_src/The_3PG_model.cpp
	3pg_sc/3pg_src/Data_io.cpp
	3pg_sc/3pg_src/MYDate.cpp
	3pg_sc/3pg_src/getopt.c
    3pg_sc/3pg_src/GDALRasterImage.cpp
)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(3pg PRIVATE GDAL::GDAL)


