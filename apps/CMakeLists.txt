find_package(GDAL CONFIG REQUIRED)

set(SOURCE_FILES
	DataInput.cpp
	DataInput.hpp
	DataOutput.cpp
	DataOutput.hpp
	GDALRasterImage.cpp 
	GDALRasterImage.hpp
	ParamStructs.hpp
	3pgModel.cpp
    	3pgModel.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} main.cpp )
target_link_libraries(${PROJECT_NAME} PRIVATE GDAL::GDAL)

#statically link
set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
target_compile_options(${PROJECT_NAME} PRIVATE /MT)
target_link_options(${PROJECT_NAME} PRIVATE /INCREMENTAL:NO /NODEFAULTLIB:MSVCRT)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	add_library(3pg2Test ${SOURCE_FILES})
	target_link_libraries(3pg2Test PRIVATE GDAL::GDAL)
endif()
