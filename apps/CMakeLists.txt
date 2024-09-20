find_package(GDAL CONFIG REQUIRED)

set(SOURCE_FILES
	Data_io.cpp
    	Data_io.hpp
	DataInput.cpp
	DataInput.hpp
	DataOutput.cpp
	DataOutput.hpp
	GDALRasterImage.cpp 
	GDALRasterImage.hpp
	ParamStructs.hpp
	The_3PG_model.cpp
    	The_3PG_model.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} 3pg.cpp )
target_link_libraries(${PROJECT_NAME} PRIVATE GDAL::GDAL)

#statically link
set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
target_compile_options(${PROJECT_NAME} PRIVATE /MT)
target_link_options(${PROJECT_NAME} PRIVATE /INCREMENTAL:NO /NODEFAULTLIB:MSVCRT)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	add_library(3pg2Test ${SOURCE_FILES})
	target_link_libraries(3pg2Test PRIVATE GDAL::GDAL)
endif()
