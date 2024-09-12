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
	util.cpp
    	util.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} 3pg.cpp )
target_link_libraries(${PROJECT_NAME} PRIVATE GDAL::GDAL)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	add_library(3pg2Test ${SOURCE_FILES} loggerForTest.cpp)
	target_link_libraries(3pg2Test PRIVATE GDAL::GDAL)
endif()
