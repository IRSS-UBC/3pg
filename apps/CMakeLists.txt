find_package(GDAL CONFIG REQUIRED)

set(SOURCE_FILES
    	3pg.cpp
	Data_io.cpp
    	Data_io.hpp
	DataInput.cpp
		DataInput.hpp
	DataOutput.cpp
		DataOutput.hpp
	GDALRasterImage.cpp
    	GDALRasterImage.hpp
	MYDate.cpp
    	MYDate.h
	ParamStructs.hpp
	The_3PG_model.cpp
    	The_3PG_model.hpp
	util.cpp
    	util.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
add_library(${PROJECT_NAME}_lib STATIC ${SOURCE_FILES})
target_link_libraries(3pg2 PRIVATE GDAL::GDAL)