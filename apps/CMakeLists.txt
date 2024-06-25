find_package(GDAL CONFIG REQUIRED)

set(SOURCE_FILES
    	3pg.cpp
	util.cpp
    	util.hpp
	The_3PG_model.cpp
    	The_3PG_model.hpp
	Data_io.cpp
    	Data_io.hpp
	MYDate.cpp
    	MYDate.h
    	GDALRasterImage.cpp
    	GDALRasterImage.hpp
	DataOutput.cpp
	DataOutput.hpp
    ParamStructs.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(3pg2 PRIVATE GDAL::GDAL)