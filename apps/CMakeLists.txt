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

if(BUILD_GOOGLE_TESTS)
	#if running tests, 3pg should .lib file which cmake can find when generating test executables
	add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES})
else()
	#otherwise, it should be an executable
	add_executable(${PROJECT_NAME} ${SOURCE_FILES} main.cpp)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VS_GLOBAL_VcpkgEnabled false)

if(WIN32)
	if(BUILD_SHARED_LIBS)
		#set to compile with multi-threaded dynamic debug
		target_compile_options(${PROJECT_NAME} PRIVATE /MDd)
		
		#link to GDAL after utilizing CMAKE's find_package()
		target_link_libraries(${PROJECT_NAME} PRIVATE GDAL::GDAL)
	else()
		#set to compile with multi-threaded (static)
		target_compile_options(${PROJECT_NAME} PRIVATE /MT)
		target_link_options(${PROJECT_NAME} PRIVATE /INCREMENTAL:NO /NODEFAULTLIB:MSVCRT)
		
		#link the static lib files required for GDAL
		target_link_libraries(${PROJECT_NAME} PRIVATE ${GDAL_STATIC_LIB_FILES})
		
		#the following are static libraries which are not part of GDAL but are required for GDAL on windows
		target_link_libraries(${PROJECT_NAME} PRIVATE
			crypt32
			wbemuuid
			Secur32
			ws2_32
			Normaliz
			Wldap32
		)
	endif()
else()
	# link Boost and Threads
	target_link_libraries(${PROJECT_NAME} PRIVATE GDAL::GDAL)
	find_package(Boost COMPONENTS system thread REQUIRED)
	find_package(Threads REQUIRED)
	target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system Boost::thread Threads::Threads)
endif()

# Debug mode
if(WITH_GOOGLE_TESTS)
	#add_library(3pgTest ${SOURCE_FILES})
	#target_link_libraries(3pgTest PRIVATE GDAL::GDAL)
	if(NOT MSVC)
		target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system Boost::thread Threads::Threads)
	endif()
endif()
