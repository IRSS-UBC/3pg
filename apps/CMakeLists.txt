IF(BUILD_SHARED_LIBS)
	find_package(GDAL CONFIG REQUIRED)
endif()	

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

add_executable(${PROJECT_NAME} ${SOURCE_FILES} main.cpp)

if(WIN32)
	set(BOOST_INCLUDE_LIBRARIES algorithm asio program_options)
	set(BOOST_ENABLE_CMAKE ON)
	FetchContent_Declare(
		Boost
		GIT_REPOSITORY https://github.com/boostorg/boost.git
		GIT_TAG boost-1.87.0
	)
	FetchContent_MakeAvailable(Boost)

	if(NOT BUILD_SHARED_LIBS)
		#statically link
		#target_link_libraries(${PROJECT_NAME} PRIVATE gdal_lib)
		target_compile_options(${PROJECT_NAME} PRIVATE /MT)
		target_link_options(${PROJECT_NAME} PRIVATE /INCREMENTAL:NO /NODEFAULTLIB:MSVCRT)
		
		#the following are static libraries required for GDAL
		#the fast majority we don't actually use but are required for compiling GDAL
		target_link_libraries(${PROJECT_NAME} PRIVATE
			aec-static 		
			charset 		
			freexl 			
			gdal 			
			geos			
			geos_c			
			geotiff			
			gif				
			iconv			
			jpeg			
			json-c			
			kmlbase			
			kmlconvenience	
			kmldom
			kmlengine
			kmlregionator
			kmlxsd
			Lerc
			libcrypto
			libcurl
			libecpg
			libecpg_compat
			libexpatMT
			libhdf5
			libhdf5_cpp
			libhdf5_hl
			libhdf5_hl_cpp
			libpgcommon
			libpgport
			libpgtypes
			libpng16
			libpq
			libsharpyuv
			libssl
			libwebp
			libwebpdecoder
			libwebpdemux
			libxml2
			lz4
			lzma
			minizip
			netcdf
			openjp2
			pcre2-8
			pcre2-16
			pcre2-32
			pcre2-posix
			proj
			qhullcpp
			qhullstatic
			qhullstatic_r
			spatialite
			sqlite3
			szip-static
			tiff
			turbojpeg
			uriparser
			zlib
			zstd
		)
		
		#the following are static libraries which are not part of GDAL but are required for GDAL on windows
		target_link_libraries(${PROJECT_NAME} PRIVATE
			crypt32
			wbemuuid
			Secur32
			ws2_32
			Normaliz
			Wldap32
		)
		
		#link boost dependencies
		target_link_libraries(${PROJECT_NAME} PRIVATE
			Boost::algorithm
			Boost::asio
			Boost::program_options
		)
	else()
	#NOP for now
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
	add_library(3pg2Test ${SOURCE_FILES})
	target_link_libraries(3pg2Test PRIVATE GDAL::GDAL)
	if(NOT MSVC)
		target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system Boost::thread Threads::Threads)
	endif()
endif()
