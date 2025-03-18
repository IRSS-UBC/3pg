# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.11)
#set compiler and options
set(CMAKE_CXX_STANDARD 20)
option(BUILD_SHARED_LIBS "build with dynamic libraries except for the standalone executable" ON)
option(BUILD_GOOGLE_TESTS "build 3pg with tests included" OFF)

#define project
project("3pg")

#dev/dynamic specific CMAKE options
if(BUILD_SHARED_LIBS)
	#getting GDAL and all it's dependencies
	set(VCPKG_ROOT ${PROJECT_SOURCE_DIR}/vcpkg)
	set(GDAL_DIR ${PROJECT_SOURCE_DIR}/vcpkg/packages/gdal_x64-windows/share/gdal)
	find_package(GDAL CONFIG REQUIRED)
	
	#add remaining include dependencies
	include_directories(${PROJECT_SOURCE_DIR}/vcpkg/installed/x64-windows/include)
	
#release/static specific CMAKE options
else()
	include_directories(${PROJECT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/include)
	link_directories(${PROJECT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/lib)
	
	#lib files required by GDAL which are within the vcpkg/instaled/x64-windows-static/lib folder
	set(GDAL_STATIC_LIB_FILES
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
endif()

#add testing specific options
if(BUILD_GOOGLE_TESTS)
	#getting gtest and all its dependencies
	set(GTEST_DIR ${PROJECT_SOURCE_DIR}/vcpkg/packages/gtest_x64-windows/share/gtest)
	find_package(GTEST CONFIG REQUIRED)
	
	#add src directory to testing include directory so we're able to include required data structures
	include_directories(${PROJECT_SOURCE_DIR}/apps)
endif()

#add subdirectories
add_subdirectory("apps")
if(BUILD_GOOGLE_TESTS)
	add_subdirectory("tests")
endif()
