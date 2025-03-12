# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
cmake_minimum_required (VERSION 3.11)
set(CMAKE_CXX_STANDARD 20)

option(BUILD_SHARED_LIBS "build with dynamic libraries except for the standalone executable" ON)
option(BUILD_GOOGLE_TESTS "build 3pg with tests included" OFF)

# Enable use of _ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project("3pg")

#if we're building a release version, add the speediest compiler optimizations
#if(NOT BUILD_SHARED_LIBS AND MSVC)
#	add_compile_options("-O2")
#endif()

if(BUILD_SHARED_LIBS)
	#for building dynamically, it's easiest to use CMAKE find_package
	
	#getting GDAL and all it's dependencies
	set(VCPKG_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg)
	set(GDAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/gdal_x64-windows/share/gdal)
	find_package(GDAL CONFIG REQUIRED)
	
	#getting remaining dependencies (boost)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows/include)
else()
	#for building statically, it's easiest to add the libs directly
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/include)
	link_directories(${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows-static/lib)
	
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add subdirectories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_subdirectory("apps")
if(BUILD_GOOGLE_TESTS)
	add_subdirectory("tests")
endif()
