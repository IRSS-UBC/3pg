# 3pg2
3pg but faster(?) and not a FloatGrid in sight.

# Installation

This iteration of 3pg2 is written in (mostly) modern C++ and uses the [CMake](https://cmake.org/) build system. It has been tested on Windows 10 and Ubuntu 18.04.

## Windows

### Prerequisites

* [Visual Studio 2022](https://visualstudio.microsoft.com/vs/) with the C++ workload installed
* [CMake](https://cmake.org/download/) (3.15 or higher)
* [Git](https://git-scm.com/download/win)

If PROJ_LIB environment variable is not set, set it to the directory containing proj.db. Must also be added to PATH environment variable (see 'Edit the system environment variables' in Windows search bar).

### Installing GDAL and Boost
*note: this way of installing gdal/boost is a bit janky, at least for boost. I figure it is a good idea to have working installation instructions, even if they are imperfect. This should be done after the repository is cloned, but before building with cmake.*
 - instructions for installing gdal using vcpkg here: https://gdal.org/download.html#vcpkg
```
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.bat for Windows # ./bootstrap-vcpkg.sh for unix-like OS's
./vcpkg integrate install
./vcpkg install gdal
```
 - note: installing gdal will take a while, for me it took 4.4 hours
 - add GDAL_DIR and VCPKG_ROOT to the environment variables:
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'
   - under user variables click 'new' and set VCPKG_ROOT as the variable name, and the path to the vcpkg folder as the variable value (for me this was C:\Github\vcpkg)
   - create another new user variable with variable name GDAL_DIR and set it as the path to the folder containing the GDALConfig.cmake file (within in the vcpkg\packages folder). For me this was C:\Github\vcpkg\packages\gdal_x64-windows\share\gdal
   - navigate the the vcpkg folder from the command line, and run the command
```
./vcpkg install boost
```
- note: installing boost will take a while, for me it took 1.7 hours. It may be possible to only install the required packges rather than all of boost, but this hasn't been tested yet.
- go to the vcpkg folder and navigate to vcpkg/installed/x64-windows/include/ (note: x64-windows may be a different folder on your installation)
- within that folder, there should be a folder called 'boost'. Copy the whole folder and paste it into the 3pg2/include directory.

### Building

1. Clone the repository
	```{bash}
	git clone
	cd 3pg2
	```
2. Create a build directory
	```{bash}
	mkdir build
	cd build
	```
3. Run CMake

	```{bash}
	cmake ..
	```
4. Run Make

	```{bash}
	cmake --build . --config Release
	```
