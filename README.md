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