# Overview
3-PG with 3-PGSpatial, now parallelized and using .tif instead of .flt images. This repository is currently being maintained by the Integrated Remote Sensing Studio (IRSS), Faculty of Forestry, University of British Columbia. For documentation, see https://francoisdt14.github.io/3PG-Help/.

If you notice a bug, or have an idea for a feature request, please use GitHub issues.

# Running 3-PG
Download and unzip the latest release from GitHub. Then, run using the run.bat script in the release folder. In Linux, you need to source the run.sh script to allow to temporarily export PROJ, if available. 

Alternatively, you may run the executable directly. If you do this, ensure that the PROJ_DATA environment variable has been set to the 'proj' folder. Runtime warnings indicating 'proj.db' not found will occur if this is not done. This is necessary because 'proj' is a dependency of GDAL, which this new version of 3-PG utilizing .tif images relies on.

To compile 3-PG from source in release mode, first check the prerequisistes, then run:
```
# Windows
cmake -D CMAKE_BUILD_TYPE=Release -B build
cmake --build build -j 4

# Ubuntu
cmake -G "Unix Makefiles" -D CMAKE_BUILD_TYPE=Release -B build
cmake --build build -j 4
```

# Developing 3-PG
### Prerequisistes
**On Windows**:
* [Visual Studio 2022](https://visualstudio.microsoft.com/vs/) with the C++ workload installed
* [CMake](https://cmake.org/download/) (3.15 or higher)
* [Git](https://git-scm.com/download/win)

**On Linux**:
* [g++](https://gcc.gnu.org/) (14 or higher)
* [CMake](https://cmake.org/download/) (3.15 or higher)
* [Git](https://git-scm.com/downloads/linux)

### Clone the repository
```
git clone https://github.com/IRSS-UBC/3pg.git
```

### Installing dependencies
**On Windows**:
 - instructions for installing gdal using vcpkg here: https://gdal.org/download.html#vcpkg
```
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install
.\vcpkg install gdal
```
 - note: installing gdal will take a while, for me it took 4.4 hours.
 - add GDAL_DIR and VCPKG_ROOT to the environment variables:
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'.
   - under user variables click 'new' and set VCPKG_ROOT as the variable name, and the path to the vcpkg folder as the variable value (for me this was C:\Github\vcpkg)
   - create another new user variable with variable name GDAL_DIR and set it as the path to the folder containing the GDALConfig.cmake file (within in the vcpkg\packages folder). For me this was C:\Github\vcpkg\packages\gdal_x64-windows\share\gdal.
   - navigate to the vcpkg folder from the command line, and run the command.
```
.\vcpkg install gtest
```

**On Linux**:
- gcc-14 is only available on Ubuntu 24.04 and higher from official repositories, so if you're using an older distribution you can compile it from source
```
sudo apt install build-essential
sudo apt install libmpfr-dev libgmp3-dev libmpc-dev -y
wget http://ftp.gnu.org/gnu/gcc/gcc-14.2.0/gcc-14.2.0.tar.gz
tar -xf gcc-14.2.0.tar.gz
cd gcc-14.2.0
./configure -v --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu --prefix=/usr/local/gcc-14.2.0 --enable-checking=release --enable-languages=c,c++ --disable-multilib --program-suffix=-14.2.0
make
sudo make install
```
And you can make it the default compiler
```
sudo update-alternatives --install /usr/bin/g++ g++ /usr/local/gcc-14.2.0/bin/g++-14.2.0 14
sudo update-alternatives --install /usr/bin/gcc gcc /usr/local/gcc-14.2.0/bin/gcc-14.2.0 14
```
- install gdal and boost via apt (note that in Linux we link gdal dynamically)
```
sudo apt install gdal-bin libgdal-dev libboost-all-dev
```

### Building project
**On Windows**:
 - navigate to the project directory and run either ./run_tests.bat or ./make_build.bat.
Running ./run_tests.bat will build the project with the testing folder and run the
unit tests. Running ./make_build.bat will build only the project.
 - if you see any 'missing *.dll' errors, see the 'missing dll errors' section of this readme for how to fix.

**On Linux**:
- running ./run_tests.sh will build the project with the testing folder and run the
unit tests. Running ./make_build.sh will build only the project.
- if you get an error like 'GLIBCXX_3.4.31/32 not found', you probably need to update libstdc++6
 ```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install --only-upgrade libstdc++6
```

### Debugging the project
 - Open 3pg.sln in the build folder using Microsoft Visual Studio. There should be a box 'solution explorer' open. In 'solution explorer' right click on the '3pg' project and select 'set as startup project'.
 - Press f5 to build & debug.
 -  If you see any 'missing *.dll' errors, see the 'missing dll errors' section of this readme for how to fix.
 - For compilations with gcc, you can debug with gdb in any IDE that supports it

### Missing dll errors
 - You may see any number of 'missing *.dll' errors. To fix these, you will need to add the windows debug bin folder from vcpkg to the 'Path' environment variable:
   - Go to the vcpkg folder, click on 'installed', then 'x64-windows', then 'debug', then 'bin'. This folder should be filled with a few .py files, and lots of .dll and .pdb files. 
   - Copy the full path (for example C:\github\vcpkg\installed\x64-windows\debug\bin).
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'.
   - double-click on 'Path'.
   - select 'new'.
   - paste the copied path, then click 'ok' twice.
- note: you will have to close visual studio and re-open it again for the missing dll errors, this is because it does not automatically reload environment variables whenever they are changed.
- after re-opening Visual Studio with the same project, try to build again.

### Testing
 - Unit testing: there are unit tests which test the correctness of the DataInput, and DataOuptut classes essential to this version of 3PG's usage. GoogleTest is used as a testing framework. The unit tests can be ran using the run_tests.bat batch file.
 - Model testing: There are python tests which test the correctness of the model, using previously verified outputs. The python tests, since they rely on relatively large images to run, are not included in this repository.
