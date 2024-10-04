# Overview
3-PG with 3-PGSpatial, now parallelized and using .tif instead of .flt images. This repository is currently being maintained by the Integrated Remote Sensing Studio (IRSS), Faculty of Forestry, University of British Columbia. For documentation, see https://francoisdt14.github.io/3PG-Help/.

If you notice a bug, or have an idea for a feature request, please use GitHub issues.

# running 3-PG
download and unzip the latest release from GitHub. Then, run using the run.bat script in the release folder. 

Alternatively, you may run the executable directly. If you do this, ensure that the PROJ_DATA environment variable has been set to the 'proj' folder. Runtime warnings indicating 'proj.db' not found will occur if this is not done. This is necessary because 'proj' is a dependency of GDAL, which this new version of 3-PG utilizing .tif images relies on.

# developing 3-PG
### Prerequisistes
* [Visual Studio 2022](https://visualstudio.microsoft.com/vs/) with the C++ workload installed
* [CMake](https://cmake.org/download/) (3.15 or higher)
* [Git](https://git-scm.com/download/win)

### clone the repository
```
git clone https://github.com/IRSS-UBC/3pg.git
```

### Installing GDAL and Boost
 - instructions for installing gdal using vcpkg here: https://gdal.org/download.html#vcpkg
```
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.bat for Windows # ./bootstrap-vcpkg.sh for unix-like OS's
./vcpkg integrate install
./vcpkg install gdal
```
 - note: installing gdal will take a while, for me it took 4.4 hours.
 - add GDAL_DIR and VCPKG_ROOT to the environment variables:
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'.
   - under user variables click 'new' and set VCPKG_ROOT as the variable name, and the path to the vcpkg folder as the variable value (for me this was C:\Github\vcpkg)
   - create another new user variable with variable name GDAL_DIR and set it as the path to the folder containing the GDALConfig.cmake file (within in the vcpkg\packages folder). For me this was C:\Github\vcpkg\packages\gdal_x64-windows\share\gdal.
   - navigate the the vcpkg folder from the command line, and run the command.
```
./vcpkg install gtest
```
```
./vcpkg install boost
```
- note: installing boost will take a while, for me it took 1.7 hours.
- go to the vcpkg folder and navigate to vcpkg/installed/x64-windows/include/ (note: x64-windows may be a different folder on your installation).
- within that folder, there should be a folder called 'boost'. Copy the whole folder and paste it into the 3pg2/include directory.

### building project
 - navigate to the project directory and run either ./run_tests.bat or ./make_build.bat.
Running ./run_tests.bat will build the project with the testing folder and run the
unit tests. Running ./make_build.bat will build only the project.
 - if you see any 'missing *.dll' errors, see the 'missing dll errors' section of this readme for how to fix.

### debugging the project
 - Open 3pg.sln in the buid folder using Microsoft Visual Studio. There should be a box 'solution explorer' open. In 'solution explorer' right click on the '3pg' project and select 'set as startup project'.
 - Press f5 to build & debug.
 -  If you see any 'missing *.dll' errors, see the 'missing dll errors' section of this readme for how to fix.

### missing dll errors
 - You may see any number of 'missing *.dll' errors. To fix these, you will need to add the windows debug bin folder from vcpkg to the 'Path' environment variable:
   - Go to the vcpkg folder, click on 'installed', then 'x64-windows', then 'debug', then 'bin'. This folder should be filled with a few .py files, and lots of .dll and .pdb files. 
   - Copy the full path (for example C:\github\vcpkg\installed\x64-windows\debug\bin).
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'.
   - double click on 'Path'.
   - select 'new'.
   - paste the copied path, then click 'ok' twice.
- note: you will have to close visual studio and re-open it again for the missing dll errors, this is because it does not automatically reload environment variables whenever they are changed.
- after re-opening Visual Studio with the same project, try to build again.

### testing
 - Unit testing: there are unit tests which test the correctness of the DataInput, and DataOuptut classes essential to this version of 3PG's usage. GoogleTest is used as a testing framework. The unit tests can be ran using the run_tests.bat batch file.
 - Model testing: There are python tests which test the correctness of the model, using previously verified outputs. The python tests, since they rely on relatively large images to run, are not included in this repository.
