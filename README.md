### 3pg2
TEST
3pg but faster(?) and not a FloatGrid in sight.


<br></br>
### Prerequisistes
* [Visual Studio 2022](https://visualstudio.microsoft.com/vs/) with the C++ workload installed
* [CMake](https://cmake.org/download/) (3.15 or higher)
* [Git](https://git-scm.com/download/win)

<br></br>
### clone the repository
```
git clone https://github.com/IRSS-UBC/3pg2.git
```

<br></br>
### Installing GDAL and Boost
*note: this way of installing gdal/boost is a bit janky, at least for boost. I figure it is a good idea to have working installation instructions, even if they are imperfect.*
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
- note: installing boost will take a while, for me it took 1.7 hours. It may be possible to only install the required packges rather than all of boost, but this hasn't been tested yet.
- go to the vcpkg folder and navigate to vcpkg/installed/x64-windows/include/ (note: x64-windows may be a different folder on your installation).
- within that folder, there should be a folder called 'boost'. Copy the whole folder and paste it into the 3pg2/include directory.

<br></br>
### building and debugging project
 - navigate to the project directory and run the following command:
```
mkdir build && cd build && cmake ..
```
 - either delete the build folder, or don't run mkdir build if the build folder already exists.
 - navigate to the build folder in file explorer, and open the '3pg.vcxproj' file.
 - Microsoft Visual Studio should open, and there should be a box 'solution explorer' open. In 'solution explorer' right click on the '3pg' project and select 'set as startup project'.
 - press f5 to build & debug.
 - at this point you may see any number of 'missing *.dll' errors. To fix these, you will need to add the windows debug bin folder from vcpkg to the 'Path' environment variable:
   - Go to the vcpkg folder, click on 'installed', then 'x64-windows' (this may be different for you depending on os), then 'debug', then 'bin'. This folder should be filled with a few .py files, and lots of .dll and .pdb files. 
   - Copy the full path (for me this is C:\Github\vcpkg\installed\x64-windows\debug\bin).
   - type 'env' into the Windows search bar and select 'Edit environment variables for your account'.
   - double click on 'Path'.
   - select 'new'.
   - paste the copied path, then click 'ok' twice.
- note: you will have to close visual studio and re-open it again for the missing dll errors, this is because it does not automatically reload environment variables whenever they are changed.
- after re-opening visual studio with the same project, try pressing f5 again to build and debug.
