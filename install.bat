git submodule update --init
cd vcpkg
call bootstrap-vcpkg.bat
call vcpkg integrate install

::install statically built libraries
call vcpkg install gdal:x64-windows-static boost-algorithm:x64-windows-static boost-asio:x64-windows-static boost-program-options:x64-windows-static

::install packages for dynamically linked version
call vcpkg install gdal boost-algorithm boost-asio boost-program-options gtest --recurse 