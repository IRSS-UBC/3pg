mkdir build 2>null
cmake -D CMAKE_BUILD_TYPE=Debug -B build
cmake --build build
cd build/tests
ctest -C Debug