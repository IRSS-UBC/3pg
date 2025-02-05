mkdir build
cmake -D CMAKE_BUILD_TYPE=Debug -B build
cmake --build build -j 4
cd build/tests
ctest -C Debug
