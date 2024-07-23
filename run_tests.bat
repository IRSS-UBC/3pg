ECHO y | rmdir build /s 2>null
mkdir build
cmake -D CMAKE_BUILD_TYPE=Debug -B build
cmake --build build
cd build/tests
ctest -C Debug
ECHO y | rmdir build /s 2>null