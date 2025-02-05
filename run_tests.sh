#!/bin/bash

mkdir build
cmake -G "Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug -B build
cmake --build build -j 4
cd build/tests
ctest -C Debug
