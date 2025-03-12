mkdir test_build
cmake -D CMAKE_BUILD_TYPE=Debug -DBUILD_GOOGLE_TESTS=ON -B test_build
cmake --build test_build -j 4
cd test_build/tests
ctest -C Debug
