ECHO y | rmdir release_build /s 2>null
mkdir release_build
cmake -DBUILD_SHARED_LIBS=OFF -DBUILD_GOOGLE_TESTS=OFF -B release_build
cmake --build release_build --target 3pg --config Release