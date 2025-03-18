ECHO y | rmdir build /s 2>null
mkdir dev_build
cmake -DBUILD_SHARED_LIBS=ON -DBUILD_GOOGLE_TESTS=OFF -B dev_build
cmake --build dev_build --target 3pg --config Debug