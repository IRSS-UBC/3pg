ECHO y | rmdir build /s 2>null
mkdir build
cmake -DBUILD_SHARED_LIBS=OFF -DBUILD_GOOGLE_TESTS=OFF -B build