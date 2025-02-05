#ensure we have gdal
find_package(GDAL CONFIG REQUIRED)

#fetch (current) latest release from googltest github using FetchContent
include(FetchContent)
FetchContent_Declare(
	googletest
	URL https://github.com/google/googletest/archive/b62593aceb9bbd4b423804855d3eba5386eb6dc8.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Tester tests
add_executable(hello_test hello_test.cpp)
target_link_libraries(hello_test PRIVATE GTest::gtest_main)
add_test(NAME hello_test COMMAND hello_test)

# DataInput tests
add_executable(DataInputTests DataInputTests.cpp)
target_link_libraries(DataInputTests PRIVATE GTest::gtest_main GDAL::GDAL 3pg2Test)
add_test(NAME DataInputTests COMMAND DataInputTests WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

# DataOutput tests
add_executable(DataOutputTests DataOutputTests.cpp)
target_link_libraries(DataOutputTests PRIVATE GTest::gtest_main GDAL::GDAL 3pg2Test)
add_test(NAME DataOutputTests COMMAND DataOutputTests WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

include(GoogleTest)
