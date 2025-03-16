enable_testing()

# Tester tests
add_executable(hello_test hello_test.cpp)
target_link_libraries(hello_test PRIVATE GTest::gtest_main)
add_test(NAME hello_test COMMAND hello_test)

# DataInput tests
add_executable(DataInputTests DataInputTests.cpp)
target_link_libraries(DataInputTests PRIVATE GTest::gtest_main GDAL::GDAL 3pg)
add_test(NAME DataInputTests COMMAND DataInputTests WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

# DataOutput tests
add_executable(DataOutputTests DataOutputTests.cpp)
target_link_libraries(DataOutputTests PRIVATE GTest::gtest_main GDAL::GDAL 3pg)
add_test(NAME DataOutputTests COMMAND DataOutputTests WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

include(GoogleTest)
