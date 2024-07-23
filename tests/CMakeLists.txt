find_package(GTest)
enable_testing()
set(TST_NAME ${CMAKE_PROJECT_NAME}_tst)
set(TST_FILES
		DataInputTest.cpp
		DataOutputTest.cpp
		GDALRasterImageTest.cpp
)

add_executable(${TST_NAME} ${TST_FILES})
target_link_libraries(
  ${TST_NAME}
  ${CMAKE_PROJECT_NAME}_lib
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(${TST_NAME})