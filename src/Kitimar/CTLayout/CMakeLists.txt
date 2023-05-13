# Library
add_library(CTLayout INTERFACE CTLayout.hpp)
target_include_directories(CTLayout INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(CTLayout INTERFACE ctre)

# Test
add_executable(TestCTLayout Test.cpp)
target_link_libraries(TestCTLayout CTLayout GTest::gtest_main)
gtest_add_tests(TARGET TestCTLayout)
