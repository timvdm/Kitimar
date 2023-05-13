# Library
add_library(Serialize INTERFACE Serialize.hpp)
target_include_directories(Serialize INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(Serialize INTERFACE Molecule CTLayout mio)

# Test
add_executable(TestSerialize Test.cpp)
target_compile_definitions(TestSerialize PRIVATE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(TestSerialize KitimarOpenBabel Serialize CTSmarts GTest::gtest_main)
gtest_add_tests(TARGET TestSerialize)

# Benchmark
add_executable(SerializeBenchmark Benchmark.cpp)
target_compile_definitions(SerializeBenchmark PRIVATE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(SerializeBenchmark KitimarOpenBabel Serialize Util CTSmarts)
