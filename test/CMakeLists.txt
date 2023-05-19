# TestData

add_library(TestData INTERFACE TestData.hpp)
target_compile_definitions(TestData INTERFACE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_include_directories(TestData INTERFACE ${CMAKE_SOURCE_DIR}/test)

add_executable(UpdateTestData UpdateTestData.cpp)
target_link_libraries(UpdateTestData PRIVATE TestData CTLayout KitimarOpenBabel KitimarRDKit)

# Test

add_library(Test INTERFACE Test.hpp)
target_link_libraries(Test INTERFACE TestData GTest::gtest_main)

# Benchmark

add_library(Benchmark INTERFACE Benchmark.hpp)
target_link_libraries(Benchmark INTERFACE TestData Util)

# Molecule

add_executable(TestMolecule TestMolecule.cpp)
target_link_libraries(TestMolecule Molecule Test)
gtest_add_tests(TARGET TestMolecule)

# CTLayout

add_executable(TestCTLayout TestCTLayout.cpp)
target_link_libraries(TestCTLayout CTLayout Test)
gtest_add_tests(TARGET TestCTLayout)

add_executable(TestSerialize TestSerialize.cpp)
target_link_libraries(TestSerialize CTLayout CTSmarts KitimarOpenBabel Test)
gtest_add_tests(TARGET TestSerialize)

# CTSmarts

add_executable(TestCTSmarts TestCTSmarts.cpp)
target_link_libraries(TestCTSmarts CTSmarts KitimarOpenBabel Test)
gtest_add_tests(TARGET TestCTSmarts)

# OpenBabel

add_executable(TestOpenBabel TestOpenBabel.cpp)
target_link_libraries(TestOpenBabel KitimarOpenBabel CTLayout Test)
gtest_add_tests(TARGET TestOpenBabel)

# RDKit

add_executable(TestRDKit TestRDKit.cpp)
target_link_libraries(TestRDKit KitimarRDKit CTLayout Test)
gtest_add_tests(TARGET TestRDKit)

# Benchmark
add_executable(BenchmarkSerialize BenchmarkSerialize.cpp)
#target_compile_definitions(SerializeBenchmark PRIVATE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(BenchmarkSerialize CTLayout CTSmarts Util KitimarOpenBabel KitimarRDKit Benchmark)
