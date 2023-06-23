# TestData

add_library(TestData INTERFACE TestData.hpp)
target_compile_definitions(TestData INTERFACE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_include_directories(TestData INTERFACE ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(TestData INTERFACE CTLayout)

add_executable(UpdateTestData UpdateTestData.cpp)
target_link_libraries(UpdateTestData TestData CTLayout)
if(KITIMAR_WITH_OPENBABEL)
    target_link_libraries(UpdateTestData KitimarOpenBabel)
endif()
if(KITIMAR_WITH_RDKIT)
    target_link_libraries(UpdateTestData KitimarRDKit)
endif()

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
target_link_libraries(TestCTLayout CTLayout Test fmt::fmt)
gtest_add_tests(TARGET TestCTLayout)

add_executable(TestSerialize TestSerialize.cpp)
target_link_libraries(TestSerialize CTLayout Test fmt::fmt)
if(KITIMAR_WITH_OPENBABEL)
    target_link_libraries(TestSerialize KitimarOpenBabel)
endif()
gtest_add_tests(TARGET TestSerialize)

# CTSmarts

add_executable(TestCTSmarts TestCTSmarts.cpp)
target_link_libraries(TestCTSmarts CTSmarts Test)
if(KITIMAR_WITH_OPENBABEL)
    target_link_libraries(TestCTSmarts KitimarOpenBabel)
endif()
gtest_add_tests(TARGET TestCTSmarts)

# OpenBabel

if(KITIMAR_WITH_OPENBABEL)
    add_executable(TestOpenBabel TestOpenBabel.cpp)
    target_link_libraries(TestOpenBabel KitimarOpenBabel CTLayout Test)
    gtest_add_tests(TARGET TestOpenBabel)
endif()

# RDKit

if(KITIMAR_WITH_RDKIT)
    add_executable(TestRDKit TestRDKit.cpp)
    target_link_libraries(TestRDKit KitimarRDKit CTLayout Test)
    gtest_add_tests(TARGET TestRDKit)
endif()

# Benchmark

add_executable(BenchmarkSerialize BenchmarkSerialize.cpp)
#target_compile_definitions(SerializeBenchmark PRIVATE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(BenchmarkSerialize CTLayout CTSmarts Util Benchmark)
if(KITIMAR_WITH_OPENBABEL)
    target_link_libraries(BenchmarkSerialize KitimarOpenBabel)
endif()
if(KITIMAR_WITH_RDKIT)
    target_link_libraries(BenchmarkSerialize KitimarRDKit)
endif()

if(KITIMAR_WITH_OPENBABEL)
    add_executable(BenchmarkCTSmarts BenchmarkCTSmarts.cpp)
    target_link_libraries(BenchmarkCTSmarts CTSmarts KitimarOpenBabel benchmark::benchmark)
endif()

