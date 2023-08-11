# TestData

add_library(TestData INTERFACE TestData.hpp)
target_compile_definitions(TestData INTERFACE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_include_directories(TestData INTERFACE ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(TestData INTERFACE Kitimar::CTLayout)

add_executable(UpdateTestData UpdateTestData.cpp)
target_link_libraries(UpdateTestData TestData Kitimar::CTLayout)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(UpdateTestData Kitimar::OpenBabel)
endif()
if(Kitimar_WITH_RDKIT)
    target_link_libraries(UpdateTestData Kitimar::RDKit)
endif()

# Test

add_library(Test INTERFACE Test.hpp)
target_link_libraries(Test INTERFACE TestData GTest::gtest_main)

# Benchmark

add_library(Benchmark INTERFACE Benchmark.hpp)
target_link_libraries(Benchmark INTERFACE TestData Kitimar::Util)

# Molecule

add_executable(TestMolecule TestMolecule.cpp)
target_link_libraries(TestMolecule Kitimar::Molecule Test)
gtest_add_tests(TARGET TestMolecule)

# CTLayout

add_executable(TestCTLayout TestCTLayout.cpp)
target_link_libraries(TestCTLayout Kitimar::CTLayout Test fmt::fmt)
gtest_add_tests(TARGET TestCTLayout)

add_executable(TestSerialize TestSerialize.cpp)
target_link_libraries(TestSerialize Kitimar::CTLayout Test fmt::fmt)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestSerialize Kitimar::OpenBabel)
endif()
gtest_add_tests(TARGET TestSerialize)

# CTSmarts

add_executable(TestCTSmarts TestCTSmarts.cpp)
target_link_libraries(TestCTSmarts Kitimar::CTSmarts Test)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestCTSmarts Kitimar::OpenBabel)
endif()
gtest_add_tests(TARGET TestCTSmarts)

# OpenBabel

if(Kitimar_WITH_OPENBABEL)
    add_executable(TestOpenBabel TestOpenBabel.cpp)
    target_link_libraries(TestOpenBabel Kitimar::OpenBabel Kitimar::CTLayout Test)
    gtest_add_tests(TARGET TestOpenBabel)
endif()

# RDKit

if(Kitimar_WITH_RDKIT)
    add_executable(TestRDKit TestRDKit.cpp)
    target_link_libraries(TestRDKit Kitimar::RDKit Kitimar::CTLayout Test)
    gtest_add_tests(TARGET TestRDKit)
endif()

# Benchmark

add_executable(BenchmarkSerialize BenchmarkSerialize.cpp)
target_link_libraries(BenchmarkSerialize Kitimar::CTLayout Kitimar::CTSmarts Kitimar::Util Benchmark)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(BenchmarkSerialize Kitimar::OpenBabel)
endif()
if(Kitimar_WITH_RDKIT)
    target_link_libraries(BenchmarkSerialize Kitimar::RDKit)
endif()

#if(Kitimar_WITH_OPENBABEL AND Kitimar_WITH_RDKIT)
if(Kitimar_WITH_OPENBABEL)
    add_executable(BenchmarkCTSmarts BenchmarkCTSmarts.cpp)
    #target_link_libraries(BenchmarkCTSmarts CTSmarts KitimarOpenBabel KitimarRDKit benchmark::benchmark)
    target_link_libraries(BenchmarkCTSmarts Kitimar::CTSmarts Kitimar::OpenBabel benchmark::benchmark)
endif()

