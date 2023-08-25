include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)

add_custom_target(RunTests ctest)

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
target_link_libraries(Test INTERFACE TestData Catch2::Catch2WithMain)



function(Kitimar_add_test name)
    add_executable(Test${name} Test${name}.cpp)
    target_link_libraries(Test${name} Test)
    catch_discover_tests(Test${name}
      TEST_PREFIX "Test${name}."
      REPORTER XML
      OUTPUT_DIR .
      OUTPUT_PREFIX "Test${name}."
      OUTPUT_SUFFIX .xml
    )
endfunction()



# Benchmark
add_library(Benchmark INTERFACE Benchmark.hpp)
target_link_libraries(Benchmark INTERFACE TestData Kitimar::Util)

# Molecule
Kitimar_add_test(Molecule)

# CTLayout
Kitimar_add_test(CTLayout)

# Serialize
Kitimar_add_test(Serialize)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestSerialize Kitimar::OpenBabel)
endif()

# CTSmarts
add_subdirectory(CTSmarts)

# OpenBabel
if(Kitimar_WITH_OPENBABEL)
    Kitimar_add_test(OpenBabel)
    target_link_libraries(TestOpenBabel Kitimar::OpenBabel Kitimar::CTLayout)    
endif()

# RDKit
if(Kitimar_WITH_RDKIT)
    Kitimar_add_test(RDKit)
    target_link_libraries(TestRDKit Kitimar::RDKit Kitimar::CTLayout)
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

if(Kitimar_WITH_OPENBABEL)
    add_executable(BenchmarkCTSmartsOld BenchmarkCTSmarts.cpp)
    target_link_libraries(BenchmarkCTSmartsOld Kitimar::CTSmarts Kitimar::OpenBabel Catch2::Catch2WithMain TestData)
endif()

