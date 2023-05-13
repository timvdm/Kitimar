# Library
add_library(CTSmarts INTERFACE
    CTSmarts.hpp
    Isomorphism.hpp
    SmartsActions.hpp
    SmartsAST.hpp
    SmartsGrammar.hpp
    Smarts.hpp
    SmartsMatch.hpp
    Util.hpp
)
target_include_directories(CTSmarts INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(CTSmarts INTERFACE Molecule ctre)

# Test
add_executable(TestCTSmarts Test.cpp)
target_compile_definitions(TestCTSmarts PRIVATE -DKITIMAR_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(TestCTSmarts CTSmarts KitimarOpenBabel GTest::gtest_main)
gtest_add_tests(TARGET TestCTSmarts)
