function(CTSmarts_add_test name)
    add_executable(TestCTSmarts${name} TestCTSmarts${name}.cpp)
    target_link_libraries(TestCTSmarts${name} Kitimar::CTSmarts Test)
    gtest_add_tests(TARGET TestCTSmarts${name})
endfunction()


CTSmarts_add_test(Parser)
CTSmarts_add_test(Match)
CTSmarts_add_test(Count)


add_executable(TestCTSmartsOld TestCTSmarts.cpp)
target_link_libraries(TestCTSmartsOld Kitimar::CTSmarts Test)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestCTSmartsOld Kitimar::OpenBabel)
endif()
gtest_add_tests(TARGET TestCTSmartsOld)
