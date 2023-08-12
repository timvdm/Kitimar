

Kitimar_add_test(CTSmartsParser)
Kitimar_add_test(CTSmartsMatch)
Kitimar_add_test(CTSmartsCount)


Kitimar_add_test(CTSmartsOld)
#add_executable(TestCTSmartsOld TestCTSmarts.cpp)
#target_link_libraries(TestCTSmartsOld Kitimar::CTSmarts Test)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestCTSmartsOld Kitimar::OpenBabel)
endif()
#gtest_add_tests(TARGET TestCTSmartsOld)
