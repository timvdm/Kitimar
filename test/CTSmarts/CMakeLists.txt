
add_library(TestCTSmarts INTERFACE TestCTSmarts.hpp)

Kitimar_add_test(CTSmartsParser)
Kitimar_add_test(CTSmartsMatch)
Kitimar_add_test(CTSmartsCount)
Kitimar_add_test(CTSmartsMap)
Kitimar_add_test(CTSmartsMaps)
Kitimar_add_test(CTSmartsCapture)
Kitimar_add_test(CTSmartsCaptures)


Kitimar_add_test(CTSmartsOld)
if(Kitimar_WITH_OPENBABEL)
    target_link_libraries(TestCTSmartsOld Kitimar::OpenBabel)
endif()
