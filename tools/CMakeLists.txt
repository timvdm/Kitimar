if(KITIMAR_WITH_OPENBABEL)
    add_executable(SmiToMock SmiToMock.cpp)
    target_link_libraries(SmiToMock Kitimar::Util Kitimar::OpenBabel)

    add_executable(ExtractPrimitiveFrequency ExtractPrimitiveFrequency.cpp)
    target_link_libraries(ExtractPrimitiveFrequency Kitimar::CTSmarts Kitimar::OpenBabel)
endif()
