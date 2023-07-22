if(KITIMAR_WITH_OPENBABEL)
    add_executable(SmiToMock SmiToMock.cpp)
    target_link_libraries(SmiToMock Util Molecule KitimarOpenBabel)

    add_executable(ExtractPrimitiveFrequency ExtractPrimitiveFrequency.cpp)
    target_link_libraries(ExtractPrimitiveFrequency CTSmarts Molecule KitimarOpenBabel)
endif()
