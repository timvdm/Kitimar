if(KITIMAR_WITH_OPENBABEL)
    add_executable(SmiToMock SmiToMock.cpp)
    target_link_libraries(SmiToMock Util Molecule KitimarOpenBabel)
endif()
