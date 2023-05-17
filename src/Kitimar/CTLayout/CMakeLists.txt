# Library
add_library(CTLayout INTERFACE
    Types.hpp
    Layout.hpp
    Object.hpp    
    Serialize.hpp
    Molecule.hpp
    CTLayout.hpp
)
target_include_directories(CTLayout INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)
target_link_libraries(CTLayout INTERFACE ctre mio)

install(TARGETS CTLayout
    EXPORT CTLayoutTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES Types.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/CTLayout)

install(EXPORT CTLayoutTargets
    FILE CTLayoutTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CTLayout
)


# Test
add_executable(TestCTLayout Test.cpp)
target_link_libraries(TestCTLayout CTLayout GTest::gtest_main)
gtest_add_tests(TARGET TestCTLayout)
