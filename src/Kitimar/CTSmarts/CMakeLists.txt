# Headers

set(CTSmarts_HDRS
    Isomorphism.hpp
    SmartsActions.hpp
    SmartsAST.hpp
    SmartsGrammar.hpp
    Smarts.hpp
    SmartsMatch.hpp
    Util.hpp
    CTSmarts.hpp
)

# Library

add_library(CTSmarts INTERFACE ${CTSmarts_HDRS})

target_include_directories(CTSmarts INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(CTSmarts INTERFACE ctre::ctre Molecule)

# Install

install(TARGETS CTSmarts
    EXPORT CTSmartsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CTSmarts_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/CTSmarts)

install(EXPORT CTSmartsTargets
    FILE CTSmartsTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
