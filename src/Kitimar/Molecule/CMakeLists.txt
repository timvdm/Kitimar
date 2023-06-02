# Headers

set(Molecule_HDRS
    Molecule.hpp
    MockMolecule.hpp
)

# Library

add_library(Molecule INTERFACE ${Molecule_HDRS})

target_include_directories(Molecule INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

# Install

install(TARGETS Molecule
    EXPORT MoleculeTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${Molecule_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/Molecule)

install(EXPORT MoleculeTargets
    FILE MoleculeTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
