# Headers

set(Molecule_HDRS
    Molecule.hpp
    Toolkit.hpp
    MockMolecule.hpp
)

# Library

add_library(KitimarMolecule INTERFACE ${Molecule_HDRS}
    MockMolecules.hpp
    NoInlineMolecule.hpp)
add_library(Kitimar::Molecule ALIAS KitimarMolecule)
set_target_properties(KitimarMolecule PROPERTIES EXPORT_NAME "Molecule")

target_include_directories(KitimarMolecule INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(KitimarMolecule INTERFACE Kitimar::Options Kitimar::Warnings)

# Install

install(TARGETS KitimarMolecule
    EXPORT KitimarMoleculeTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${Molecule_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/Molecule)

install(EXPORT KitimarMoleculeTargets
    FILE KitimarMoleculeTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
