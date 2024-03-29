set(Boost_USE_MULTITHREADED OFF)
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system serialization)
find_package(RDKit REQUIRED)

# Headers

set(KitimarRDKit_HDRS
    RDKit.hpp
)

add_library(KitimarRDKit INTERFACE ${KitimarRDKit_HDRS})
add_library(Kitimar::RDKit ALIAS KitimarRDKit)
set_target_properties(KitimarRDKit PROPERTIES EXPORT_NAME "RDKit")

target_include_directories(KitimarRDKit INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
    ${RDKIT_INCLUDE_DIR}
)

target_compile_definitions(KitimarRDKit INTERFACE -DKITIMAR_WITH_RDKIT)
target_include_directories(KitimarRDKit INTERFACE ${RDKit_INCLUDE_DIRS})
target_link_libraries(KitimarRDKit INTERFACE Kitimar::Molecule
    RDKit::GraphMol_static RDKit::SmilesParse_static RDKit::FileParsers_static)


# Install

install(TARGETS KitimarRDKit
    EXPORT KitimarRDKitTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${KitimarRDKit_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/RDKit)

install(EXPORT KitimarRDKitTargets
    FILE KitimarRDKitTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
