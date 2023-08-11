find_package(OpenBabel3 REQUIRED)

# Headers

set(KitimarOpenBabel_HDRS
    OpenBabel.hpp
)

# Library

add_library(KitimarOpenBabel INTERFACE ${KitimarOpenBabel_HDRS})
add_library(Kitimar::OpenBabel ALIAS KitimarOpenBabel)
set_target_properties(KitimarOpenBabel PROPERTIES EXPORT_NAME "OpenBabel")

target_compile_definitions(KitimarOpenBabel INTERFACE -DKITIMAR_WITH_OPENBABEL)

target_include_directories(KitimarOpenBabel INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(KitimarOpenBabel INTERFACE Kitimar::Molecule openbabel)

# Install

install(TARGETS KitimarOpenBabel
    EXPORT KitimarOpenBabelTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${KitimarOpenBabel_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/OpenBabel)

install(EXPORT KitimarOpenBabelTargets
    FILE KitimarOpenBabelTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
