# Headers

set(CTLayout_HDRS
    Source.hpp
    Sink.hpp
    Value.hpp
    Struct.hpp
    Vector.hpp
    Molecule.hpp
    CTLayout.hpp

    #Types.hpp
    #Layout.hpp
    #Object.hpp
    #Serialize.hpp


)

# Library

add_library(KitimarCTLayout INTERFACE ${CTLayout_HDRS})
add_library(Kitimar::CTLayout ALIAS KitimarCTLayout)
set_target_properties(KitimarCTLayout PROPERTIES EXPORT_NAME "CTLayout")

target_include_directories(KitimarCTLayout INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(KitimarCTLayout INTERFACE ctre::ctre mio::mio Kitimar::Util Kitimar::Molecule)

# Install

install(TARGETS KitimarCTLayout
    EXPORT KitimarCTLayoutTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CTLayout_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/CTLayout)

install(EXPORT KitimarCTLayoutTargets
    FILE KitimarCTLayoutTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
