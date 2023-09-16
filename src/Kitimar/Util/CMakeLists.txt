# Headers

set(KitimarUtil_HDRS
    Util.hpp
    Timer.hpp
)

# Library

add_library(KitimarUtil INTERFACE ${KitimarUtil_HDRS})
add_library(Kitimar::Util ALIAS KitimarUtil)
set_target_properties(KitimarUtil PROPERTIES EXPORT_NAME "Util")

target_include_directories(KitimarUtil INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(KitimarUtil INTERFACE fmt::fmt)

# Install

install(TARGETS KitimarUtil
    EXPORT KitimarUtilTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${KitimarUtil_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/Util)

install(EXPORT KitimarUtilTargets
    FILE KitimarUtilTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
