# Headers

set(Util_HDRS
    Util.hpp
    Timer.hpp
)

# Library

add_library(Util INTERFACE ${Util_HDRS})

target_include_directories(Util INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(Util INTERFACE fmt::fmt)

# Install

install(TARGETS Util
    EXPORT UtilTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${Util_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/Util)

install(EXPORT UtilTargets
    FILE UtilTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
