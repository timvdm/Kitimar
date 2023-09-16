add_custom_target(CTSmartsSingleHeader
    COMMAND quom
    --include_guard "KITIMAR_.*_HPP"
    "${CMAKE_SOURCE_DIR}/src/Kitimar/CTSmarts/CTSmarts.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/CTSmartsSingle.hpp"
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/CTSmartsSingle.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar)
