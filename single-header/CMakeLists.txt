add_custom_target(CTSmartsSingleHeader
    COMMAND quom "${CMAKE_SOURCE_DIR}/src/Kitimar/CTSmarts/CTSmarts.hpp" "${CMAKE_CURRENT_SOURCE_DIR}/CTSmartsSingle.hpp"
)
