# Headers

set(CTSmarts_HDRS
    CTSmarts.hpp
    Smarts.hpp
    MatchExpr.hpp
    Isomorphism.hpp
    AST.hpp
    Util.hpp
    Parser/Grammar.hpp
    Parser/Actions.hpp
    Algorithm/EdgeList.hpp
    Algorithm/VertexDegree.hpp
    Algorithm/IncidentList.hpp
    Algorithm/DfsSearch.hpp
    Algorithm/DfsSearchEvents.hpp
    Algorithm/DfsEdgeList.hpp
    Algorithm/CycleMembership.hpp
    Algorithm/DfsBondList.hpp
    Algorithm/CtllSort.hpp
    Optimizer/PrimitiveFrequency.hpp
    Optimizer/ExpressionFrequency.hpp
    Optimizer/OptimizeExpression.hpp
)

# Library

add_library(CTSmarts INTERFACE ${CTSmarts_HDRS})

target_include_directories(CTSmarts INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(CTSmarts INTERFACE ctre::ctre Molecule)
if(KITIMAR_WITH_COROUTINES)
    target_link_libraries(CTSmarts INTERFACE cppcoro)
endif()

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
