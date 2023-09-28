# Headers

set(CTSmarts_HDRS
    CTSmarts.hpp
    Smarts.hpp
    MatchExpr.hpp
    Mapping.hpp
    Config.hpp
    Isomorphism.hpp
    AST/AST.hpp
    AST/Primitives.hpp
    AST/Operations.hpp
    AST/BasicSmarts.hpp
    AST/ContainsExpr.hpp
    AST/ReplaceExpr.hpp
    Util/Util.hpp
    Util/Ctll.hpp
    Util/CtllSort.hpp
    Parser/Grammar.hpp
    Parser/Actions.hpp
    Graph/DfsBond.hpp
    Graph/ValueEdgeList.hpp
    Graph/TypeEdgeList.hpp
    Graph/ValueVertexDegree.hpp
    Graph/TypeVertexDegree.hpp
    Graph/ValueIncidentList.hpp
    Graph/TypeIncidentList.hpp
    Graph/DfsSearchEventType.hpp
    Graph/ValueDfsSearch.hpp
    Graph/TypeDfsSearch.hpp
    Graph/ValueDfsSearchEvents.hpp
    Graph/TypeDfsSearchEvents.hpp
    Graph/ValueDfsEdgeList.hpp
    Graph/TypeDfsEdgeList.hpp
    Graph/ValueCycleMembership.hpp
    Graph/TypeCycleMembership.hpp
    Graph/ValueDfsBondList.hpp
    Graph/TypeDfsBondList.hpp
    Filter/UnconditionalFilter.hpp
    Filter/NumAtomBondFilter.hpp
    Filter/ElementFilter.hpp
    Filter/FilterPolicy.hpp
    Optimizer/PrimitiveFrequency.hpp
    Optimizer/ExpressionFrequency.hpp
    Optimizer/OptimizeExpression.hpp
    Optimizer/AtomFrequency.hpp
    Optimizer/OptimizeIncidentList.hpp
    Optimizer/Optimizer.hpp
    API/Util.hpp
    API/Match.hpp
    API/Count.hpp
    API/Map.hpp
    API/Maps.hpp
    API/Capture.hpp
    API/Captures.hpp
    API/Find.hpp
)

# Library

add_library(KitimarCTSmarts INTERFACE ${CTSmarts_HDRS})
add_library(Kitimar::CTSmarts ALIAS KitimarCTSmarts)
set_target_properties(KitimarCTSmarts PROPERTIES EXPORT_NAME "CTSmarts")

target_include_directories(KitimarCTSmarts INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>"
)

target_link_libraries(KitimarCTSmarts INTERFACE ctre::ctre Kitimar::Molecule)

target_precompile_headers(KitimarCTSmarts INTERFACE ${CTSmarts_HDRS})

# Install

install(TARGETS KitimarCTSmarts
    EXPORT KitimarCTSmartsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CTSmarts_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Kitimar/CTSmarts)

install(EXPORT KitimarCTSmartsTargets
    FILE KitimarCTSmartsTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)
