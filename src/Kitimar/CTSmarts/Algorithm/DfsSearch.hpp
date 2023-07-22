#pragma once

#include <array>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename Visitor, typename IncidentListT, int SourceIndex, int AdjIndex>
    constexpr void dfsSearchHelper(std::array<bool, SmartsT::numAtoms> &visitedVertices,
                                   std::array<bool, SmartsT::numBonds> &visitedEdges,
                                   Visitor &visitor)
    {
        constexpr auto sourceDegree = IncidentListT::degrees.data[SourceIndex];
        if constexpr (AdjIndex < sourceDegree) {
            constexpr auto edgeIndex = IncidentListT::get(SourceIndex, AdjIndex);

            if (visitedEdges[edgeIndex]) {
                dfsSearchHelper<SmartsT, Visitor, IncidentListT, SourceIndex, AdjIndex + 1>(visitedVertices, visitedEdges, visitor);
                return;
            }

            constexpr auto edge = IncidentListT::edges.data[edgeIndex];
            constexpr auto targetIndex = edge.source == SourceIndex ? edge.target : edge.source;
            auto isNewComponent = !visitedVertices[SourceIndex];
            auto isClosure = visitedVertices[targetIndex];

            visitor.visit(edgeIndex, SourceIndex, targetIndex, isNewComponent, isClosure);

            visitedEdges[edgeIndex] = true;
            visitedVertices[SourceIndex] = true;
            visitedVertices[targetIndex] = true;

            // dfs
            if (!isClosure)
                dfsSearchHelper<SmartsT, Visitor, IncidentListT, targetIndex, 0>(visitedVertices, visitedEdges, visitor);

            visitor.backtrack(edgeIndex, targetIndex, isClosure);

            // next incident bond
            dfsSearchHelper<SmartsT, Visitor, IncidentListT, SourceIndex, AdjIndex + 1>(visitedVertices, visitedEdges, visitor);
        } else if constexpr (SourceIndex == 0) {
            visitor.backtrack(SourceIndex);
        }
    }

    template<typename SmartsT, typename Visitor, typename IncidentListT>
    constexpr void dfsSearch(SmartsT, Visitor &visitor, IncidentListT)
    {
        std::array<bool, SmartsT::numAtoms> visitedAtoms = {};
        std::array<bool, SmartsT::numBonds> visitedBonds = {};
        dfsSearchHelper<SmartsT, Visitor, IncidentListT, 0, 0>(visitedAtoms, visitedBonds, visitor);
    }

}
