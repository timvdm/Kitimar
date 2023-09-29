#pragma once

#include <array>
#include <concepts>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<int SourceIndex, int AdjIndex>
        consteval void valueDfsSearch(auto smarts, auto incidentList, auto &visitor,
                                      auto &visitedVertices, auto &visitedEdges) noexcept
        {
            constexpr auto sourceDegree = incidentList.degrees.data[SourceIndex];
            if constexpr (AdjIndex < sourceDegree) {
                constexpr auto edgeIndex = incidentList.get(SourceIndex, AdjIndex);

                if (visitedEdges[edgeIndex]) {
                    valueDfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, visitedVertices, visitedEdges);
                    return;
                }

                constexpr auto edge = incidentList.edges.data[edgeIndex];
                constexpr auto targetIndex = edge.source == SourceIndex ? edge.target : edge.source;
                auto isNewComponent = !visitedVertices[SourceIndex];
                auto isClosure = visitedVertices[targetIndex];

                visitor.visit(edgeIndex, SourceIndex, targetIndex, isNewComponent, isClosure);

                visitedEdges[edgeIndex] = true;
                visitedVertices[SourceIndex] = true;
                visitedVertices[targetIndex] = true;

                // dfs
                if (!isClosure)
                    valueDfsSearch<targetIndex, 0>(smarts, incidentList, visitor, visitedVertices, visitedEdges);

                visitor.backtrack(edgeIndex, targetIndex, isClosure);

                // next incident bond
                valueDfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, visitedVertices, visitedEdges);
            } else if constexpr (SourceIndex == 0) {
                visitor.backtrack(SourceIndex);
            }
        }

    } // namespace impl

    struct ValueDFSVisitorBase
    {
        constexpr void visit([[maybe_unused]] int edge, [[maybe_unused]] int source, [[maybe_unused]] int target,
            [[maybe_unused]] bool isNewComponent, [[maybe_unused]] bool isClosure) noexcept {}
        constexpr void backtrack([[maybe_unused]]int edge, [[maybe_unused]] int target, [[maybe_unused]] bool isClosure) noexcept {}
        constexpr void backtrack([[maybe_unused]] int source) noexcept {}
    };

    template<int SourceIndex = 0>
    consteval void valueDfsSearch(auto smarts, auto incidentList, auto &visitor) noexcept
    {
        std::array<bool, smarts.numAtoms> visitedAtoms = {};
        std::array<bool, smarts.numBonds> visitedBonds = {};
        impl::valueDfsSearch<SourceIndex, 0>(smarts, incidentList, visitor, visitedAtoms, visitedBonds);
    }

} // namespace Kitimar::CTSmarts
