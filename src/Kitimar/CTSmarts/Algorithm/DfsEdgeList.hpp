#pragma once

#include "EdgeList.hpp"

namespace Kitimar::CTSmarts {

    struct DfsEdge : Edge
    {
        int index = -1;
        bool closure = false;
    };

    template<typename SmartsT>
    struct DfsEdgeListVisitor
    {

        consteval DfsEdgeListVisitor() noexcept {}
        consteval DfsEdgeListVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            edges[nextEdgeIndex++] = DfsEdge{{source, target}, edge, isClosure};
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept {}

        constexpr void backtrack(int source) noexcept {}


        std::array<DfsEdge, SmartsT::numBonds> edges;
        int nextEdgeIndex = 0;
    };

    template<typename SmartsT, typename IncidentListT>
    consteval auto makeDfsEdgeList()
    {
        DfsEdgeListVisitor<SmartsT> visitor;
        dfsSearch(SmartsT{}, visitor, IncidentListT{});
        return visitor.edges;
    }

    template<typename SmartsT, typename IncidentListT>
    struct DfsEdgeList
    {
        static constexpr inline auto data = makeDfsEdgeList<SmartsT, IncidentListT>();

        static consteval inline DfsEdge get(int index) noexcept
        {
            return data[index];
        }

        consteval DfsEdgeList() noexcept {}
        consteval DfsEdgeList(SmartsT, IncidentListT) noexcept {}
    };

}
