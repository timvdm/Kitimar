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

        constexpr DfsEdgeListVisitor() noexcept {}
        constexpr DfsEdgeListVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            edges[nextEdgeIndex++] = DfsEdge{{source, target}, edge, isClosure};
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept {}

        constexpr void backtrack(int source) noexcept {}


        std::array<DfsEdge, SmartsT::numBonds> edges;
        int nextEdgeIndex = 0;
    };

    template<typename SmartsT, typename AdjacencyListT>
    constexpr auto makeDfsEdgeList()
    {
        DfsEdgeListVisitor<SmartsT> visitor;
        dfsSearch(SmartsT{}, visitor, AdjacencyListT{});
        return visitor.edges;
    }


    template<typename SmartsT, typename AdjacencyListT>
    struct DfsEdgeList
    {
        static constexpr inline auto data = makeDfsEdgeList<SmartsT, AdjacencyListT>();

        static constexpr inline DfsEdge get(int index) noexcept
        {
            return data[index];
        }

        constexpr DfsEdgeList() noexcept {}
        constexpr DfsEdgeList(SmartsT, AdjacencyListT) noexcept {}
    };

}
