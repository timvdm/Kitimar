#pragma once

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    struct Edge
    {
        int source = -1;
        int target = -1;

        constexpr bool operator<=>(const Edge&) const = default;
    };

    template<typename SmartsT>
    consteval auto makeEdgeList(auto bonds) noexcept
    {
        if constexpr (ctll::empty(bonds))
            return std::array<Edge, SmartsT::numBonds>{};
        else {
            auto [bond, tail] = ctll::pop_and_get_front(bonds);
            auto edges = makeEdgeList<SmartsT>(tail);
            edges[bond.index] = Edge{bond.source, bond.target};
            return edges;
        }
    }

    template<typename SmartsT>
    struct EdgeList
    {
        static constexpr inline auto data = makeEdgeList<SmartsT>(SmartsT::bonds);

        consteval EdgeList() noexcept {}
        consteval EdgeList(SmartsT) noexcept {}
    };

}
