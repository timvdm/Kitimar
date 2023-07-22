#pragma once

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    struct Edge
    {
        int source = -1;
        int target = -1;
    };

    template<typename SmartsT>
    consteval auto makeEdgeList(ctll::empty_list) noexcept
    {
        return std::array<Edge, SmartsT::numBonds>{};
    }

    template<typename SmartsT, typename Bond, typename ...Bonds>
    consteval auto makeEdgeList(ctll::list<Bond, Bonds...> bonds) noexcept
    {
        auto [bond, tail] = ctll::pop_and_get_front(bonds);
        auto edges = makeEdgeList<SmartsT>(tail);
        auto bondIndex = SmartsT::numBonds - ctll::size(bonds);
        edges[bondIndex] = Edge{bond.source, bond.target};
        return edges;
    }

    template<typename SmartsT>
    struct EdgeList
    {
        static constexpr inline auto data = makeEdgeList<SmartsT>(SmartsT::bonds);

        consteval EdgeList() noexcept {}
        consteval EdgeList(SmartsT) noexcept {}
    };

}
