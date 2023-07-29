#pragma once

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    struct Edge
    {
        int index = -1;
        int source = -1;
        int target = -1;

        constexpr bool operator<=>(const Edge&) const noexcept = default;
    };

    namespace impl {

        consteval auto makeEdgeList(auto smarts, auto bonds) noexcept
        {
            if constexpr (ctll::empty(bonds))
                return std::array<Edge, smarts.numBonds>{};
            else {
                auto [bond, tail] = ctll::pop_and_get_front(bonds);
                auto edges = makeEdgeList(smarts, tail);
                edges[bond.index] = Edge{bond.index, bond.source, bond.target};
                return edges;
            }
        }

    } // namespace impl

    template<typename SmartsT>
    struct EdgeList
    {
        static constexpr inline auto data = impl::makeEdgeList(SmartsT{}, SmartsT::bonds);

        consteval EdgeList() noexcept {}
        consteval EdgeList(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
