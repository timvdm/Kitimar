#pragma once

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    struct ValueEdge
    {
        int index = -1;
        int source = -1;
        int target = -1;

        constexpr auto operator<=>(const ValueEdge&) const noexcept = default;
    };

    namespace impl {

        consteval auto makeValueEdgeList(auto smarts, auto bonds) noexcept
        {
            if constexpr (ctll::empty(bonds))
                return std::array<ValueEdge, smarts.numBonds>{};
            else {
                auto [bond, tail] = ctll::pop_and_get_front(bonds);
                auto edges = makeValueEdgeList(smarts, tail);
                edges[bond.index] = ValueEdge{bond.index, bond.source, bond.target};
                return edges;
            }
        }

    } // namespace impl

    template<typename SmartsT>
    struct ValueEdgeList
    {
        static constexpr inline auto data = impl::makeValueEdgeList(SmartsT{}, SmartsT::bonds);

        consteval ValueEdgeList() noexcept {}
        consteval ValueEdgeList(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
