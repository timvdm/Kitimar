#pragma once

#include "DfsBond.hpp"

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<int DfsEdgeIndex>
        consteval auto makeValueDfsBondListHelper(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            if constexpr (DfsEdgeIndex == smarts.numBonds) {
                return ctll::empty_list{};
            } else {
                constexpr auto edge = dfsEdgeList.data[DfsEdgeIndex];
                constexpr auto expr = get<edge.index>(smarts.bonds).expr;
                auto source = get<edge.source>(smarts.atoms);
                auto target = get<edge.target>(smarts.atoms);
                constexpr auto isCyclic = cycleMembership.edges[edge.index];

                auto dfsBond = DfsBond<decltype(source), decltype(target), std::remove_const_t<decltype(expr)>, isCyclic, edge.closure>{};

                return ctll::push_front(dfsBond, makeValueDfsBondListHelper<DfsEdgeIndex + 1>(smarts, dfsEdgeList, cycleMembership));
            }
        }

        consteval auto makeValueDfsBondList(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            return makeValueDfsBondListHelper<0>(smarts, dfsEdgeList, cycleMembership);
        }

    } // namespace impl

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct ValueDfsBondList
    {
        static constexpr inline auto data = impl::makeValueDfsBondList(SmartsT{}, DfsEdgeListT{}, CycleMembershipT{});

        consteval ValueDfsBondList() noexcept {}
        consteval ValueDfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
