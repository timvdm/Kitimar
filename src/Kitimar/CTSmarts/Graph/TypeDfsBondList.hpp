#pragma once

#include "DfsBond.hpp"

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<int DfsEdgeIndex>
        consteval auto makeTypeDfsBondListHelper(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            if constexpr (DfsEdgeIndex == smarts.numBonds) {
                return ctll::empty_list{};
            } else {
                constexpr auto edge = get<DfsEdgeIndex>(dfsEdgeList.data);
                constexpr auto expr = get<edge.index>(smarts.bonds).expr;
                auto source = get<edge.source>(smarts.atoms);
                auto target = get<edge.target>(smarts.atoms);
                constexpr auto isCyclic = cycleMembership.edges[edge.index];

                auto dfsBond = DfsBond<decltype(source), decltype(target), std::remove_const_t<decltype(expr)>, isCyclic, edge.closure>{};

                return ctll::push_front(dfsBond, makeTypeDfsBondListHelper<DfsEdgeIndex + 1>(smarts, dfsEdgeList, cycleMembership));
            }
        }

        consteval auto makeTypeDfsBondList(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            return makeTypeDfsBondListHelper<0>(smarts, dfsEdgeList, cycleMembership);
        }

    } // namespace impl

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct TypeDfsBondList
    {
        static constexpr inline auto data = impl::makeTypeDfsBondList(SmartsT{}, DfsEdgeListT{}, CycleMembershipT{});

        consteval TypeDfsBondList() noexcept {}
        consteval TypeDfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
