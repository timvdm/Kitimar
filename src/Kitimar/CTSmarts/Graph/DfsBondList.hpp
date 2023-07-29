#pragma once

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<typename Source, typename Target, typename Expr, bool IsCyclic, bool IsRingClosure>
    struct DfsBond
    {
        static constexpr inline auto source = Source{};
        static constexpr inline auto target = Target{};
        static constexpr inline auto expr = Expr();
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
    };

    namespace impl {

        template<int DfsEdgeIndex>
        consteval auto makeDfsBondListHelper(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            if constexpr (DfsEdgeIndex == smarts.numBonds) {
                return ctll::empty_list{};
            } else {
                constexpr auto edge = dfsEdgeList.data[DfsEdgeIndex];
                constexpr auto expr = get<edge.index>(smarts.bonds).expr;
                auto source = get<edge.source>(smarts.atoms);
                auto target = get<edge.target>(smarts.atoms);
                constexpr auto isCyclic = cycleMembership.edges[edge.index];

                auto dfsBond = DfsBond<decltype(source), decltype(target), decltype(expr), isCyclic, edge.closure>{};

                return ctll::push_front(dfsBond, makeDfsBondListHelper<DfsEdgeIndex + 1>(smarts, dfsEdgeList, cycleMembership));
            }
        }

        consteval auto makeDfsBondList(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            return makeDfsBondListHelper<0>(smarts, dfsEdgeList, cycleMembership);
        }

    } // namespace impl

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct DfsBondList
    {
        static constexpr inline auto data = impl::makeDfsBondList(SmartsT{}, DfsEdgeListT{}, CycleMembershipT{});

        consteval DfsBondList() noexcept {}
        consteval DfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
