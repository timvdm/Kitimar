#pragma once

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<typename Source, typename Target, typename BondExpr, bool IsCyclic, bool IsRingClosure>
    struct DfsBond
    {
        static constexpr inline auto source = Source{};
        static constexpr inline auto target = Target{};
        static constexpr inline auto expr = BondExpr();
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
    };

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT, int DfsEdgeIndex>
    constexpr auto makeDfsBondListHelper()
    {
        if constexpr (DfsEdgeIndex == SmartsT::numBonds) {
            return ctll::empty_list{};
        } else {
            constexpr auto edge = DfsEdgeListT::get(DfsEdgeIndex);

            constexpr auto expr = get<edge.index>(SmartsT::bonds).expr;
            auto source = get<edge.source>(SmartsT::atoms);
            auto target = get<edge.target>(SmartsT::atoms);
            constexpr auto isCyclic = CycleMembershipT::edges[edge.index];

            auto dfsBond = DfsBond<decltype(source), decltype(target), decltype(expr), isCyclic, edge.closure>{};

            return ctll::push_front(dfsBond, makeDfsBondListHelper<SmartsT, DfsEdgeListT, CycleMembershipT, DfsEdgeIndex + 1>());
        }
    }

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    constexpr auto makeDfsBondList()
    {
        return makeDfsBondListHelper<SmartsT, DfsEdgeListT, CycleMembershipT, 0>();
    }


    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct DfsBondList
    {
        static constexpr inline auto data = makeDfsBondList<SmartsT, DfsEdgeListT, CycleMembershipT>();

        static constexpr inline auto get(int index) noexcept
        {
            return data[index];
        }

        constexpr DfsBondList() noexcept {}
        constexpr DfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };

}
