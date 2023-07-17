#pragma once

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<int Source, int Target, bool IsCyclic, bool IsRingClosure, typename SourceExpr, typename TargetExpr, typename BondExpr>
    struct DfsBond
    {
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
        static constexpr inline auto sourceExpr = SourceExpr();
        static constexpr inline auto targetExpr = TargetExpr();
        static constexpr inline auto bondExpr = BondExpr();
    };

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT, int DfsEdgeIndex>
    constexpr auto makeDfsBondListHelper()
    {
        if constexpr (DfsEdgeIndex == SmartsT::numBonds) {
            return ctll::empty_list{};
        } else {
            constexpr auto edge = DfsEdgeListT::get(DfsEdgeIndex);

            constexpr auto bondExpr = get<edge.index>(SmartsT::bonds).expr;
            auto sourceExpr = get<edge.source>(SmartsT::atoms);
            auto targetExpr = get<edge.target>(SmartsT::atoms);
            constexpr auto isCyclic = CycleMembershipT::edges[edge.index];

            auto dfsBond = DfsBond<edge.source, edge.target, isCyclic, edge.closure, decltype(sourceExpr), decltype(targetExpr), decltype(bondExpr)>{};

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