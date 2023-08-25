#pragma once


namespace Kitimar::CTSmarts {

    template<typename VertexDegree, typename DfsBondList>
    struct SmartsQuery
    {
        static constexpr inline auto degrees = VertexDegree{}.data;
        static constexpr inline auto bonds = DfsBondList{}.data;

        consteval SmartsQuery() noexcept = default;
        consteval SmartsQuery(VertexDegree, DfsBondList) noexcept {};
    };

    struct NoOptimizer
    {
        template<typename SmartsT>
        static constexpr auto create(SmartsT) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = EdgeList{smarts};
            constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
            constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};
            constexpr auto dfsEdges = DfsEdgeList{smarts, incidentList};
            constexpr auto cycleMembership = CycleMembership{smarts, incidentList};
            constexpr auto dfsBonds = DfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

    // FIXME: optimize atom expressions...
    struct FullOptimizer
    {
        template<typename SmartsT>
        static constexpr auto create(SmartsT) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = EdgeList{smarts};
            constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
            constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};
            constexpr auto atomFreq = AtomFrequency{smarts};
            constexpr auto optimizedIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq};
            constexpr auto sourceIndex = std::ranges::min_element(atomFreq.data) - std::begin(atomFreq.data);
            constexpr auto dfsEdges = DfsEdgeList{smarts, optimizedIncidentList, Number<sourceIndex>{}};
            constexpr auto cycleMembership = CycleMembership{smarts, optimizedIncidentList};
            constexpr auto dfsBonds = DfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

}
