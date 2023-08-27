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

    enum class SeedType {
        Molecule, // start from any atom/bond
        Atom, // start from first atom
        Bond // start from first bond
    };

    template<SeedType T>
    using SeedTypeTag = std::integral_constant<SeedType, T>;

    struct NoOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
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
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = EdgeList{smarts};
            constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
            constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};
            constexpr auto atomFreq = AtomFrequency{smarts};
            constexpr auto optimizedIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};
            constexpr auto sourceIndex = (SeedT == SeedType::Molecule) ? std::ranges::min_element(atomFreq.data) - std::begin(atomFreq.data) : 0;
            constexpr auto dfsEdges = DfsEdgeList{smarts, optimizedIncidentList, Number<sourceIndex>{}};
            constexpr auto cycleMembership = CycleMembership{smarts, optimizedIncidentList};
            constexpr auto dfsBonds = DfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

}
