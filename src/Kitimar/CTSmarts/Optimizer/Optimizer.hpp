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

#if _MSC_FULL_VER < 193732824 // MSVC 19.37.32824

    struct NoOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = ValueEdgeList{smarts};
            constexpr auto vertexDegree = ValueVertexDegree{smarts, edgeList};
            constexpr auto incidentList = ValueIncidentList{smarts, edgeList, vertexDegree};
            constexpr auto dfsEdges = ValueDfsEdgeList{smarts, incidentList};
            constexpr auto cycleMembership = ValueCycleMembership{smarts, incidentList};
            constexpr auto dfsBonds = ValueDfsBondList{smarts, dfsEdges, cycleMembership};
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
            constexpr auto edgeList = ValueEdgeList{smarts};
            constexpr auto vertexDegree = ValueVertexDegree{smarts, edgeList};
            constexpr auto incidentList = ValueIncidentList{smarts, edgeList, vertexDegree};
            constexpr auto atomFreq = AtomFrequency{smarts};
            constexpr auto optimizedIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};
            constexpr auto sourceIndex = (SeedT == SeedType::Molecule) ? std::ranges::min_element(atomFreq.data) - std::begin(atomFreq.data) : 0;
            constexpr auto dfsEdges = ValueDfsEdgeList{smarts, optimizedIncidentList, Number<sourceIndex>{}};
            constexpr auto cycleMembership = ValueCycleMembership{smarts, optimizedIncidentList};
            constexpr auto dfsBonds = ValueDfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

#else

    struct NoOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = TypeEdgeList{smarts};
            constexpr auto vertexDegree = TypeVertexDegree{smarts, edgeList};
            constexpr auto incidentList = TypeIncidentList{smarts, edgeList, vertexDegree};
            constexpr auto dfsEdges = TypeDfsEdgeList{smarts, incidentList};
            constexpr auto cycleMembership = TypeCycleMembership{smarts, incidentList};
            constexpr auto dfsBonds = TypeDfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

    // FIXME: optimize atom expressions...

    using FullOptimizer = NoOptimizer; // FIXME: TypeOptimizeIncidentList
    /*
    struct FullOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = TypeEdgeList{smarts};
            constexpr auto vertexDegree = TypeVertexDegree{smarts, edgeList};
            constexpr auto incidentList = TypeIncidentList{smarts, edgeList, vertexDegree};
            constexpr auto atomFreq = AtomFrequency{smarts};
            constexpr auto optimizedIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};
            constexpr auto sourceIndex = (SeedT == SeedType::Molecule) ? std::ranges::min_element(atomFreq.data) - std::begin(atomFreq.data) : 0;
            constexpr auto dfsEdges = TypeDfsEdgeList{smarts, optimizedIncidentList, Number<sourceIndex>{}};
            constexpr auto cycleMembership = TypeCycleMembership{smarts, optimizedIncidentList};
            constexpr auto dfsBonds = TypeDfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };
    */

#endif

}
