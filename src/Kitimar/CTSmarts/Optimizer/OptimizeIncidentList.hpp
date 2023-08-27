#pragma once

#include <algorithm>
#include <type_traits>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename IncidentListT, typename VertexFrequencyT, bool SeedBond>
    consteval auto makeOptimizeIncidentList()
    {
        auto incident = IncidentListT::data;

        for (auto i = 0; i < static_cast<int>(SmartsT::numAtoms); ++i) {
            auto offset = i * IncidentListT::stride;
            auto end = offset + IncidentListT::degrees.data[i];
            if (SeedBond)
                ++offset;
            std::ranges::sort(std::begin(incident) + offset, std::begin(incident) + end, [i] (auto index1, auto index2) {
                auto edge1 = IncidentListT::edges.data[index1];
                auto edge2 = IncidentListT::edges.data[index2];
                auto target1 = i == edge1.source ? edge1.target : edge1.source;
                auto target2 = i == edge2.source ? edge2.target : edge2.source;
                return VertexFrequencyT::data[target1] < VertexFrequencyT::data[target2];
            });
        }

        return incident;
    }

    template<typename SmartsT, typename IncidentListT, typename VertexFrequencyT, bool SeedBond>
    struct OptimizeIncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = makeOptimizeIncidentList<SmartsT, IncidentListT, VertexFrequencyT, SeedBond>();
        static constexpr inline auto edges = IncidentListT::edges;
        static constexpr inline auto degrees = IncidentListT::degrees;
        static constexpr inline auto stride = IncidentListT::stride;

        static consteval auto get(int AtomIndex, int AdjIndex)
        {
            return data[stride * AtomIndex + AdjIndex];
        }

        consteval OptimizeIncidentList() noexcept {}
        consteval OptimizeIncidentList(SmartsT, IncidentListT, VertexFrequencyT, std::integral_constant<bool, SeedBond>) noexcept {}
    };

}
