#pragma once

#include <array>
#include <algorithm>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto makeIncidentList(auto smarts, auto edgeList, auto degrees) noexcept
        {
            constexpr auto stride = *std::ranges::max_element(degrees.data);
            std::array<int, smarts.numAtoms * stride> incident = {};
            std::ranges::fill(incident, -1);

            std::array<int, smarts.numAtoms> sizes = {};
            for (auto i = 0UL; i < smarts.numBonds; ++i) {
                auto edge = edgeList.data[i];
                auto source = edge.source;
                auto target = edge.target;
                incident[stride * source + sizes[source]] = i;
                incident[stride * target + sizes[target]] = i;
                ++sizes[source];
                ++sizes[target];
            }

            return incident;
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
    struct IncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = impl::makeIncidentList(SmartsT{}, EdgeListT{}, VertexDegreeT{});
        static constexpr inline auto edges = EdgeListT{};
        static constexpr inline auto degrees = VertexDegreeT{};
        static constexpr inline auto stride = *std::ranges::max_element(VertexDegreeT::data);

        static consteval auto get(int VertexIndex, int IncidentIndex) noexcept
        {
            return data[stride * VertexIndex + IncidentIndex];
        }

        consteval IncidentList() noexcept {}
        consteval IncidentList(SmartsT, EdgeListT, VertexDegreeT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
