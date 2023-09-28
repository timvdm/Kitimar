#pragma once

#include <array>
#include <algorithm>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto makeTypeIncidentList(auto smarts, auto edges, auto degrees) noexcept
        {
            constexpr std::size_t stride = *std::ranges::max_element(degrees.data);
            if constexpr (ctll::empty(edges)) {
                auto incident = std::array<int, ctll::size(smarts.atoms) * stride>{};
                std::ranges::fill(incident, -1);
                std::array<int, smarts.numAtoms> sizes = {};
                return std::make_tuple(incident, sizes);
            } else {
                auto [edge, tail] = ctll::pop_and_get_front(edges);
                auto [incident, sizes] = makeTypeIncidentList(smarts, tail, degrees);
                auto source = edge.source;
                auto target = edge.target;
                incident[stride * source + degrees.data[source] - sizes[source] - 1] = edge.index;
                incident[stride * target + degrees.data[target] - sizes[target] - 1] = edge.index;
                ++sizes[source];
                ++sizes[target];
                return std::make_tuple(incident, sizes);
            }



            /*
            constexpr std::size_t stride = *std::ranges::max_element(VertexDegreeT::data);
            std::array<int, SmartsT::numAtoms * stride> incident = {};
            std::ranges::fill(incident, -1);

            std::array<int, SmartsT::numAtoms> sizes = {};
            for (auto i = 0UL; i < SmartsT::numBonds; ++i) {
                auto edge = EdgeListT::data[i];
                auto source = edge.source;
                auto target = edge.target;
                incident[stride * source + sizes[source]] = i;
                incident[stride * target + sizes[target]] = i;
                ++sizes[source];
                ++sizes[target];
            }

            return incident;
            */
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
    struct TypeIncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = std::get<0>(impl::makeTypeIncidentList(SmartsT{}, EdgeListT::data, VertexDegreeT{}));
        static constexpr inline auto edges = EdgeListT{};
        static constexpr inline auto degrees = VertexDegreeT{};
        static constexpr inline auto stride = *std::ranges::max_element(VertexDegreeT::data);

        static consteval auto get(int VertexIndex, int IncidentIndex) noexcept
        {
            return data[stride * VertexIndex + IncidentIndex];
        }

        consteval TypeIncidentList() noexcept {}
        consteval TypeIncidentList(SmartsT, EdgeListT, VertexDegreeT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
