#pragma once

#include <ranges>

namespace Kitimar::CTSmarts {


    template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
    constexpr auto makeIncidentList()
    {
        constexpr auto stride = VertexDegreeT::max();
        std::array<int, SmartsT::numAtoms * stride> adj = {};

        for (auto &i : adj)
            i = -1; // FIXME: needed?

        std::array<int, SmartsT::numAtoms> sizes = {};
        for (auto i = 0; i < SmartsT::numBonds; ++i) {
            auto edge = EdgeListT::get(i);
            auto source = edge.source;
            auto target = edge.target;
            adj[stride * source + sizes[source]] = i;
            adj[stride * target + sizes[target]] = i;
            ++sizes[source];
            ++sizes[target];
        }

        return adj;
    }

    template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
    struct IncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
            static constexpr inline auto data = makeIncidentList<SmartsT, EdgeListT, VertexDegreeT>();
        static constexpr inline auto edges = EdgeListT{};
        static constexpr inline auto degrees = VertexDegreeT{};
        static constexpr inline auto stride = VertexDegreeT::max();

        static constexpr auto get(int AtomIndex, int AdjIndex)
        {
            return data[stride * AtomIndex + AdjIndex];
        }

        constexpr IncidentList() noexcept {}
        constexpr IncidentList(SmartsT, EdgeListT, VertexDegreeT) noexcept {}
    };

}
