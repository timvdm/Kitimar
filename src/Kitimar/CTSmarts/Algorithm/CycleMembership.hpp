#pragma once

#include <array>
#include <tuple>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename IncidentListT>
    struct CycleMembershipVisitor
    {

        consteval CycleMembershipVisitor() noexcept {}
        consteval CycleMembershipVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            path[depth++] = edge;

            if (isClosure) {
                for (auto i = depth - 1; i >= 0; --i) {                    
                    auto e = IncidentListT::edges.data[path[i]];
                    edges[path[i]] = true;
                    vertices[e.source] = true;
                    vertices[e.target] = true;
                    if (i < depth - 1)
                        if (e.source == target || e.target == target)
                            break;
                }
            }
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept
        {
            --depth;
        }

        constexpr void backtrack(int source) noexcept {}

        std::array<int, SmartsT::numBonds> path = {};
        std::array<bool, SmartsT::numAtoms> vertices = {};
        std::array<bool, SmartsT::numBonds> edges = {};
        int depth = 0;
    };

    template<typename SmartsT, typename IncidentListT>
    consteval auto makeCycleMembership()
    {
        CycleMembershipVisitor<SmartsT, IncidentListT> visitor;
        dfsSearch(SmartsT{}, visitor, IncidentListT{});
        return std::make_tuple(visitor.vertices, visitor.edges);
    }

    template<typename SmartsT, typename EdgeListT>
    struct CycleMembership
    {
        static constexpr inline auto data = makeCycleMembership<SmartsT, EdgeListT>();
        static constexpr inline auto vertices = std::get<0>(data);
        static constexpr inline auto edges = std::get<1>(data);

        consteval CycleMembership() noexcept {}
        consteval CycleMembership(SmartsT, EdgeListT) noexcept {}
    };

}
