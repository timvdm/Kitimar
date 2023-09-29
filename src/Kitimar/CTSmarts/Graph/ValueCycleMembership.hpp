#pragma once

#include <array>
#include <tuple>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename IncidentListT>
        struct ValueCycleMembershipVisitor
        {

            consteval ValueCycleMembershipVisitor() noexcept {}
            consteval ValueCycleMembershipVisitor(SmartsT, IncidentListT) noexcept {}

            constexpr void visit(int edge, [[maybe_unused]] int source, int target, [[maybe_unused]] bool isNewComponent, bool isClosure) noexcept
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

            constexpr void backtrack([[maybe_unused]] int edge, [[maybe_unused]] int target, [[maybe_unused]] bool isClosure) noexcept
            {
                --depth;
            }

            constexpr void backtrack([[maybe_unused]] int source) noexcept {}

            std::array<int, SmartsT::numBonds> path = {};
            std::array<bool, SmartsT::numAtoms> vertices = {};
            std::array<bool, SmartsT::numBonds> edges = {};
            int depth = 0;
        };

        consteval auto makeValueCycleMembership(auto smarts, auto incidentList) noexcept
        {
            ValueCycleMembershipVisitor visitor{smarts, incidentList};
            valueDfsSearch(smarts, incidentList, visitor);
            return std::make_tuple(visitor.vertices, visitor.edges);
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct ValueCycleMembership
    {
        static constexpr inline auto data = impl::makeValueCycleMembership(SmartsT{}, EdgeListT{});
        static constexpr inline auto vertices = std::get<0>(data);
        static constexpr inline auto edges = std::get<1>(data);

        consteval ValueCycleMembership() noexcept {}
        consteval ValueCycleMembership(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
