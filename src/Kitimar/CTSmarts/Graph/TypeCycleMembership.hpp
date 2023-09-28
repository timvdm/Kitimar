#pragma once

#include <array>
#include <tuple>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename IncidentListT>
        struct TypeCycleMembershipVisitor : TypeDFSVisitorBase
        {

            consteval TypeCycleMembershipVisitor() noexcept {}
            consteval TypeCycleMembershipVisitor(SmartsT, IncidentListT) noexcept {}

            template<int Target, bool First = true>
            constexpr auto addPath(auto path, auto &vertices, auto &edges)
            {
                if constexpr (!ctll::empty(path)) {

                    //auto [edgeIndex, tail] = ctll::pop_and_get_front(path); // does not work with MSVC
                    auto edgeIndex = ctll::front(path);
                    auto tail = ctll::pop_front(path);
                    constexpr auto e = get<edgeIndex.value>(IncidentListT::edges.data);
                    edges[edgeIndex.value] = true;
                    vertices[e.source] = true;
                    vertices[e.target] = true;
                    constexpr auto done = !First && (e.source == Target || e.target == Target);
                    if constexpr (!ctll::empty(tail) && !done)
                        addPath<Target, false>(tail, vertices, edges);
                }
            }

            template<int Edge, int Source, int Target, bool IsNewComponent, bool IsClosure>
            constexpr auto visit(auto state) noexcept
            {
                auto [vertices, edges, path] = state;
                auto path2 = ctll::push_front(Number<Edge>{}, path);

                if constexpr (IsClosure)
                    addPath<Target>(path2, vertices, edges);

                return std::make_tuple(vertices, edges, path2);
            }

            template<int Edge, int Target, bool IsClosure>
            constexpr auto backtrack(auto state) noexcept
            {
                auto [vertices, edges, path] = state;
                return std::make_tuple(vertices, edges, ctll::pop_front(path));
            }

            template<int Source>
            constexpr auto backtrack(auto state) noexcept
            {
                return state;
            }
        };

        template<int SourceIndex = 0>
        consteval auto makeTypeCycleMembership(auto smarts, auto incidentList) noexcept
        {
            TypeCycleMembershipVisitor visitor{smarts, incidentList};
            auto state = std::make_tuple(std::array<bool, smarts.numAtoms>{}, std::array<bool, smarts.numBonds>{}, ctll::empty_list{});
            auto state2 = typeDfsSearch<SourceIndex>(smarts, incidentList, visitor, state);
            return std::make_tuple(std::get<0>(state2), std::get<1>(state2));
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct TypeCycleMembership
    {
        static constexpr inline auto data = impl::makeTypeCycleMembership(SmartsT{}, EdgeListT{});
        static constexpr inline auto vertices = std::get<0>(data);
        static constexpr inline auto edges = std::get<1>(data);

        consteval TypeCycleMembership() noexcept {}
        consteval TypeCycleMembership(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
