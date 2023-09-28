#pragma once

#include "../Util/Util.hpp"

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<int SourceIndex, int AdjIndex>
        consteval auto typeDfsSearch(auto smarts, auto incidentList, auto visitor, auto state,
                                     auto visitedVertices, auto visitedEdges) noexcept
        {
            constexpr auto sourceDegree = incidentList.degrees.data[SourceIndex];
            if constexpr (AdjIndex < sourceDegree) {
                constexpr auto edgeIndex = incidentList.get(SourceIndex, AdjIndex);

                if constexpr (ctll::exists_in(Number<edgeIndex>{}, visitedEdges)) {
                    return typeDfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, state, visitedVertices, visitedEdges);
                } else {
                    constexpr auto edge = get<edgeIndex>(incidentList.edges.data);
                    constexpr auto targetIndex = edge.source == SourceIndex ? edge.target : edge.source;
                    constexpr auto isNewComponent = !ctll::exists_in(Number<SourceIndex>{}, visitedVertices);
                    constexpr auto isClosure = ctll::exists_in(Number<targetIndex>{}, visitedVertices);

                    auto state1 = visitor.template visit<edgeIndex, SourceIndex, targetIndex, isNewComponent, isClosure>(state);
                    auto vv1 = ctll::add_item(Number<targetIndex>{}, ctll::add_item(Number<SourceIndex>{}, visitedVertices));
                    auto ve1 = ctll::add_item(Number<edgeIndex>{}, visitedEdges);

                    // dfs
                    if constexpr (!isClosure) {
                        auto [state2, vv2, ve2] = typeDfsSearch<targetIndex, 0>(smarts, incidentList, visitor, state1, vv1, ve1);
                        auto state3 = visitor.template backtrack<edgeIndex, targetIndex, isClosure>(state2);
                        return typeDfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, state3, vv2, ve2);
                    } else {
                        auto state2 = visitor.template backtrack<edgeIndex, targetIndex, isClosure>(state1);
                        return typeDfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, state2, vv1, ve1);
                    }
                }
            } else if constexpr (SourceIndex == 0) {
                return visitor.template backtrack<SourceIndex>(state);
            } else {
                return std::make_tuple(state, visitedVertices, visitedEdges);
            }
        }

    } // namespace impl

    struct TypeDFSVisitorBase
    {
        template<int Edge, int Source, int Target, bool IsNewComponent, bool IsClosure>
        constexpr auto visit(auto state) noexcept
        {
            return state;
        }
        template<int Edge, int Target, bool IsClosure>
        constexpr auto backtrack(auto state) noexcept
        {
            return state;
        }
        template<int Source>
        constexpr auto backtrack(auto state) noexcept
        {
            return state;
        }
    };

    template<int SourceIndex = 0>
    consteval auto typeDfsSearch(auto smarts, auto incidentList, auto visitor, auto state) noexcept
    {        
        return impl::typeDfsSearch<SourceIndex, 0>(smarts, incidentList, visitor, state, ctll::empty_list{}, ctll::empty_list{});
    }



} // namespace Kitimar::CTSmarts
