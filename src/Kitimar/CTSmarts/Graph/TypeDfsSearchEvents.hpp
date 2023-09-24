#pragma once

#include "DfsSearchEventType.hpp"

#include <array>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>
#endif

namespace Kitimar::CTSmarts {

    template<DfsSearchEventType Type, int Index, bool Flag>
    struct TypeDfsSearchEvent
    {
        static constexpr inline auto type = Type;
        static constexpr inline auto index = Index;
        static constexpr inline auto flag = Flag;
    };

    namespace impl {

        template<typename SmartsT>
        struct TypeDfsSearchEventsVisitor
        {
            consteval TypeDfsSearchEventsVisitor() noexcept {}
            consteval TypeDfsSearchEventsVisitor(SmartsT) noexcept {}

            template<int Source, bool IsNewComponent>
            constexpr auto visitSource(auto state) noexcept
            {
                if constexpr (IsNewComponent)
                    return ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, Source, true>{}, state);
                else
                    return state;
            }

            template<int Target, bool IsClosure>
            constexpr auto visitTarget(auto state) noexcept
            {
                if constexpr (!IsClosure)
                    return ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, Target, false>{}, state);
                else
                    return state;
            }

            template<int Edge, int Source, int Target, bool IsNewComponent, bool IsClosure>
            constexpr auto visit(auto state) noexcept
            {
                auto state1 = visitSource<Source, IsNewComponent>(state);
                auto state2 = ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, Edge, IsClosure>{}, state1);
                return visitTarget<Target, IsClosure>(state2);
            }

            template<int Target, bool IsClosure>
            constexpr auto backtrackTarget(auto state) noexcept
            {
                if constexpr (!IsClosure)
                    return ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, Target, false>{}, state);
                else
                    return state;
            }

            template<int Edge, int Target, bool IsClosure>
            constexpr auto backtrack(auto state) noexcept
            {
                auto state1 = backtrackTarget<Target, IsClosure>(state);
                return ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, Edge, IsClosure>{}, state1);
            }

            template<int Source>
            constexpr auto backtrack(auto state) noexcept
            {
                return ctll::push_front(TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, Source, true>{}, state);
            }
        };

        consteval auto makeTypeDfsSearchEvents(auto smarts, auto incidentList) noexcept
        {
            TypeDfsSearchEventsVisitor visitor{smarts};
            return ctll::rotate(typeDfsSearch(smarts, incidentList, visitor, ctll::empty_list{}));
        }

    } // namespace impl


    template<typename SmartsT, typename IncidentListT>
    struct TypeDfsSearchEvents
    {
        static constexpr inline auto events = impl::makeTypeDfsSearchEvents(SmartsT{}, IncidentListT{});

        consteval TypeDfsSearchEvents() noexcept {}
        consteval TypeDfsSearchEvents(SmartsT, IncidentListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
