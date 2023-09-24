#pragma once

#include "DfsSearchEventType.hpp"

#include <array>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>
#endif

namespace Kitimar::CTSmarts {

    struct ValueDfsSearchEvent
    {
        DfsSearchEventType type = DfsSearchEventType::Invalid;
        int index = -1;
        bool flag = false;

        constexpr auto operator<=>(const ValueDfsSearchEvent&) const noexcept = default;
    };

    #ifdef KITIMAR_WITH_IOSTREAM

    inline std::ostream& operator<<(std::ostream &os, const ValueDfsSearchEvent &event)
    {
        switch (event.type) {
            case DfsSearchEventType::VisitVertex:
                os << "VisitVertex( index = " << event.index << ", isNewComponent = " << event.flag << " )";
                break;
            case DfsSearchEventType::VisitEdge:
                os << "VisitEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                break;
            case DfsSearchEventType::BacktrackVertex:
                os << "BacktrackVertex( index = " << event.index << ", isEndComponent = " << event.flag << " )";
                break;
            case DfsSearchEventType::BacktrackEdge:
                os << "BacktrackEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                break;
            default:
                os << "InvalidEvent";
                break;
        }

        return os;
    }

    #endif // KITIMAR_WITH_IOSTREAM

    namespace impl {

        template<typename SmartsT>
        struct ValueDfsSearchEventsVisitor
        {
            using Events = std::array<ValueDfsSearchEvent, 2 * (SmartsT::numAtoms + SmartsT::numBonds)>;

            consteval ValueDfsSearchEventsVisitor() noexcept {}
            consteval ValueDfsSearchEventsVisitor(SmartsT) noexcept {}

            constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
            {
                if (isNewComponent)
                    events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, source, true};
                events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, edge, isClosure};
                if (!isClosure)
                    events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, target, false};
            }

            constexpr void backtrack(int edge, int target, bool isClosure) noexcept
            {
                if (!isClosure)
                    events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, target, false};
                events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, edge, isClosure};
            }

            constexpr void backtrack(int source) noexcept
            {
                events[nextEventIndex++] = ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, source, true};
            }

            Events events = {};
            int nextEventIndex = 0;
        };

        consteval auto makeValueDfsSearchEvents(auto smarts, auto incidentList) noexcept
        {
            ValueDfsSearchEventsVisitor visitor{smarts};
            valueDfsSearch(smarts, incidentList, visitor);
            return visitor.events;
        }

    } // namespace impl


    template<typename SmartsT, typename IncidentListT>
    struct ValueDfsSearchEvents
    {
        static constexpr inline auto events = impl::makeValueDfsSearchEvents(SmartsT{}, IncidentListT{});

        consteval ValueDfsSearchEvents() noexcept {}
        consteval ValueDfsSearchEvents(SmartsT, IncidentListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
