#pragma once

#include <array>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>
#endif

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT>
        struct DfsSearchEventsVisitor
        {
            struct Event
            {
                enum Type {
                    Invalid,
                    VisitVertex, // index = vertex index, flag = isNewComponent
                    VisitEdge, // index = edge index, flag = isClosure
                    BacktrackVertex, // index = vertex index, flag = isEndComponent
                    BacktrackEdge // index = edge index, flag = isClosure
                };

                Type type = Invalid;
                int index = -1;
                bool flag = false;
            };

            #ifdef KITIMAR_WITH_IOSTREAM
            friend std::ostream& operator<<(std::ostream &os, const Event &event)
            {
                switch (event.type) {
                    case Event::VisitVertex:
                        os << "VisitVertex( index = " << event.index << ", isNewComponent = " << event.flag << " )";
                        break;
                    case Event::VisitEdge:
                        os << "VisitEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                        break;
                    case Event::BacktrackVertex:
                        os << "BacktrackVertex( index = " << event.index << ", isEndComponent = " << event.flag << " )";
                        break;
                    case Event::BacktrackEdge:
                        os << "BacktrackEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                        break;
                    default:
                        os << "InvalidEvent";
                        break;
                }

                return os;
            }
            #endif // KITIMAR_WITH_IOSTREAM

            using Events = std::array<Event, 2 * (SmartsT::numAtoms + SmartsT::numBonds)>;

            consteval DfsSearchEventsVisitor() noexcept {}
            consteval DfsSearchEventsVisitor(SmartsT) noexcept {}

            constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
            {
                if (isNewComponent)
                    events[nextEventIndex++] = Event{Event::VisitVertex, source, true};
                events[nextEventIndex++] = Event{Event::VisitEdge, edge, isClosure};
                if (!isClosure)
                    events[nextEventIndex++] = Event{Event::VisitVertex, target, false};
            }

            constexpr void backtrack(int edge, int target, bool isClosure) noexcept
            {
                if (!isClosure)
                    events[nextEventIndex++] = Event{Event::BacktrackVertex, target, false};
                events[nextEventIndex++] = Event{Event::BacktrackEdge, edge, isClosure};
            }

            constexpr void backtrack(int source) noexcept
            {
                events[nextEventIndex++] = Event{Event::BacktrackVertex, source, true};
            }

            Events events = {};
            int nextEventIndex = 0;
        };

        consteval auto makeDfsSearchEvents(auto smarts, auto incidentList) noexcept
        {
            DfsSearchEventsVisitor visitor{smarts};
            dfsSearch(smarts, incidentList, visitor);
            return visitor.events;
        }

    } // namespace impl


    template<typename SmartsT, typename IncidentListT>
    struct DfsSearchEvents
    {
        static constexpr inline auto events = impl::makeDfsSearchEvents(SmartsT{}, IncidentListT{});

        consteval DfsSearchEvents() noexcept {}
        consteval DfsSearchEvents(SmartsT, IncidentListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
