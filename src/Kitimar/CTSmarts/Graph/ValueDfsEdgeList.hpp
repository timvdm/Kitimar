#pragma once

#include "ValueEdgeList.hpp"
#include "ValueDfsSearch.hpp"
#include "../Util/Util.hpp"

namespace Kitimar::CTSmarts {

    struct ValueDfsEdge : ValueEdge
    {
        bool closure = false;
    };

    namespace impl {

        template<typename SmartsT>
        struct ValueDfsEdgeListVisitor : ValueDFSVisitorBase
        {
            consteval ValueDfsEdgeListVisitor() noexcept {}
            consteval ValueDfsEdgeListVisitor(SmartsT) noexcept {}

            constexpr void visit(int edge, int source, int target, [[maybe_unused]] bool isNewComponent, bool isClosure) noexcept
            {
                edges[nextEdgeIndex++] = ValueDfsEdge{{edge, source, target}, isClosure};
            }

            std::array<ValueDfsEdge, SmartsT::numBonds> edges;
            int nextEdgeIndex = 0;
        };

        template<int SourceIndex = 0>
        consteval auto makeValueDfsEdgeList(auto smarts, auto incidentList) noexcept
        {
            ValueDfsEdgeListVisitor visitor{smarts};
            valueDfsSearch<SourceIndex>(smarts, incidentList, visitor);
            return visitor.edges;
        }

    } // namespace impl

    template<typename SmartsT, typename IncidentListT, int SourceIndex = 0>
    struct ValueDfsEdgeList
    {
        static constexpr inline auto data = impl::makeValueDfsEdgeList<SourceIndex>(SmartsT{}, IncidentListT{});

        consteval ValueDfsEdgeList() noexcept {}
        consteval ValueDfsEdgeList(SmartsT, IncidentListT, Number<SourceIndex> = {}) noexcept {}
    };

} // namespace Kitimar::CTSmarts
