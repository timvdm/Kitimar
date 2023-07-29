#pragma once

#include "EdgeList.hpp"
#include "DfsSearch.hpp"
#include "../Util/Util.hpp"

namespace Kitimar::CTSmarts {

    struct DfsEdge : Edge
    {
        bool closure = false;
    };

    namespace impl {

        template<typename SmartsT>
        struct DfsEdgeListVisitor : DFSVisitorBse
        {
            consteval DfsEdgeListVisitor() noexcept {}
            consteval DfsEdgeListVisitor(SmartsT) noexcept {}

            constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
            {
                edges[nextEdgeIndex++] = DfsEdge{{edge, source, target}, isClosure};
            }

            std::array<DfsEdge, SmartsT::numBonds> edges;
            int nextEdgeIndex = 0;
        };

        template<int SourceIndex = 0>
        consteval auto makeDfsEdgeList(auto smarts, auto incidentList) noexcept
        {
            DfsEdgeListVisitor visitor{smarts};
            dfsSearch<SourceIndex>(smarts, incidentList, visitor);
            return visitor.edges;
        }

    } // namespace impl

    template<typename SmartsT, typename IncidentListT, int SourceIndex = 0>
    struct DfsEdgeList
    {
        static constexpr inline auto data = impl::makeDfsEdgeList<SourceIndex>(SmartsT{}, IncidentListT{});

        consteval DfsEdgeList() noexcept {}
        consteval DfsEdgeList(SmartsT, IncidentListT, Number<SourceIndex> = {}) noexcept {}
    };

} // namespace Kitimar::CTSmarts
