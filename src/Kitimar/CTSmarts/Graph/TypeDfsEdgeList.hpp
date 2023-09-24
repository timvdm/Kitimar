#pragma once

#include "TypeEdgeList.hpp"
#include "TypeDfsSearch.hpp"
#include "../Util/Util.hpp"

namespace Kitimar::CTSmarts {

    template<int Index, int Source, int Target, bool Closure>
    struct TypeDfsEdge : TypeEdge<Index, Source, Target>
    {
        static constexpr inline auto closure = Closure;
    };

    namespace impl {

        template<typename SmartsT>
        struct TypeDfsEdgeListVisitor : TypeDFSVisitorBase
        {
            consteval TypeDfsEdgeListVisitor() noexcept {}
            consteval TypeDfsEdgeListVisitor(SmartsT) noexcept {}

            template<int Edge, int Source, int Target, bool IsNewComponent, bool IsClosure>
            constexpr auto visit(auto state) noexcept
            {
                return ctll::push_front(TypeDfsEdge<Edge, Source, Target, IsClosure>{}, state);
            }
        };

        template<int SourceIndex = 0>
        consteval auto makeTypeDfsEdgeList(auto smarts, auto incidentList) noexcept
        {
            TypeDfsEdgeListVisitor visitor{smarts};
            return ctll::rotate(typeDfsSearch<SourceIndex>(smarts, incidentList, visitor, ctll::empty_list{}));
        }

    } // namespace impl

    template<typename SmartsT, typename IncidentListT, int SourceIndex = 0>
    struct TypeDfsEdgeList
    {
        static constexpr inline auto data = impl::makeTypeDfsEdgeList<SourceIndex>(SmartsT{}, IncidentListT{});

        consteval TypeDfsEdgeList() noexcept {}
        consteval TypeDfsEdgeList(SmartsT, IncidentListT, Number<SourceIndex> = {}) noexcept {}
    };

} // namespace Kitimar::CTSmarts
