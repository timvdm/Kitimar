#pragma once

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    template<int Index, int Source, int Target>
    struct TypeEdge
    {
        static constexpr inline auto index = Index;
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;

        //constexpr auto operator<=>(const Edge&) const noexcept = default;
    };

    namespace impl {

        consteval auto makeTypeEdgeList(auto smarts, auto bonds) noexcept
        {
            if constexpr (ctll::empty(bonds))
                return ctll::empty_list{};
            else {
                auto [bond, tail] = ctll::pop_and_get_front(bonds);
                auto edges = makeTypeEdgeList(smarts, tail);
                return ctll::push_front(TypeEdge<bond.index, bond.source, bond.target>{}, edges);
            }
        }

    } // namespace impl

    template<typename SmartsT>
    struct TypeEdgeList
    {
        static constexpr inline auto data = impl::makeTypeEdgeList(SmartsT{}, SmartsT::bonds);

        consteval TypeEdgeList() noexcept {}
        consteval TypeEdgeList(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
