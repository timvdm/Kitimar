#pragma once

#include <array>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto makeTypeVertexDegree(auto smarts, auto edges) noexcept
        {
            if constexpr (ctll::empty(edges))
                return std::array<int, smarts.numAtoms>{};
            else {
                auto [edge, tail] = ctll::pop_and_get_front(edges);
                auto degrees = makeTypeVertexDegree(smarts, tail);
                ++degrees[edge.source];
                ++degrees[edge.target];
                return degrees;
            }
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct TypeVertexDegree
    {
        static constexpr inline auto data = impl::makeTypeVertexDegree(SmartsT{}, EdgeListT::data);

        consteval TypeVertexDegree() noexcept {}
        consteval TypeVertexDegree(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
