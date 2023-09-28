#pragma once

#include <array>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename EdgeListT>
        consteval auto makeVertexDegree(SmartsT, EdgeListT) noexcept
        {
            std::array<int, SmartsT::numAtoms> degrees = {};
            for (const auto &edge : EdgeListT::data) {
                ++degrees[edge.source];
                ++degrees[edge.target];
            }
            return degrees;
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct VertexDegree
    {
        static constexpr inline auto data = impl::makeVertexDegree(SmartsT{}, EdgeListT{});

        consteval VertexDegree() noexcept {}
        consteval VertexDegree(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
