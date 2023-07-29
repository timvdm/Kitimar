#pragma once

#include <array>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto makeVertexDegree(auto smarts, auto edgeList) noexcept
        {
            std::array<int, smarts.numAtoms> degrees = {};
            for (const auto &edge : edgeList.data) {
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
