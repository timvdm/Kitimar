#pragma once

#include <array>
#include <algorithm>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename EdgeListT>
    consteval auto makeVertexDegree() noexcept
    {
        std::array<int, SmartsT::numAtoms> degrees = {};
        for (const auto &edge : EdgeListT::data) {
            ++degrees[edge.source];
            ++degrees[edge.target];
        }
        return degrees;
    }

    template<typename SmartsT, typename EdgeListT>
    struct VertexDegree
    {
        static constexpr inline auto data = makeVertexDegree<SmartsT, EdgeListT>();

        static consteval auto max() noexcept
        {
            return *std::ranges::max_element(data);
        }

        consteval VertexDegree() noexcept {}
        consteval VertexDegree(SmartsT, EdgeListT) noexcept {}
    };

}
