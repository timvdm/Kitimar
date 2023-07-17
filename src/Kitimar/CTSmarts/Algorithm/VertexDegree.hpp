#pragma once

#include <algorithm>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename EdgeListT>
    constexpr auto makeVertexDegree() noexcept
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

        static constexpr inline int get(int index) noexcept
        {
            return data[index];
        }

        static constexpr auto max() noexcept
        {
            return *std::ranges::max_element(data);
        }

        constexpr VertexDegree() noexcept {}
        constexpr VertexDegree(SmartsT, EdgeListT) noexcept {}
    };

}
