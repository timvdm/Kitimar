#pragma once

#include <array>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename EdgeListT>
        consteval auto makeValueVertexDegree(SmartsT, EdgeListT) noexcept
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
    struct ValueVertexDegree
    {
        static constexpr inline auto data = impl::makeValueVertexDegree(SmartsT{}, EdgeListT{});

        consteval ValueVertexDegree() noexcept {}
        consteval ValueVertexDegree(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
