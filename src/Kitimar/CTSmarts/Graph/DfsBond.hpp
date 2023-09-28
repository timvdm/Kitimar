#pragma once

namespace Kitimar::CTSmarts {

    template<typename Source, typename Target, typename Expr, bool IsCyclic, bool IsRingClosure>
    struct DfsBond
    {
        static constexpr inline auto source = Source{};
        static constexpr inline auto target = Target{};
        static constexpr inline auto expr = Expr{};
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
    };

} // namespace Kitimar::CTSmarts
