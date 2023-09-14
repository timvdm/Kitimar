#pragma once

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<typename Expr>
    struct Not
    {
        static constexpr inline auto expr = Expr{};

        constexpr Not() noexcept {}
        constexpr Not(Expr) noexcept {}
    };

    template<typename ...Expr>
    struct And
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr And() noexcept {}
        constexpr And(Expr...) noexcept {}
        constexpr And(ctll::list<Expr...>) noexcept {}
    };

    template<typename ...Expr>
    struct Or
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr Or() noexcept {}
        constexpr Or(Expr...) noexcept {}
        constexpr Or(ctll::list<Expr...>) noexcept {}
    };

} // namespace Kitimar::CTSmarts
