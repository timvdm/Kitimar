#pragma once

#include "Primitives.hpp"
#include "Operations.hpp"
#include "BasicSmarts.hpp"

namespace Kitimar::CTSmarts {

    namespace impl {

        // Primitive
        template<typename Query, typename Expr>
        constexpr auto containsExprHelper(Query, Expr) noexcept
        {
            return std::is_same_v<Query, Expr>;
        }

        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, ctll::list<Expr...>) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, Or<Expr...>) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, And<Expr...>) noexcept;

        // Not
        template<typename Query, typename Expr>
        constexpr auto containsExprHelper(Query, Not<Expr>) noexcept
        {
            if constexpr (std::is_same_v<Query, Not<Expr>>)
                return true;
            else
                return containsExprHelper(Query{}, Expr{});
        }

        // Or
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, Or<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<Query, Or<Expr...>>)
                return true;
            else
                return (containsExpr(Query{}, Expr{}) || ...);
        }

        // And
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, And<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<Query, And<Expr...>>)
                return true;
            else
                return (containsExpr(Query{}, Expr{}) || ...);
        }

        // Atom
        template<typename Query, int Index, typename Expr>
        constexpr auto containsExprHelper(Query, Atom<Index, Expr>) noexcept
        {
            return containsExprHelper(Query{}, Expr{});
        }

        // ctll::list
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, ctll::list<Expr...>) noexcept
        {
            return (containsExpr(Query{}, Expr{}) || ...);
        }

    } // namespace impl

    template<typename Query, typename Expr>
    constexpr auto containsExpr(Query, Expr) noexcept
    {
        return impl::containsExprHelper(Query{}, Expr{});
    }

} // namespace Kitimar::CTSmarts
