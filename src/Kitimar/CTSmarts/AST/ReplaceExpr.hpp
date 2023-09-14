#pragma once

#include "Primitives.hpp"
#include "Operations.hpp"
#include "BasicSmarts.hpp"

namespace Kitimar::CTSmarts {

    namespace impl {

        // Primitive
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Expr) noexcept
        {
            if constexpr (std::is_same_v<Expr, OldExpr>)
                return NewExpr{};
            else
                return Expr{};
        }

        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, ctll::list<Expr...> l) noexcept;

        // Not
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Not<Expr>) noexcept
        {
            return Not{replaceExprHelper(OldExpr{}, NewExpr{}, Expr{})};
        }

        // Or
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Or<Expr...> op) noexcept
        {
            return Or{replaceExprHelper(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // And
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, And<Expr...> op) noexcept
        {
            return And{replaceExprHelper(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // Atom
        template<typename OldExpr, typename NewExpr, int Index, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Atom<Index, Expr>) noexcept
        {
            return Atom<Index, decltype(replaceExprHelper(OldExpr{}, NewExpr{}, Expr{}))>{};
        }

        // ctll::list
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, ctll::list<Expr...> l) noexcept
        {
            return decltype(transform(l, [] (auto expr) { return replaceExprHelper(OldExpr{}, NewExpr{}, expr); })){};
        }

    } // namespace impl

    template<typename OldExpr, typename NewExpr, typename Expr>
    constexpr auto replaceExpr(OldExpr, NewExpr, Expr) noexcept
    {
        return impl::replaceExprHelper(OldExpr{}, NewExpr{}, Expr{});
    }

} // namespace Kitimar::CTSmarts
