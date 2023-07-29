#pragma once

#include "PrimitiveFrequency.hpp"
#include "../Util/CtllSort.hpp"

namespace Kitimar::CTSmarts {

    struct ProjExprFrequency;

    template<typename Expr>
    consteval double expressionFrequency(Not<Expr> op) noexcept
    {
        return 1 - expressionFrequency(op.expr);
    }

    template<typename ...Expr>
    consteval double expressionFrequency(And<Expr...> op) noexcept
    {        
        return expressionFrequency(selectLast<ProjExprFrequency, std::greater<>>(op.expr));
    }

    template<typename ...Expr>
    consteval double expressionFrequency(Or<Expr...> op) noexcept
    {        
        return expressionFrequency(selectLast<ProjExprFrequency>(op.expr));
    }

    consteval double expressionFrequency(...) noexcept { return 0; }

    struct ProjExprFrequency
    {
        consteval auto operator()(auto expr)
        {
            return expressionFrequency(expr);
        }
    };

} // namespace Kitimar::CTSmarts
