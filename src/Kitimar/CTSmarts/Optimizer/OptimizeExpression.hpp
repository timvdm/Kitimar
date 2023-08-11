#pragma once

#include "ExpressionFrequency.hpp"

namespace Kitimar::CTSmarts {

    struct BestCase {};
    struct WorstCase {};

    // Not

    template<typename Goal = BestCase, typename Expr>
    consteval auto optimizeExpression(Not<Expr> op)
    {
        using NotGoal = std::conditional_t<std::is_same_v<Goal, BestCase>, WorstCase, BestCase>;
        return Not(optimizeExpression<NotGoal>(op.expr));
    }

    // And

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpression(And<Expr...> op)
    {        
        using Compare = std::conditional_t<std::is_same_v<Goal, BestCase>, std::less<>, std::greater<>>;
        return And(ctllSort<ProjExprFrequency, Compare>(optimizeExpressions<Goal>(op.expr)));
    }

    // Or

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpression(Or<Expr...> op)
    {
        using Compare = std::conditional_t<std::is_same_v<Goal, BestCase>, std::greater<>, std::less<>>;
        return Or(ctllSort<ProjExprFrequency, Compare>(optimizeExpressions<Goal>(op.expr)));
    }

    // Primitive

    template<typename Goal = BestCase>
    consteval auto optimizeExpression(auto expr)
    {
        return expr;
    }

    // ctll::list<Expr>

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpressions(ctll::list<Expr...> expressions) noexcept
    {
        if constexpr (!ctll::empty(expressions)) {
            auto [expr, tail] = ctll::pop_and_get_front(expressions);
            return ctll::push_front(optimizeExpression<Goal>(expr), optimizeExpressions<Goal>(tail));
        } else
            return ctll::empty_list{};
    }

} // namespace Kitimar::CTSmarts
