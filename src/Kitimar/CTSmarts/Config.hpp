#pragma once

#include "Mapping.hpp"
#include "Optimizer/Optimizer.hpp"

namespace Kitimar::CTSmarts {

    enum class DefaultImplicitH
    {
        AtLeastOne, // Daylight, OpenSMARTS, RDKit
        ExactlyOne // OpenBabel
    };

    enum class Specialize : int
    {
        None,
        Atom  = 1,
        Bond  = 2,
        Chain = 4,
        Star  = 8,
        All = Atom | Bond | Chain | Star
    };

    inline constexpr bool operator&(Specialize lhs, Specialize rhs) noexcept
    {
        return static_cast<int>(lhs) & static_cast<int>(rhs);
    }

    namespace impl {

        //
        // replaceExpr
        //

        // Primitive
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, Expr) noexcept
        {
            if constexpr (std::is_same_v<Expr, OldExpr>)
                return NewExpr{};
            else
                return Expr{};
        }

        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, ctll::list<Expr...> l) noexcept;

        // Not
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, Not<Expr>) noexcept
        {
            return Not{replaceExpr(OldExpr{}, NewExpr{}, Expr{})};
        }

        // Or
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, Or<Expr...> op) noexcept
        {
            return Or{replaceExpr(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // And
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, And<Expr...> op) noexcept
        {
            return And{replaceExpr(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // Atom
        template<typename OldExpr, typename NewExpr, int Index, typename Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, Atom<Index, Expr>) noexcept
        {
            return Atom<Index, decltype(replaceExpr(OldExpr{}, NewExpr{}, Expr{}))>{};
        }

        // ctll::list
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExpr(OldExpr, NewExpr, ctll::list<Expr...> l) noexcept
        {
            return decltype(transform(l, [] (auto expr) { return replaceExpr(OldExpr{}, NewExpr{}, expr); })){};
        }

        //
        // containsExpr
        //

        // Primitive
        template<typename Query, typename Expr>
        constexpr auto containsExpr(Query, Expr) noexcept
        {
            return std::is_same_v<Query, Expr>;
        }

        template<typename Query, typename ...Expr> constexpr auto containsExpr(Query, ctll::list<Expr...>) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExpr(Query, Or<Expr...> op) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExpr(Query, And<Expr...> op) noexcept;

        // Not
        template<typename Query, typename Expr>
        constexpr auto containsExpr(Query, Not<Expr>) noexcept
        {
            return containsExpr(Query{}, Expr{});
        }

        // Or
        template<typename Query, typename ...Expr>
        constexpr auto containsExpr(Query, Or<Expr...> op) noexcept
        {
            return (containsExpr(Query{}, Expr{}) || ...);
        }

        // And
        template<typename Query, typename ...Expr>
        constexpr auto containsExpr(Query, And<Expr...> op) noexcept
        {
            return (containsExpr(Query{}, Expr{}) || ...);
        }

        // Atom
        template<typename Query, int Index, typename Expr>
        constexpr auto containsExpr(Query, Atom<Index, Expr>) noexcept
        {
            return containsExpr(Query{}, Expr{});
        }

        // ctll::list
        template<typename Query, typename ...Expr>
        constexpr auto containsExpr(Query, ctll::list<Expr...> l) noexcept
        {
            return (containsExpr(Query{}, Expr{}) || ...);
        }

    } // namespace impl



    template<DefaultImplicitH H, Specialize S, typename OptimizerT, template<std::integral, int N> class MapT>
    struct Config
    {
        // Matching behavior
        static constexpr inline auto defaultImplicitH = H;


        // Performance
        static constexpr inline auto specialize = S;

        using Optimizer = OptimizerT;

        template<std::integral Index, int N>
        using Map = MapT<Index, N>;

        static constexpr auto transformSmarts(auto smarts) noexcept
        {
            if constexpr (H == DefaultImplicitH::AtLeastOne)
                return smarts;
            else
                return BasicSmarts<decltype(impl::replaceExpr(HasImplicitH{}, ImplicitH<1>{}, smarts.atoms)), decltype(smarts.bonds), decltype(smarts.classes)>{};
        }
    };

    using DefaultConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::All, FullOptimizer, InverseMap>;

    using NoOptimizeConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::None, NoOptimizer, InverseMap>;




} // namespace Kitimar::CTSmarts
