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
                return BasicSmarts<decltype(replaceExpr(HasImplicitH{}, ImplicitH<1>{}, smarts.atoms)), decltype(smarts.bonds), decltype(smarts.classes)>{};
        }
    };

    using DefaultConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::All, FullOptimizer, InverseMap>;

    using NoOptimizeConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::None, NoOptimizer, InverseMap>;




} // namespace Kitimar::CTSmarts
