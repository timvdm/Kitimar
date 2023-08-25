#pragma once

#include "Mapping.hpp"
#include "Optimizer/Optimizer.hpp"

namespace Kitimar::CTSmarts {

    template<bool Specialize, typename OptimizerT, template<std::integral, int N> class MapT>
    struct Config
    {
        static constexpr inline auto specialize = Specialize;

        using Optimizer = OptimizerT;

        template<std::integral Index, int N>
        using Map = MapT<Index, N>;
    };

    static constexpr inline auto Specialize = true;
    static constexpr inline auto NoSpecialize = false;

    using DefaultConfig = Config<Specialize, FullOptimizer, InverseMap>;
    //using DefaultConfig = Config<Specialize, FullOptimizer, LookupMap>;

    using NoOptimizeConfig = Config<NoSpecialize, NoOptimizer, InverseMap>;
    //using NoSpecializeConfig = Config<NoSpecialize, FullOptimizer, InverseMap>;

} // namespace Kitimar::CTSmarts
