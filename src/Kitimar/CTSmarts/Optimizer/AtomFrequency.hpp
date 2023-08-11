#pragma once

#include "ExpressionFrequency.hpp"

namespace Kitimar::CTSmarts {

    struct ProjAtomFrequency
    {        
        consteval auto operator()(auto atom)
        {
            return expressionFrequency(atom.expr);
        }
    };


    template<typename SmartsT>
    consteval auto makeAtomFrequency(auto atoms) noexcept
    {
        if constexpr (ctll::empty(atoms))
            return std::array<double, SmartsT::numAtoms>{};
        else {
            auto [atom, tail] = ctll::pop_and_get_front(atoms);
            auto freq = makeAtomFrequency<SmartsT>(tail);
            freq[atom.index] = expressionFrequency(atom.expr);
            return freq;
        }
    }

    template<typename SmartsT>
    struct AtomFrequency
    {
        static constexpr inline auto data = makeAtomFrequency<SmartsT>(SmartsT::atoms);

        consteval AtomFrequency() noexcept {}
        consteval AtomFrequency(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts
