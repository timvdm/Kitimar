#pragma once

#include "../../Molecule/Molecule.hpp"

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto enabledFilters(auto smarts, auto filters) noexcept
        {
            if constexpr (!ctll::empty(filters)) {
                auto [filter, tail] = ctll::pop_and_get_front(filters);
                if constexpr (filter.enable(smarts))
                    return ctll::push_front(filter, enabledFilters(smarts, tail));
                else
                    return enabledFilters(smarts, tail);
            } else
                return ctll::empty_list{};
        }

        template<typename ...Filter>
        constexpr bool rejectMolecule(auto smarts, Molecule::Molecule auto &mol, ctll::list<Filter...> filters) noexcept
        {
            return (Filter::reject(smarts, mol) || ...);
        }

    } // namespace impl


    template<typename SmartsT, typename ...Filter>
    struct FilterPolicyHelper
    {
        static constexpr inline auto filters = impl::enabledFilters(SmartsT{}, ctll::list<Filter...>{});

        static constexpr bool reject(Molecule::Molecule auto &mol) noexcept
        {
            return impl::rejectMolecule(SmartsT{}, mol, filters);
        }

        consteval FilterPolicyHelper() noexcept {}
        consteval FilterPolicyHelper(SmartsT, ctll::list<Filter...>) noexcept {}
    };


    template<typename ...Filter>
    struct FilterPolicy
    {
        static constexpr inline auto filters = ctll::list<Filter...>{};
    };



} // namespace Kitimar::CTSmarts
