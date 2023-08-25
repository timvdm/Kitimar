#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <vector>
#include <algorithm>
#include <functional>
#include <type_traits>

namespace Kitimar::CTSmarts {

    namespace impl {

        /**
         * @brief Convert atom index map to atom map using capture set.
         */
        auto toCapture(Molecule::Molecule auto &mol, auto smarts,
                       const auto &captureSet, bool found, const auto &map)
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto NumCaptures = std::tuple_size_v<std::remove_cvref_t<decltype(captureSet)>>;
            if constexpr (NumCaptures) {
                std::array<Atom, NumCaptures> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < NumCaptures; ++i)
                        atoms[i] = get_atom(mol, map[captureSet[i]]);
                return atoms;
            } else {
                std::array<Atom, smarts.numAtoms> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < smarts.numAtoms; ++i)
                        atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
                return atoms;
            }
        }

        auto toCaptures(Molecule::Molecule auto &mol, const auto &iso,
                        const auto &captureSet, const auto &maps) noexcept
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto NumCaptures = std::tuple_size_v<std::remove_cvref_t<decltype(captureSet)>>;
            static constexpr auto M = NumCaptures ? NumCaptures : decltype(iso.smarts)::numAtoms;
            auto r = maps | std::views::transform([&] (const auto &map) {
                return toCapture(mol, iso.smarts, captureSet, true, map);
            });
            return std::vector<std::array<Atom, M>>{r.begin(), r.end()};
        }

        auto captureMatchAtoms(Molecule::Molecule auto &mol, auto smarts,
                               const auto &captureSet, bool found, const auto &map)
        {
            return std::tuple_cat(std::make_tuple(found), toCapture(mol, smarts, captureSet, found, map));
        }

        constexpr std::size_t captureHash(auto &mol, const auto &capture)
        {
            auto atoms = std::vector<bool>(num_atoms(mol));
            for (const auto &atom : capture)
                atoms[get_index(mol, atom)] = true;
            return std::hash<std::vector<bool>>()(atoms);
        }

        template<typename Proj = std::identity>
        constexpr auto numInversions(const auto &permutation, Proj proj = {})
        {
            std::size_t n = 0;
            for (auto i = 0; i < permutation.size(); ++i)
                for (auto j = i + 1; j < permutation.size(); ++j)
                    if (std::invoke(proj, permutation[i]) > std::invoke(proj, permutation[j]))
                        ++n;
            return n;
        }

        constexpr bool singleAtomMatch(auto smarts, auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr);
        }

        constexpr bool singleBondMatchHelper(auto smarts, auto &mol, const auto &bond, const auto &source, const auto &target)
        {
            return matchAtomExpr(mol, source, get<0>(smarts.atoms).expr) && matchAtomExpr(mol, target, get<1>(smarts.atoms).expr);
        }

        // 0 -> no match
        // 1 -> source is SMARTS atom 0, target is SMARTS atom 1
        // 2 -> source is SMARTS atom 1, target is SMARTS atom 0
        constexpr int singleBondMatch(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return 1;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return 2;
            return 0;
        }

        constexpr int singleBondCount(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            auto n = 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                ++n;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                ++n;
            return n;
        }

        template<typename Map>
        constexpr auto singleBondMap(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return std::make_tuple(false, Map{});

            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return std::make_tuple(true, Map{get_index(mol, get_source(mol, bond)),
                                                 get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return std::make_tuple(true, Map{get_index(mol, get_target(mol, bond)),
                                                 get_index(mol, get_source(mol, bond))});

            return std::make_tuple(false, Map{});
        }

        template<typename Map>
        constexpr void singleBondMaps(auto smarts, auto &mol, const auto &bond, std::vector<Map> &maps)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                maps.push_back({get_index(mol, get_source(mol, bond)),
                                get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                maps.push_back({get_index(mol, get_target(mol, bond)),
                                get_index(mol, get_source(mol, bond))});
        }

        constexpr auto singleBondCapture(auto smarts, auto &mol, const auto &bond, int singleBondMatchType)
        {
            constexpr auto cap = captureMapping(smarts);
            if constexpr (cap.size() == 1) {
                switch (singleBondMatchType) {
                    case 1:
                        return cap[0] == 0 ? std::make_tuple(true, get_source(mol, bond)) :
                                             std::make_tuple(true, get_target(mol, bond));
                    case 2:
                        return cap[0] == 0 ? std::make_tuple(true, get_target(mol, bond)) :
                                             std::make_tuple(true, get_source(mol, bond));
                    default:
                        return std::make_tuple(false, null_atom(mol));
                }
            } else {
                if (!singleBondMatchType)
                    return std::make_tuple(false, null_atom(mol), null_atom(mol));

                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                if constexpr (cap.size())
                    if (cap[0] > cap[1])
                        std::swap(source, target);

                if (singleBondMatchType == 1)
                    return std::make_tuple(true, source, target);
                return std::make_tuple(true, target, source);
            }
        }

        template<typename AtomMap>
        constexpr void singleBondCaptures(auto smarts, auto &mol, const auto &bond, const auto &captureSet, std::vector<AtomMap> &maps, bool unique)
        {
            using IndexMap = std::array<decltype(get_index(mol, get_atom(mol, 0))), 2>;
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;

            if (singleBondMatchHelper(smarts, mol, bond, source, target)) {
                auto map = IndexMap{get_index(mol, get_source(mol, bond)),
                                    get_index(mol, get_target(mol, bond))};
                maps.push_back(impl::toCapture(mol, smarts, captureSet, true, map));
                if (unique)
                    return;
            }

            if (singleBondMatchHelper(smarts, mol, bond, target, source)) {
                auto map = IndexMap{get_index(mol, get_target(mol, bond)),
                                    get_index(mol, get_source(mol, bond))};
                maps.push_back(impl::toCapture(mol, smarts, captureSet, true, map));
            }
        }

        constexpr auto centralAtomMap(auto smarts, auto &mol, const auto &atom)
        {
            std::array<int, smarts.numAtoms> map;
        }

    } // namespace impl


// function_search(mol) -> function(mol, Search)

#define CTSMARTS_API_SEARCH(function, Search, search) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig> \
    constexpr auto function##_##search(Molecule::Molecule auto &mol) \
    { return function<SMARTS, SearchType::Search, Config>(mol); }

#define CTSMARTS_API_UNIQUE(function) CTSMARTS_API_SEARCH(function, Unique, unique)
#define CTSMARTS_API_ALL(function) CTSMARTS_API_SEARCH(function, All, all)

// function_arg_search(mol, arg) -> function_arg(mol, arg, Search)

#define CTSMARTS_API_ARG_SEARCH(function, arg, Search, search) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig> \
    constexpr auto function##_##arg##_##search(Molecule::Molecule auto &mol, const auto &arg) \
    { return function##_##arg<SMARTS, SearchType::Search, Config>(mol, arg); }

#define CTSMARTS_API_ATOM_UNIQUE(function) CTSMARTS_API_ARG_SEARCH(function, atom, Unique, unique)
#define CTSMARTS_API_BOND_UNIQUE(function) CTSMARTS_API_ARG_SEARCH(function, bond, Unique, unique)
#define CTSMARTS_API_ATOM_ALL(function) CTSMARTS_API_ARG_SEARCH(function, atom, All, all)
#define CTSMARTS_API_BOND_ALL(function) CTSMARTS_API_ARG_SEARCH(function, bond, All, all)

// caller(mol, arg) -> callee(mol, arg)

#define CTSMARTS_API_OVERLOAD_ARG(caller, callee, Arg, arg) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol> \
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) \
    constexpr auto caller(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Arg arg) \
    { return callee<SMARTS, Config>(mol, arg); }

#define CTSMARTS_API_OVERLOAD_ATOM(function) CTSMARTS_API_OVERLOAD_ARG(function, function##_atom, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND(function) CTSMARTS_API_OVERLOAD_ARG(function, function##_bond, Bond, bond)

#define CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(function) CTSMARTS_API_OVERLOAD_ARG(function##_unique, function##_atom_##unique, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_UNIQUE(function) CTSMARTS_API_OVERLOAD_ARG(function##_unique, function##_bond_##unique, Bond, bond)
#define CTSMARTS_API_OVERLOAD_ATOM_ALL(function) CTSMARTS_API_OVERLOAD_ARG(function##_all, function##_atom_##all, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_ALL(function) CTSMARTS_API_OVERLOAD_ARG(function##_all, function##_bond_##all, Bond, bond)

// function(mol, arg, search) -> function_atom(mol, arg, search)

#define CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Arg, arg) \
    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol> \
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) \
    constexpr auto function(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Arg arg, SearchTypeTag<M> searchType = {}) \
    { return function##_##arg<SMARTS, M, Config>(mol, arg, searchType); }

#define CTSMARTS_API_OVERLOAD_ATOM_SEARCH(function) CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_SEARCH(function) CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Bond, bond)

} // mamespace Kitimar::CTSmarts
