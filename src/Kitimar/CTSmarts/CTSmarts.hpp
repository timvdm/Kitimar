#pragma once

#include "Isomorphism.hpp"
#include "../Molecule/MockMolecule.hpp" // for testing single header using compiler explorer

namespace Kitimar::CTSmarts {

    namespace detail {

        template<typename Smarts, auto N>
        auto captureAtoms(Molecule::Molecule auto &mol, Smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
        {
            using Atom = decltype(get_atom(mol, 0));
            if constexpr (N) {
                std::array<Atom, N> atoms = {};
                std::cout << "atoms.size(): " << atoms.size() << std::endl;
                std::cout << "map.size(): " << map.size() << std::endl;
                std::cout << "cap.size(): " << cap.size() << std::endl;
                if (map.empty())
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < N; ++i) {
                        std::cout << "i: " << i << std::endl;
                        std::cout << "    cap[i]: " << cap[i] << std::endl;
                        std::cout << "    map[cap[i]]: " << map[cap[i]] << std::endl;
                        atoms[i] = get_atom(mol, map[cap[i]]);
                    }
                return atoms;
            } else {
                std::array<Atom, Smarts::numAtoms> atoms = {};
                if (map.empty())
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < Smarts::numAtoms; ++i)
                        atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
                return atoms;
            }
        }


        template<typename Smarts, auto N>
        auto captureMatchAtoms(Molecule::Molecule auto &mol, Smarts smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
        {
            return std::tuple_cat(std::make_tuple(!map.empty()), captureAtoms(mol, smarts, map, cap));
        }

        template<auto N>
        auto copyCapture(Molecule::Molecule auto &mol, const auto &iso, const std::array<int, N> &cap, const auto &caps) noexcept
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto M = N ? N : iso.smarts.numAtoms;
            std::vector<std::array<Atom, M>> v;
            auto r = caps | std::views::transform([&] (const auto &map) {
                return captureAtoms(mol, iso.smarts, map, cap);
            });
            std::ranges::copy(r, std::back_inserter(v));
            return v;
        }

    } // namespace detail


    //
    // CTSmarts::contains<"SMARTS">(mol) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool contains(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        return iso.match(mol);
    }

    //
    // CTSmarts::atom<"SMARTS">(mol, atom) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool atom(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.numAtoms == 1) {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms));
        } else {
            auto smarts = Smarts<SMARTS>{};
            auto iso = Isomorphism{smarts, Single};
            return iso.match(mol, atom);
        }
    }

    //
    // CTSmarts::bond<"SMARTS">(mol, bond) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool bond(Molecule::Molecule auto &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        auto source = get_source(mol, bond);
        auto target = get_target(mol, bond);
        if constexpr (smarts.numAtoms == 2 && smarts.numBonds == 1) {
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return false;
            if (matchAtomExpr(mol, source, get<0>(smarts.atoms)))
                return matchAtomExpr(mol, target, get<1>(smarts.atoms));
            return matchAtomExpr(mol, source, get<1>(smarts.atoms)) &&
                   matchAtomExpr(mol, target, get<0>(smarts.atoms));
        } else {
            auto smarts = Smarts<SMARTS>{};
            auto iso = Isomorphism{smarts, Single};
            if (iso.match(mol, source))
                return true;
            return iso.match(mol, target);
        }
    }

    //
    // CTSmarts::count<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::integeral
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    constexpr auto count(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::contains()");
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        return iso.count(mol);
    }

    //
    // CTSmarts::single<"SMARTS">(mol) -> std::vector<int>
    //

    template<ctll::fixed_string SMARTS>
    constexpr auto single(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        return iso.single(mol);
    }

    //
    // CTSmarts::single<"SMARTS">(mol, atom) -> std::vector<int>
    //

    template<ctll::fixed_string SMARTS>
    constexpr auto single(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        return iso.single(mol, atom);
    }

    //
    // CTSmarts::multi<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::vector<int>>
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    constexpr auto multi(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        return iso.all(mol);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS>
    auto capture(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(iso.smarts);
        auto map = iso.single(mol);
        return detail::captureMatchAtoms(mol, iso.smarts, map, cap);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS>
    auto capture(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(iso.smarts);
        auto map = iso.single(mol, atom);
        return detail::captureMatchAtoms(mol, iso.smarts, map, cap);
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        static constexpr auto cap = captureMapping(iso.smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, iso.smarts, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol));
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        static constexpr auto cap = captureMapping(iso.smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol, atom) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, iso.smarts, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol, atom));
    }

} // namespace Kitimar::CTSmarts
