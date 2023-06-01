#pragma once

#include "Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    template<typename Molecule, typename Smarts, auto N>
    auto captureAtoms(Molecule &mol, Smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
    {
        using Atom = decltype(get_atom(mol, 0));
        if constexpr (N) {
            std::array<Atom, N> atoms = {};
            for (auto i = 0; i < N; ++i)
                atoms[i] = get_atom(mol, map[cap[i]]); // FIXME: null atoms...
            return atoms;
        } else {
            std::array<Atom, Smarts::numAtoms> atoms = {};
            for (auto i = 0; i < Smarts::numAtoms; ++i)
                atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
            return atoms;
        }
    }


    template<typename Molecule, typename Smarts, auto N>
    auto captureMatchAtoms(Molecule &mol, Smarts smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
    {
        return std::tuple_cat(std::make_tuple(!map.empty()), captureAtoms(mol, smarts, map, cap));
    }


    //
    // Molecule
    //


    //
    // ctsmarts::match<SMARTS>(mol) -> bool
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr bool match(Molecule &mol)
    {
        auto iso = SingleIsomorphism(Smarts<SMARTS>());
        return iso.match(mol);
    }

    //
    // ctsmarts::single<SMARTS>(mol) -> Mapping
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto single(Molecule &mol)
    {
        auto iso = SingleIsomorphism(Smarts<SMARTS>());
        return iso.single(mol);
    }

    //
    // ctsmarts::count<SMARTS>(mol) -> integer
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto count(Molecule &mol)
    {
        auto iso = UniqueIsomorphism(Smarts<SMARTS>());
        return iso.count(mol);
    }

    //
    // ctsmarts::unique<SMARTS>(mol) -> list<Mapping>
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto unique(Molecule &mol)
    {
        auto iso = UniqueIsomorphism(Smarts<SMARTS>());
        return iso.all(mol);
    }

    //
    // ctsmarts::countAll<SMARTS>(mol) -> integer
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto countAll(Molecule &mol)
    {
        auto iso = AllIsomorphism(Smarts<SMARTS>());
        return iso.count(mol);
    }

    //
    // ctsmarts::all<SMARTS>(mol) -> list<Mapping>
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto all(Molecule &mol)
    {
        auto iso = AllIsomorphism(Smarts<SMARTS>());
        return iso.all(mol);
    }

    //
    // ctsmarts::capture<SMARTS>(mol) -> tuple<bool, Atom...>
    //
    template <ctll::fixed_string SMARTS, typename Molecule>
    auto capture(Molecule &mol)
    {
        auto smarts = Smarts<SMARTS>();
        constexpr auto cap = captureMapping(smarts);
        SingleIsomorphism iso(smarts);
        auto map = iso.single(mol);

        return captureMatchAtoms(mol, smarts, map, cap);
    }

    //
    // ctsmarts::captures<SMARTS>(mol) -> std::range<std::tuple<bool, Atom...>>
    //
    template <ctll::fixed_string SMARTS, typename Molecule>
    auto captures(Molecule &mol)
    {
        auto smarts = Smarts<SMARTS>();
        constexpr auto cap = captureMapping(smarts);
        UniqueIsomorphism iso(smarts);
        return iso.all(mol) | std::views::transform([&mol, smarts, cap] (const auto &map) {
            return captureAtoms(mol, smarts, map, cap);
        });
    }


    //
    // Atom
    //


    //
    // ctsmarts::captures<SMARTS>(mol, atom) -> std::range<std::tuple<bool, Atom...>>
    //
    template <ctll::fixed_string SMARTS>
    auto captures(auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>();
        constexpr auto cap = captureMapping(smarts);
        UniqueIsomorphism iso(smarts);
        return iso.all(mol, atom) | std::views::transform([&mol, smarts, cap] (const auto &map) {
            return captureAtoms(mol, smarts, map, cap);
        });
    }


    //
    // API
    //

    /*


    Molecule mol;


    ctsmarts::match  <SMARTS> (mol)         -> boolean
    ctsmarts::count  <SMARTS> (mol, unique) -> integer
    ctsmarts::single <SMARTS> (mol)         -> mapping
    ctsmarts::unique <SMARTS> (mol)         -> list<mapping>
    ctsmarts::all    <SMARTS> (mol)         -> list<mapping>


    ctsmarts::match  <SMARTS> ([mol,] atom)          -> boolean
    ctsmarts::count  <SMARTS> ([mol,] atom, unique)  -> integer
    ctsmarts::single <SMARTS> ([mol,] atom)          -> mapping
    ctsmarts::unique <SMARTS> ([mol,] atom)          -> list<mapping>
    ctsmarts::all    <SMARTS> ([mol,] atom)          -> list<mapping>

    ctsmarts::capture  <SMARTS>(mol)                 -> tuple<boolean, atom...>
    ctsmarts::captures <SMARTS>(mol, unique)         -> list<tuple<boolean, atom...>>

    ctsmarts::capture  <SMARTS>([mol,] atom)         -> tuple<boolean, atom...>
    ctsmarts::captures <SMARTS>([mol,] atom, unique) -> list<tuple<boolean, atom...>>


    // FIXME captureAll captureAll


    auto [match, a, b] = ctsmarts<"[*:1]=[:2]">::capture(mol) -> tuple(boolean, atom...)

    auto [match, a, b] = ctsmarts<"[*:1]=[:2]">::capture(mol, atom) -> tuple(boolean, atom...)


    for (auto [a, b] : ctsmarts<"[*:1]=[:2]">::captures(mol, atom))


    */















} // namespace Kitimar::CTSmarts
