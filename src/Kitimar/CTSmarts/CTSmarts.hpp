#pragma once

#include "Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    namespace detail {

        template<typename Molecule, typename Smarts, auto N>
        auto captureAtoms(Molecule &mol, Smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
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


        template<typename Molecule, typename Smarts, auto N>
        auto captureMatchAtoms(Molecule &mol, Smarts smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
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
    // Molecule
    //

    //
    // ctsmarts::match<SMARTS>(mol) -> bool
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr bool match(Molecule &mol)
    {
        auto iso = SingleIsomorphism<SMARTS>{};
        return iso.match(mol);
    }

    //
    // ctsmarts::single<SMARTS>(mol) -> std::vector<int>
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto single(Molecule &mol)
    {
        auto iso = SingleIsomorphism<SMARTS>{};
        return iso.single(mol);
    }

    //
    // ctsmarts::count<SMARTS>(mol) -> std::integeral
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto count(Molecule &mol)
    {
        auto iso = UniqueIsomorphism<SMARTS>{};
        return iso.count(mol);
    }

    //
    // ctsmarts::unique<SMARTS>(mol) -> std::vector<std::vector<int>>
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto unique(Molecule &mol)
    {
        auto iso = UniqueIsomorphism<SMARTS>{};
        return iso.all(mol);
    }

    //
    // ctsmarts::countAll<SMARTS>(mol) -> std::integeral
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto countAll(Molecule &mol)
    {
        auto iso = AllIsomorphism<SMARTS>{};
        return iso.count(mol);
    }

    //
    // ctsmarts::all<SMARTS>(mol) -> std::vector<std::vector<int>>
    //
    template<ctll::fixed_string SMARTS, typename Molecule>
    constexpr auto all(Molecule &mol)
    {
        auto iso = AllIsomorphism<SMARTS>{};
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
        auto iso = SingleIsomorphism<SMARTS>{};
        auto map = iso.single(mol);
        return detail::captureMatchAtoms(mol, smarts, map, cap);
    }

    template<MapType T>
    using CaptureTag = std::integral_constant<MapType, T>;

    static constexpr auto CaptureUnique = CaptureTag<MapType::Unique>{};
    static constexpr auto CaptureAll = CaptureTag<MapType::All>{};

    //
    // ctsmarts::captures<SMARTS>(mol) -> std::range<std::array<Atom...>>
    //
    template <ctll::fixed_string SMARTS, MapType CaptureType = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, CaptureTag<CaptureType> = {})
    {
        auto iso = Isomorphism<SMARTS, CaptureType>{};
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
    // Atom
    //

    //
    // ctsmarts::captures<SMARTS>(mol, atom) -> std::range<std::array<Atom...>>
    //
    template <ctll::fixed_string SMARTS, MapType CaptureType = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, const auto &atom, CaptureTag<CaptureType> = {})
    {
        auto iso = Isomorphism<SMARTS, CaptureType>{};
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
