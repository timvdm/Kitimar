#pragma once

#include <concepts>
#include <ranges>

namespace Kitimar::Molecule {



    template<typename R, typename AtomBond>
    concept AtomBondRange = std::ranges::input_range<R> &&
                        std::convertible_to<std::remove_cvref_t<std::ranges::range_value_t<R>>, AtomBond>;


    template<typename Mol>
    concept IsAtomList = requires (Mol &mol)
    {
        { num_atoms(mol) } -> std::convertible_to<std::size_t>;
        { get_atom(mol, 0) };
        { get_atoms(mol) } -> AtomBondRange<std::remove_cvref_t<decltype(get_atom(mol, 0))>>;
        { get_index(mol, get_atom(mol, 0)) } -> std::convertible_to<std::size_t>;
        { null_atom(mol) };
    };

    template<typename Mol>
    concept IsBondList = requires (Mol &mol)
    {
        { num_bonds(mol) } -> std::convertible_to<std::size_t>;
        { get_bond(mol, 0) };
        { get_bonds(mol) } -> AtomBondRange<std::remove_cvref_t<decltype(get_bond(mol, 0))>>;
        { get_index(mol, get_bond(mol, 0)) } -> std::convertible_to<std::size_t>;
        { get_source(mol, get_bond(mol, 0)) } -> std::convertible_to<decltype(get_atom(mol, 0))>;
        { get_target(mol, get_bond(mol, 0)) } -> std::convertible_to<decltype(get_atom(mol, 0))>;
        { null_bond(mol) };
    };

    template<typename Mol>
    concept MoleculeGraph = IsAtomList<Mol> &&
                            IsBondList<Mol>;

    template<typename Mol>
    concept IncidentBondList = requires (Mol &mol)
    {
        { get_degree(mol, get_atom(mol, 0)) } -> std::integral;
        { get_bonds(mol, get_atom(mol, 0)) } -> AtomBondRange<std::remove_cvref_t<decltype(get_bond(mol, 0))>>;
    };

    template<typename Mol>
    concept AdjacentAtomList = requires (Mol &mol)
    {
        { get_degree(mol, get_atom(mol, 0)) } -> std::integral;
        { get_nbrs(mol, get_atom(mol, 0)) } -> AtomBondRange<std::remove_cvref_t<decltype(get_atom(mol, 0))>>;
    };

    template<typename Mol>
    concept ElementLayer = requires (Mol &mol)
    {
        { get_element(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept IsotopeLayer = requires (Mol &mol)
    {
        { get_isotope(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept ChargeLayer = requires (Mol &mol)
    {
        { get_charge(mol, get_atom(mol, 0)) } -> std::signed_integral;
    };

    template<typename Mol>
    concept BondOrderLayer = requires (Mol &mol)
    {
        { get_order(mol, get_bond(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept ImplicitHydrogensLayer = requires (Mol &mol)
    {
        { get_implicit_hydrogens(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept MobileHydrogensLayer = requires (Mol &mol)
    {
        { get_mobile_hydrogens(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept AromaticLayer = requires (Mol &mol)
    {
        { is_aromatic_atom(mol, get_atom(mol, 0)) } -> std::same_as<bool>;
        { is_aromatic_bond(mol, get_bond(mol, 0)) } -> std::same_as<bool>;
    };



    template<typename Mol>
    concept Molecule = MoleculeGraph<Mol> &&
                       IncidentBondList<Mol> &&
                       AdjacentAtomList<Mol> &&
                       ElementLayer<Mol> &&
                       IsotopeLayer<Mol>;



    constexpr auto get_nbr(const auto &mol, auto bond, auto atom) noexcept -> decltype(get_atom(mol, 0))
    {
        auto source = get_source(mol, bond);
        return source != atom ? source : get_target(mol, bond);
    }

    constexpr auto get_bond(const auto &mol, auto source, auto target) noexcept -> decltype(get_bond(mol, 0))
    {
        auto sourceIndex1 = get_index(mol, source);
        auto targetIndex1 = get_index(mol, target);
        for (auto bond : get_bonds(mol, source)) {
            auto sourceIndex2 = get_index(mol, get_source(mol, bond));
            auto targetIndex2 = get_index(mol, get_target(mol, bond));
            if ((sourceIndex1 == sourceIndex2 && targetIndex1 == targetIndex2) ||
                (sourceIndex1 == targetIndex2 && targetIndex1 == sourceIndex2))
                return bond;
        }
        return null_bond(mol);
    }




} // namespace Kitimar::Molecule
