#pragma once

#include "AST/AST.hpp"
#include "Config.hpp"

#include "../Molecule/Molecule.hpp"

namespace Kitimar::CTSmarts {

    //
    // Operators
    //

    template<typename Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Not<Expr>)
    {
        return !matchAtomExpr(mol, atom, Expr{});
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Or<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr{}) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, And<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr{}) && ...);
    }

    template<typename Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Not<Expr>)
    {
        return !matchBondExpr(mol, bond, Expr{});
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Or<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr{}) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, And<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr{}) && ...);
    }

    //
    // Atoms
    //

    constexpr bool matchAtomExpr(const auto&, const auto&, AnyAtom)
    {
        return true;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAromatic)
    {
        return is_aromatic_atom(mol, atom);
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAliphatic)
    {
        return !is_aromatic_atom(mol, atom);
    }

    template<int Element>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AliphaticAtom<Element>)
    {
        return get_element(mol, atom) == Element && !is_aromatic_atom(mol, atom);
    }

    template<int Element>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AromaticAtom<Element>)
    {
        return get_element(mol, atom) == Element && is_aromatic_atom(mol, atom);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Isotope<N>)
    {
        return get_isotope(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Element<N>)
    {
        return get_element(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Degree<N>)
    {
        return get_degree(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Valence<N>)
    {
        return get_valence(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Connectivity<N>)
    {
        return get_degree(mol, atom) + get_implicit_hydrogens(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, TotalH<N>)
    {
        return get_total_hydrogens(mol, atom) == N;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, HasImplicitH)
    {
        return get_implicit_hydrogens(mol, atom) > 0;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, ImplicitH<N>)
    {
        return get_implicit_hydrogens(mol, atom) == N;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Cyclic)
    {
        return is_ring_atom(mol, atom);
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Acyclic)
    {
        return !is_ring_atom(mol, atom);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingCount<N>)
    {
        return get_ring_count(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingSize<N>)
    {
        return is_in_ring_size(mol, atom, N);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingConnectivity<N>)
    {
        return get_ring_degree(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Charge<N>)
    {
        return get_charge(mol, atom) == N;
    }

    //
    // Recursive
    //

    namespace impl {
        template<typename SmartsT, typename Config, Molecule::Molecule Mol>
        constexpr bool match_atom(const Mol &mol, const auto &atom);
    }

    template<Molecule::Molecule Mol, typename Atoms, typename Bonds, typename Classes>
    constexpr bool matchAtomExpr(const Mol &mol, const auto &atom, BasicSmarts<Atoms, Bonds, Classes>)
    {
        return impl::match_atom<BasicSmarts<Atoms, Bonds, Classes>, DefaultConfig, Mol>(mol, atom);
    }

    //
    // Bonds
    //

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, ImplicitBond)
    {
        return get_order(mol, bond) == 1 || is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto&, const auto&, AnyBond)
    {
        return true;
    }

    template<int Order>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, BondOrder<Order>)
    {
        return get_order(mol, bond) == Order && !is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, AromaticBond)
    {
        return is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, RingBond)
    {
        return is_ring_bond(mol, bond);
    }

} // namespace ctsmarts
