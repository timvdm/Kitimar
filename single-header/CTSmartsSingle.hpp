#pragma once

#include <concepts>
#include <ranges>

namespace Kitimar::Molecule {

    /*
     *
     * AtomSet
     *
     *
     *     
     *
     *
     *
     *

     *
     * CyclicMolecule
     * - is_cylic / is_aclyclic
     *
     *
     *
     */

    //template<typename Mol>
    //using atom_t = std::remove_cvref_t<decltype(get_atom)

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

    constexpr auto get_nbr(const auto &mol, auto bond, auto atom) noexcept
    {
        auto source = get_source(mol, bond);
        return source != atom ? source : get_target(mol, bond);
    }

} // namespace Kitimar::Molecule

#include <ctll/list.hpp>

#include <tuple>
#include <stdexcept>

namespace Kitimar::CTSmarts {

    inline constexpr void PRE(auto cond)
    {
        if (!cond) throw std::runtime_error("");
    }

    template<int Value>
    struct Number
    {
        static constexpr inline auto value = Value;
    };

    template<int Value>
    struct Char
    {
        static constexpr inline auto value = Value;
    };

    //
    // ctll::list extensions
    //

    template<int Index, typename ...Ts>
    constexpr auto get(ctll::list<Ts...>)
    {
        return std::get<Index>(std::tuple<Ts...>());
    }

    template<int Size, typename T, typename ...Ts>
    constexpr auto resize(ctll::list<Ts...> = {})
    {
        if constexpr (ctll::size(ctll::list<Ts...>()) == Size)
            return ctll::list<Ts...>();
        else
            return resize<Size, T>(ctll::list<T, Ts...>());
    }

    template<int Index, typename U, typename T, typename ...Ts, typename ...Us>
    constexpr auto replace(ctll::list<T, Ts...>, ctll::list<Us...> = {})
    {
        if constexpr (Index)
            return replace<Index - 1, U>(ctll::list<Ts...>(), ctll::list<Us..., T>());
        else
            return ctll::list<Us..., U, Ts...>();
    }

    template<typename Separator, typename ...Ts, typename ...Us, typename ...Vs>
    constexpr auto split(ctll::list<Ts...> list, ctll::list<Us...> parts = ctll::empty_list(), ctll::list<Vs...> part = ctll::empty_list())
    {
        if constexpr (ctll::empty(list)) {
            static_assert(!ctll::empty(part));
            return ctll::list<Us..., ctll::list<Vs...>>();
        } else {
            auto [head, tail] = ctll::pop_and_get_front(list);
            if constexpr (std::is_same_v<decltype(head), Separator>) {
                static_assert(!ctll::empty(part));
                auto parts2 = ctll::list<Us..., ctll::list<Vs...>>();
                return split<Separator>(tail, parts2, ctll::empty_list());
            } else {
                auto part2 = ctll::list<Vs..., decltype(head)>();
                return split<Separator>(tail, parts, part2);
            }
        }
    }

    template<int ...N>
    constexpr auto to_array(ctll::list<Number<N>...>)
    {
        return std::array<int, sizeof...(N)>({ N... });
    }

    //
    // Helper functions to create runtime variable from compile time type in ctll::list
    //

    namespace detail {

        template<std::size_t I, typename R, typename F>
        inline constexpr R with_integral_constant(F f)
        {
            return static_cast<F>(f)(std::integral_constant<std::size_t, I>{});
        }

    } // namespace detail

    template<std::size_t N, typename R = void, typename F>
    inline constexpr R with_n(int n, F &&f)
    {
        constexpr auto invoke_array = [] <std::size_t...I> (std::index_sequence<I...>) {
            return std::array{ detail::with_integral_constant<I, R, F&&>... };
        }(std::make_index_sequence<N>{});

        return invoke_array.at(n)(std::forward<F>(f));
    }

} // namespace ctsmarts

#include <ctll/list.hpp>
#include <ctll/grammars.hpp>

namespace Kitimar::CTSmarts {

    struct SmartsGrammar
    {
        //
        // Symbols
        //

        struct atom {};
        struct atom_B {};
        struct atom_C {};
        struct atom_expr {}; // FIXME: rename to bracket_atom
        struct atom_expr2 {}; // FIXME: rename to atom_expr
        struct atom_exprA {}; // atom_expr_A
        struct atom_exprB {};
        struct atom_exprC {};
        struct atom_exprD {};
        struct atom_exprE {};
        struct atom_exprF {};
        struct atom_exprG {};
        struct atom_exprH {};
        struct atom_exprI {};
        struct atom_exprK {};
        struct atom_exprL {};
        struct atom_exprM {};
        struct atom_exprN {};
        struct atom_exprO {};
        struct atom_exprP {};
        struct atom_exprR {};
        struct atom_exprS {};
        struct atom_exprT {};
        struct atom_exprX {};
        struct atom_exprY {};
        struct atom_exprZ {};
        struct atom_expr_a {};
        struct atom_expr_s {};
        struct atom_expr_isotope {};
        struct atom_expr_isotope2 {};
        struct atom_expr_element {};
        struct atom_expr_element2 {};
        struct atom_expr_degree {};
        struct atom_expr_valence {};
        struct atom_expr_valence2 {};
        struct atom_expr_connectivity {};
        struct atom_expr_total_h {};
        struct atom_expr_impl_h {};
        struct atom_expr_ring_count {};
        struct atom_expr_ring_size {};
        struct atom_expr_ring_size2 {};
        struct atom_expr_ring_connectivity {};
        struct atom_expr_ring_connectivity2 {};
        struct atom_expr_neg_charge {};
        struct atom_expr_pos_charge {};
        struct atom_expr_chiral {};
        struct atom_expr_class {};
        struct atom_expr_class2 {};
        struct bond_expr {};
        struct bond_expr2 {};
        struct chain {};
        struct chain_up_down {};
        struct ring_bond {};
        struct ring_bond2 {};

        //
        // Actions
        //

        // Create atom AST elements
        struct make_any_atom : ctll::action {}; // '*'
        struct make_any_aliphatic : ctll::action {}; // 'A'
        struct make_any_aromatic : ctll::action {}; // 'a'
        struct make_aliphatic : ctll::action {}; // single letter symbols
        struct make_aromatic : ctll::action {}; // single letter symbols

        struct make_isotope : ctll::action {};
        struct make_element : ctll::action {};
        struct make_degree : ctll::action {};
        struct make_valence : ctll::action {};
        struct make_connectivity : ctll::action {};
        struct make_total_h : ctll::action {};
        struct make_impl_h : ctll::action {};
        struct make_cyclic : ctll::action {};
        struct make_acyclic : ctll::action {};
        struct make_ring_count : ctll::action {};
        struct make_ring_size : ctll::action {};
        struct make_ring_connectivity : ctll::action {};
        struct make_neg_charge : ctll::action {};
        struct make_pos_charge : ctll::action {};
        struct make_chiral : ctll::action {};
        struct make_class : ctll::action {};

        struct make_bond_primitive : ctll::action {};
        //struct make_up_or_down_bond : ctll::action {};

        // Create bond AST elements
        struct next_atom : ctll::action {};

        // Branches
        struct push_prev : ctll::action {};
        struct pop_prev : ctll::action {};
        struct reset_prev : ctll::action {};
        struct set_bond_type : ctll::action {};
        struct handle_ring_bond : ctll::action {};

        // Operators
        /*
        struct toggle_not : ctll::action {};
        struct reset_not : ctll::action {};
        struct set_and_high : ctll::action {};
        struct set_or : ctll::action {};
        struct set_and_low : ctll::action {};
        struct set_no_op : ctll::action {};
        */

        struct make_atom_not : ctll::action {};
        struct make_atom_and_high : ctll::action {};
        struct make_atom_or : ctll::action {};
        struct make_atom_and_low : ctll::action {};

        struct make_bond_not : ctll::action {};
        struct make_bond_and_high : ctll::action {};
        struct make_bond_or : ctll::action {};
        struct make_bond_and_low : ctll::action {};

        // Helpers
        struct push_char : ctll::action {};
        struct pop_char : ctll::action {};

        struct start_number : ctll::action {};
        struct push_number : ctll::action {};

        struct error_empty_bracket : ctll::action {}; // '[]'

        using _start = atom;

        //
        // Organic atoms
        //

        // aliphatic atom: 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
        static constexpr auto rule(atom, ctll::set<'N', 'O', 'P', 'S', 'F', 'I'>) -> ctll::push<ctll::anything, make_aliphatic, next_atom, chain>;
        static constexpr auto rule(atom, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_B>;
        static constexpr auto rule(atom, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_C>;
        static constexpr auto rule(atom_B, ctll::term<'r'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain>;
        static constexpr auto rule(atom_B, ctll::neg_set<'r'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain>;
        static constexpr auto rule(atom_B, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        static constexpr auto rule(atom_C, ctll::term<'l'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain>;
        static constexpr auto rule(atom_C, ctll::neg_set<'l'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain>;
        static constexpr auto rule(atom_C, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        // aromatic atom: 'b' | 'c' | 'n' | 'o' | 's' | 'p'
        static constexpr auto rule(atom, ctll::set<'b', 'c','n','o','p','s'>) -> ctll::push<ctll::anything, make_aromatic, next_atom, chain>;
        // any atom: '*'
        static constexpr auto rule(atom, ctll::term<'*'>) -> ctll::push<ctll::anything, make_any_atom, next_atom, chain>;
        // any aromatic: 'a'
        static constexpr auto rule(atom, ctll::term<'a'>) -> ctll::push<ctll::anything, make_any_aromatic, next_atom, chain>;
        // any aliphatic: 'A'
        static constexpr auto rule(atom, ctll::term<'A'>) -> ctll::push<ctll::anything, make_any_aliphatic, next_atom, chain>;
        // bracket atom: '[' atom_expression+ ']'
        static constexpr auto rule(atom, ctll::term<'['>) -> ctll::push<ctll::anything, atom_expr>;

        //
        // Atom expressions (i.e. [ ... ] )
        //

        using digit_chars = ctll::range<'0', '9'>;
        using not_digit_chars = ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;

        // check for [], [!] or [...!]
        static constexpr auto rule(atom_expr, ctll::set<']'>) -> ctll::push<error_empty_bracket, ctll::reject>;
        static constexpr auto rule(atom_expr, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr>;
        static constexpr auto rule(atom_expr, ctll::neg_set<']', '!'>) -> ctll::push<atom_expr2>;

        // operations + end ]
        static constexpr auto rule(atom_expr2, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr>;
        static constexpr auto rule(atom_expr2, ctll::term<'&'>) -> ctll::push<ctll::anything, make_atom_and_high, atom_expr>;
        static constexpr auto rule(atom_expr2, ctll::term<','>) -> ctll::push<ctll::anything, make_atom_or, atom_expr>;
        static constexpr auto rule(atom_expr2, ctll::term<';'>) -> ctll::push<ctll::anything, make_atom_and_low, atom_expr>;
        static constexpr auto rule(atom_expr2, ctll::term<']'>) -> ctll::push<ctll::anything, next_atom, chain>;

        // isotope: NUMBER
        static constexpr auto rule(atom_expr2, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_isotope>;
        static constexpr auto rule(atom_expr_isotope, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_isotope>;
        static constexpr auto rule(atom_expr_isotope, not_digit_chars) -> ctll::push<make_isotope, atom_expr2>;

        // element: '#' | '#' NUMBER
        static constexpr auto rule(atom_expr2, ctll::term<'#'>) -> ctll::push<ctll::anything, atom_expr_element>;
        static constexpr auto rule(atom_expr_element, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_element2>;
        static constexpr auto rule(atom_expr_element, not_digit_chars) -> ctll::push<make_element, atom_expr2>;
        static constexpr auto rule(atom_expr_element2, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_element2>;
        static constexpr auto rule(atom_expr_element2, not_digit_chars) -> ctll::push<make_element, atom_expr2>;

        // valence: 'v' | 'v' NUMBER
        static constexpr auto rule(atom_expr2, ctll::term<'v'>) -> ctll::push<ctll::anything, atom_expr_valence>;
        static constexpr auto rule(atom_expr_valence, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_valence2>;
        static constexpr auto rule(atom_expr_valence, not_digit_chars) -> ctll::push<make_valence, atom_expr2>;
        static constexpr auto rule(atom_expr_valence2, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_valence2>;
        static constexpr auto rule(atom_expr_valence2, not_digit_chars) -> ctll::push<make_valence, atom_expr2>;

        // implicit hydrogens: 'h' | 'h' DIGIT
        static constexpr auto rule(atom_expr2, ctll::term<'h'>) -> ctll::push<ctll::anything, atom_expr_impl_h>;
        static constexpr auto rule(atom_expr_impl_h, digit_chars) -> ctll::push<ctll::anything, make_impl_h, atom_expr2>;
        static constexpr auto rule(atom_expr_impl_h, not_digit_chars) -> ctll::push<make_impl_h, atom_expr2>;

        // cyclic: 'r'
        // acyclic: 'r0'
        // ring size: 'r' NUMBER
        // FIXME: r1 & r2  = error??
        static constexpr auto rule(atom_expr2, ctll::term<'r'>) -> ctll::push<ctll::anything, atom_expr_ring_size>;
        static constexpr auto rule(atom_expr_ring_size, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2>;
        static constexpr auto rule(atom_expr_ring_size, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_size2>;
        static constexpr auto rule(atom_expr_ring_size, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2>;
        static constexpr auto rule(atom_expr_ring_size2, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_size2>;
        static constexpr auto rule(atom_expr_ring_size2, not_digit_chars) -> ctll::push<make_ring_size, atom_expr2>;

        // cyclic: 'x'
        // acyclic: 'x0'
        // ring connectivity: 'x' NUMBER
        // FIXME: x1 = error?
        static constexpr auto rule(atom_expr2, ctll::term<'x'>) -> ctll::push<ctll::anything, atom_expr_ring_connectivity>;
        static constexpr auto rule(atom_expr_ring_connectivity, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2>;
        static constexpr auto rule(atom_expr_ring_connectivity, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_connectivity2>;
        static constexpr auto rule(atom_expr_ring_connectivity, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2>;
        static constexpr auto rule(atom_expr_ring_connectivity2, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_connectivity2>;
        static constexpr auto rule(atom_expr_ring_connectivity2, not_digit_chars) -> ctll::push<make_ring_connectivity, atom_expr2>;

        // charge: '-' | '-' DIGIT | '+' | '+' DIGIT
        static constexpr auto rule(atom_expr2, ctll::term<'-'>) -> ctll::push<ctll::anything, atom_expr_neg_charge>;
        static constexpr auto rule(atom_expr2, ctll::term<'+'>) -> ctll::push<ctll::anything, atom_expr_pos_charge>;
        static constexpr auto rule(atom_expr_neg_charge, digit_chars) -> ctll::push<ctll::anything, make_neg_charge, atom_expr2>;
        static constexpr auto rule(atom_expr_neg_charge, not_digit_chars) -> ctll::push<make_neg_charge, atom_expr2>;
        static constexpr auto rule(atom_expr_pos_charge, digit_chars) -> ctll::push<ctll::anything, make_pos_charge, atom_expr2>;
        static constexpr auto rule(atom_expr_pos_charge, not_digit_chars) -> ctll::push<make_pos_charge, atom_expr2>;

        // atom class: ':' NUMBER
        static constexpr auto rule(atom_expr2, ctll::term<':'>) -> ctll::push<ctll::anything, atom_expr_class>;
        static constexpr auto rule(atom_expr_class, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_class2>;
        static constexpr auto rule(atom_expr_class, not_digit_chars) -> ctll::reject; //ctll::push<make_class, set_and_high, atom_expr2>;
        static constexpr auto rule(atom_expr_class2, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_class2>;
        static constexpr auto rule(atom_expr_class2, not_digit_chars) -> ctll::push<make_class, atom_expr2>;

        // symbol: 'U' | 'V' | 'W'
        static constexpr auto rule(atom_expr2, ctll::set<'U', 'V', 'W'>) -> ctll::push<ctll::anything, make_aliphatic, atom_expr2>;
        // symbol: 'b' | 'c' | 'n' | 'o' | 'p'
        static constexpr auto rule(atom_expr2, ctll::set<'b', 'c', 'n', 'o', 'p'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2>;
        // any atom: '*'
        static constexpr auto rule(atom_expr2, ctll::set<'*'>) -> ctll::push<ctll::anything, make_any_atom, atom_expr2>;

        static constexpr auto rule(atom_expr2, ctll::set<'A'>) -> ctll::push<ctll::anything, push_char, atom_exprA>;
        static constexpr auto rule(atom_expr2, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_exprB>;
        static constexpr auto rule(atom_expr2, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_exprC>;
        static constexpr auto rule(atom_expr2, ctll::set<'D'>) -> ctll::push<ctll::anything, push_char, atom_exprD>;
        static constexpr auto rule(atom_expr2, ctll::set<'E'>) -> ctll::push<ctll::anything, push_char, atom_exprE>;
        static constexpr auto rule(atom_expr2, ctll::set<'F'>) -> ctll::push<ctll::anything, push_char, atom_exprF>;
        static constexpr auto rule(atom_expr2, ctll::set<'G'>) -> ctll::push<ctll::anything, push_char, atom_exprG>;
        static constexpr auto rule(atom_expr2, ctll::set<'H'>) -> ctll::push<ctll::anything, push_char, atom_exprH>;
        static constexpr auto rule(atom_expr2, ctll::set<'I'>) -> ctll::push<ctll::anything, push_char, atom_exprI>;
        static constexpr auto rule(atom_expr2, ctll::set<'K'>) -> ctll::push<ctll::anything, push_char, atom_exprK>;
        static constexpr auto rule(atom_expr2, ctll::set<'L'>) -> ctll::push<ctll::anything, push_char, atom_exprL>;
        static constexpr auto rule(atom_expr2, ctll::set<'M'>) -> ctll::push<ctll::anything, push_char, atom_exprM>;
        static constexpr auto rule(atom_expr2, ctll::set<'N'>) -> ctll::push<ctll::anything, push_char, atom_exprN>;
        static constexpr auto rule(atom_expr2, ctll::set<'O'>) -> ctll::push<ctll::anything, push_char, atom_exprO>;
        static constexpr auto rule(atom_expr2, ctll::set<'P'>) -> ctll::push<ctll::anything, push_char, atom_exprP>;
        static constexpr auto rule(atom_expr2, ctll::set<'R'>) -> ctll::push<ctll::anything, push_char, atom_exprR>;
        static constexpr auto rule(atom_expr2, ctll::set<'S'>) -> ctll::push<ctll::anything, push_char, atom_exprS>;
        static constexpr auto rule(atom_expr2, ctll::set<'T'>) -> ctll::push<ctll::anything, push_char, atom_exprT>;
        static constexpr auto rule(atom_expr2, ctll::set<'X'>) -> ctll::push<ctll::anything, push_char, atom_exprX>;
        static constexpr auto rule(atom_expr2, ctll::set<'Y'>) -> ctll::push<ctll::anything, push_char, atom_exprY>;
        static constexpr auto rule(atom_expr2, ctll::set<'Z'>) -> ctll::push<ctll::anything, push_char, atom_exprZ>;
        static constexpr auto rule(atom_expr2, ctll::set<'a'>) -> ctll::push<ctll::anything, push_char, atom_expr_a>;
        static constexpr auto rule(atom_expr2, ctll::set<'s'>) -> ctll::push<ctll::anything, push_char, atom_expr_s>;

        using Symbol1 = ctll::push<pop_char, make_aliphatic, atom_expr2>;
        using Symbol2 = ctll::push<ctll::anything, make_aliphatic, atom_expr2>;

        // symbol: 'Ac' | 'Ag' | 'Al' | 'Am' | 'Ar' | 'As' | 'At' | 'Au'
        // any aliphatic: 'A'
        static constexpr auto rule(atom_exprA,     ctll::set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprA, ctll::neg_set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> ctll::push<pop_char, make_any_aliphatic, atom_expr2>;
        // symbol: 'B' | 'Ba' | 'Be' | 'Bh' | 'Bi' | 'Bk' | 'Br'
        static constexpr auto rule(atom_exprB,     ctll::set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprB, ctll::neg_set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol1;
        // C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu
        static constexpr auto rule(atom_exprC,     ctll::set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprC, ctll::neg_set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol1;
        // symbol: 'Db' | 'Ds' | 'Dy'
        // degree: 'D' | 'D' NUMBER
        static constexpr auto rule(atom_exprD,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_degree>;
        static constexpr auto rule(atom_exprD,     ctll::set<'b', 's', 'y'>) -> Symbol2;
        static constexpr auto rule(atom_exprD, ctll::neg_set<'b', 's', 'y', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_degree, atom_expr2>;
        static constexpr auto rule(atom_expr_degree,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_degree>;
        static constexpr auto rule(atom_expr_degree, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_degree, atom_expr2>;
        // symbol: 'Er' | 'Es' | 'Eu'
        static constexpr auto rule(atom_exprE,     ctll::set<'r', 's', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprE, ctll::neg_set<'r', 's', 'u'>) -> ctll::reject;
        // symbol: 'F' | 'Fe' | 'Fm' | 'Fr'
        static constexpr auto rule(atom_exprF,     ctll::set<'e', 'm', 'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprF, ctll::neg_set<'e', 'm', 'r'>) -> Symbol1;
        // symbol: 'Ga' | 'Gd' | 'Ge'
        static constexpr auto rule(atom_exprG,     ctll::set<'a', 'd', 'e'>) -> Symbol2;
        static constexpr auto rule(atom_exprG, ctll::neg_set<'a', 'd', 'e'>) -> ctll::reject;
        // symbol: 'H' | 'He' | 'Hf' | 'Hg' | 'Ho' | 'Hs'
        // total hydrogens: 'H' | 'H' DIGIT
        static constexpr auto rule(atom_exprH,     ctll::set<'e', 'f', 'g', 'o', 's'>) -> Symbol2;
        static constexpr auto rule(atom_exprH,   ctll::range<'0', '9'>) -> ctll::list<ctll::anything, pop_char, make_total_h, atom_expr2>;
        static constexpr auto rule(atom_exprH, ctll::neg_set<'e', 'f', 'g', 'o', 's', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::list<pop_char, make_total_h, atom_expr2>;
        // symbol: 'I' | 'In' | 'Ir'
        static constexpr auto rule(atom_exprI,     ctll::set<'n', 'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprI, ctll::neg_set<'n', 'r'>) -> Symbol1;
        // symbol: 'K' | 'Kr'
        static constexpr auto rule(atom_exprK,    ctll::term<'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprK, ctll::neg_set<'r'>) -> Symbol1;
        // La Li Lr Lu
        static constexpr auto rule(atom_exprL,     ctll::set<'a', 'i', 'r', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprL, ctll::neg_set<'a', 'i', 'r', 'u'>) -> ctll::reject;
        // Md Mg Mn Mo Mt
        static constexpr auto rule(atom_exprM,     ctll::set<'d', 'g', 'n', 'o', 't'>) -> Symbol2;
        static constexpr auto rule(atom_exprM, ctll::neg_set<'d', 'g', 'n', 'o', 't'>) -> ctll::reject;
        // N Na Nb Nd Ne Ni No Np
        static constexpr auto rule(atom_exprN,     ctll::set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol2;
        static constexpr auto rule(atom_exprN, ctll::neg_set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol1;
        // O Os
        static constexpr auto rule(atom_exprO,    ctll::term<'s'>) -> Symbol2;
        static constexpr auto rule(atom_exprO, ctll::neg_set<'s'>) -> Symbol1;
        // P Pa Pb Pd Pm Po Pr Pt Pu
        static constexpr auto rule(atom_exprP,     ctll::set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprP, ctll::neg_set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol1;
        // symbol: 'Ra' | 'Rb' | 'Re' | 'Rf' | 'Rg' | 'Rh' | 'Rn' | 'Ru'
        // cyclic: 'R'
        // acyclic: 'R0'
        // ring count: 'R' NUMBER
        static constexpr auto rule(atom_exprR,   ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2>;
        static constexpr auto rule(atom_exprR,   ctll::range<'1', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_ring_count>;
        static constexpr auto rule(atom_exprR,     ctll::set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u'>) -> Symbol2;
        static constexpr auto rule(atom_exprR, ctll::neg_set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_cyclic, atom_expr2>;
        static constexpr auto rule(atom_expr_ring_count,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_ring_count>;
        static constexpr auto rule(atom_expr_ring_count, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_ring_count, atom_expr2>;
        // S Sb Sc Se Sg Si Sm Sn Sr
        static constexpr auto rule(atom_exprS,     ctll::set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprS, ctll::neg_set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol1;
        // Ta Tb Tc Te Th Ti Tl Tm
        static constexpr auto rule(atom_exprT,     ctll::set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> Symbol2;
        static constexpr auto rule(atom_exprT, ctll::neg_set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> ctll::reject;
        // symbol: 'Xe'
        // connectivity: 'X' | 'X' NUMBER
        static constexpr auto rule(atom_exprX, digit_chars) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_connectivity>;
        static constexpr auto rule(atom_exprX,    ctll::term<'e'>) -> Symbol2;
        static constexpr auto rule(atom_exprX, ctll::neg_set<'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_connectivity, atom_expr2>;
        static constexpr auto rule(atom_expr_connectivity, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_connectivity>;
        static constexpr auto rule(atom_expr_connectivity, not_digit_chars) -> ctll::push<make_connectivity, atom_expr2>;
        // symbol: 'Y' | 'Yb'
        static constexpr auto rule(atom_exprY,    ctll::term<'b'>) -> Symbol2;
        static constexpr auto rule(atom_exprY, ctll::neg_set<'b'>) -> Symbol1;
        // symbol: 'Zn' | 'Zr'
        static constexpr auto rule(atom_exprZ,     ctll::set<'n', 'r'>) -> Symbol2;
        static constexpr auto rule(atom_exprZ, ctll::neg_set<'n', 'r'>) -> ctll::reject;

        // 'a' | 'as'
        static constexpr auto rule(atom_expr_a,     ctll::set<'s'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2>;
        static constexpr auto rule(atom_expr_a, ctll::neg_set<'s'>) -> ctll::push<pop_char, make_any_aromatic, atom_expr2>;
        // 's' | 'se'
        static constexpr auto rule(atom_expr_s,     ctll::set<'e'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2>;
        static constexpr auto rule(atom_expr_s, ctll::neg_set<'e'>) -> ctll::push<pop_char, make_aromatic, atom_expr2>;

        //
        // Bond expressions
        //

        static constexpr auto rule(bond_expr, ctll::set<'-', '=', '#', '$', ':', '~', '@'>) -> ctll::push<ctll::anything, make_bond_primitive, bond_expr>;
        static constexpr auto rule(bond_expr, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, push_char, bond_expr2>;

        static constexpr auto rule(bond_expr, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr>;
        static constexpr auto rule(bond_expr, ctll::term<'&'>) -> ctll::push<ctll::anything, make_bond_and_high, bond_expr>;
        static constexpr auto rule(bond_expr, ctll::term<','>) -> ctll::push<ctll::anything, make_bond_or, bond_expr>;
        static constexpr auto rule(bond_expr, ctll::term<';'>) -> ctll::push<ctll::anything, make_bond_and_low, bond_expr>;

        static constexpr auto rule(bond_expr, ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '!', '&', ',', ';'>) -> ctll::push<chain>; // FIXME: ring bonds?

        static constexpr auto rule(bond_expr2, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, pop_char, bond_expr>;
        static constexpr auto rule(bond_expr2, ctll::neg_set<'?'>) -> ctll::push<make_bond_primitive, pop_char, bond_expr>;

        //
        // Chain expressions
        //

        static constexpr auto rule(chain, ctll::epsilon) -> ctll::epsilon;        
        static constexpr auto rule(chain, ctll::set<'-', '=', '#', '$', ':', '~', '@', '/', '\\'>) -> ctll::push<bond_expr>;
        //static constexpr auto rule(chain, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, make_bond_primitive, chain_up_down>;
        //static constexpr auto rule(chain_up_down, ctll::neg_set<'?'>) -> ctll::push<chain>;
        //static constexpr auto rule(chain_up_down, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, set_and_high, reset_not, atom>;

        static constexpr auto rule(chain, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr>;
        //static constexpr auto rule(chain, ctll::term<'&'>) -> ctll::push<ctll::anything, set_and_high, bond_expr>;
        //static constexpr auto rule(chain, ctll::term<','>) -> ctll::push<ctll::anything, set_or, bond_expr>;
        //static constexpr auto rule(chain, ctll::term<';'>) -> ctll::push<ctll::anything, set_and_low, bond_expr>;

        static constexpr auto rule(chain, ctll::term<'.'>) -> ctll::push<ctll::anything, reset_prev, chain>;
        static constexpr auto rule(chain, ctll::term<'('>) -> ctll::push<ctll::anything, push_prev, chain>;
        static constexpr auto rule(chain, ctll::term<')'>) -> ctll::push<ctll::anything, pop_prev, chain>;

        // ring bond: DIGIT | '%' DIGIT DIGIT
        static constexpr auto rule(chain, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain>;
        static constexpr auto rule(chain, ctll::term<'%'>) -> ctll::push<ctll::anything, ring_bond>;
        static constexpr auto rule(ring_bond, digit_chars) -> ctll::push<ctll::anything, start_number, ring_bond2>;
        static constexpr auto rule(ring_bond2, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain>;

        // FIXME: add ops  etc
        using not_chain_chars = ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;
        static constexpr auto rule(chain, not_chain_chars) -> ctll::push<atom>;

        // chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom

    };

} // namespace ctsmarts

namespace Kitimar::CTSmarts {

    //
    // Operators
    //

    template<typename Expr>
    struct Not
    {
        constexpr Not() noexcept {}
        constexpr Not(Expr) noexcept {}
    };

    template<typename ...Expr>
    struct And
    {
        constexpr And() noexcept {}
        constexpr And(Expr...) noexcept {}
    };

    template<typename ...Expr>
    struct Or
    {
        constexpr Or() noexcept {}
        constexpr Or(Expr...) noexcept {}
    };

    //
    // Atom primitives
    //

    // '*'
    struct AnyAtom {};

    // 'A'
    struct AnyAliphatic {};

    // 'a'
    struct AnyAromatic {};

    // 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I' | ...
    template<int Element>
    struct AliphaticAtom {};

    // 'b' | 'c' | 'n' | 'o' | 's' | 'p' | 'se' | 'as'
    template<int Element>
    struct AromaticAtom {};

    // NUMBER
    template<int Mass>
    struct Isotope {};

    // '#' NUMBER
    template<int AtomicNumber>
    struct Element {};

    // 'D' | 'D' NUMBER
    template<int N>
    struct Degree {};

    // 'v' | 'v' NUMBER
    template<int N>
    struct Valence {};

    // 'X' | 'X' NUMBER
    template<int N>
    struct Connectivity {};

    // 'H' | 'H' DIGIT
    template<int N>
    struct TotalH {};

    // 'h' | 'h' DIGIT
    template<int N>
    struct ImplicitH {};

    // 'R' | 'r' | 'x'
    struct Cyclic {};

    // 'R0' | 'r0' | 'x0'
    struct Acyclic {};

    // 'R' NUMBER
    template<int N>
    struct RingCount {};

    // 'r' NUMBER
    template<int N>
    struct RingSize {};

    // 'x' NUMBER
    template<int N>
    struct RingConnectivity {};

    // '-' | '-' DIGIT | '+' | '+' DIGIT
    template<int N>
    struct Charge {};

    // '@' '?'? | '@@' '?' |
    // '@TH1' '?'? | '@TH2' '?'? |
    // '@AL1' '?'? | '@AL2' '?'? |
    // '@SP1' '?'? | '@SP2' '?'? | '@SP3' '?'? |
    // '@TB1' '?'? | '@TB2' '?'? | '@TB3' '?'? | ... | '@TB20' '?'?
    // '@OH1' '?'? | '@OH2' '?'? | '@OH3' '?'? | ... | '@OH30' '?'?
    // '@TH?' | '@SP?' | '@AL?' | '@TB?' | '@OH?'
    template<int N>
    struct Chiral {};

    // ':' NUMBER
    template<int N>
    struct Class {};

    //
    // Bond primitives
    //

    // single/aromatic depending on context
    struct ImplicitBond {};

    // '-' | '=' | '#' | '$'
    template<int Order>
    struct BondOrder {};

    // '~'
    struct AnyBond {};

    // ':'
    struct AromaticBond {};

    // '@'
    struct RingBond {};

    // '/'
    struct UpBond {};

    // '\'
    struct DownBond {};

    // '/?' | '\?'
    struct UpOrDownBond {};

    //
    // Graph
    //

    template<int Source, int Target, typename Expr>
    struct Bond
    {
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto expr = Expr();        
    };

} // namespace ctsmarts

#include <cctype>

namespace Kitimar::CTSmarts {

    template<int N, int AtomIndex, typename BondExpr>
    struct RingBondHelper
    {
        static constexpr inline auto n = N;
        static constexpr inline auto atomIndex = AtomIndex;
        static constexpr inline auto bondExpr = BondExpr();
    };

    template<auto N, auto AtomIndex>
    struct ClassHelper
    {
        static constexpr inline auto n = N;
        static constexpr inline auto atomIndex = AtomIndex;
    };

    // Operation tags
    //struct NoOpTag {};
    struct NotTag {};
    struct AndHighTag {};
    struct OrTag {};
    struct AndLowTag {};

    // Error tags

    struct NoErrorTag {};

    struct EmptyBracketAtomTag {}; // '[]'
    struct OpenBracketAtomTag {}; // '[C'
    struct MissingSymbolTag {}; // '[13]' '[+]'
    struct InvalidAtomExprTag {};
    struct InvalidChiralTag {};
    struct InvalidClassTag {}; // '[C:]'
    struct InvalidAtomPrimitiveTag {}; // '[Q]'
    struct InvalidBondPrimitiveTag {}; // 'C^C'
    struct InvalidChiralValenceTag {};
    struct InvalidChiralHydrogenCountTag {};
    struct HydrogenHydrogenCountTag {}; // '[HH1]'

    struct UnmatchedBranchOpeningTag {}; // 'CC(C'
    struct UnmatchedBranchClosingTag {}; // 'CC)C'
    struct InvalidRingBondTag {}; // 'C%' 'C%1'
    struct ConflicingRingBondTag {}; // 'C-1CC=1'
    struct LoopRingBondTag {}; // 'C11'
    struct DuplicateRingBondTag {}; // 'C12CC12'

    struct LeadingDotTag {}; // '.C'
    struct TrailingDotTag {}; // 'C.'
    struct LeadingBondTag {}; // '-C'
    struct TrailingBondTag {}; // 'C-'

    struct MissingLeftOperandTag {}; // '[&C]'
    struct MissingRightOperandTag {}; // '[C&]' '[!]'

    template<int NextIndex, typename PrevIndexT, typename AtomExprT, typename BondExprT, typename RingBondsT, typename ClassesT, typename Error>
    struct SmartsParams
    {
        using PrevIndex = PrevIndexT;
        using AtomExpr = AtomExprT;
        using BondExpr = BondExprT;
        using RingBonds = RingBondsT;
        using Classes = ClassesT;

        static constexpr inline auto nextIndex = NextIndex;
        static constexpr inline auto prevIndex = PrevIndexT();
        static constexpr inline auto atomExpr = AtomExprT();
        static constexpr inline auto bondExpr = BondExprT();
        static constexpr inline auto ringBonds = RingBondsT();
        static constexpr inline auto classes = ClassesT();        
        static constexpr inline auto error = Error();

        constexpr auto pushPrevIndex() const
        {
            constexpr auto index = ctll::front(prevIndex).value;
            return SmartsParams<NextIndex, decltype(ctll::push_front(Number<index>(), prevIndex)), AtomExpr, BondExpr, RingBonds, Classes, Error>();
        }
        constexpr auto popPrevIndex() const
        {
            static_assert(!ctll::empty(prevIndex));
            return SmartsParams<NextIndex, decltype(ctll::pop_front(prevIndex)), AtomExpr, BondExpr, RingBonds, Classes, Error>();
        }
        constexpr auto nextAtom() const
        {
            if constexpr (ctll::empty(prevIndex))
                return SmartsParams<NextIndex + 1, decltype(ctll::push_front(Number<NextIndex>(), prevIndex)), AtomExpr, ctll::empty_list, RingBonds, Classes, Error>();
            else
                return SmartsParams<NextIndex + 1, decltype(ctll::pop_front_and_push_front(Number<NextIndex>(), prevIndex)), AtomExpr, ctll::empty_list, RingBonds, Classes, Error>();
        }
        template<typename Expr>
        constexpr auto setAtomExpr(Expr) const
        {
            return SmartsParams<NextIndex, PrevIndex, Expr, BondExpr, RingBonds, Classes, Error>();
        }
        template<typename Expr>
        constexpr auto setBondExpr(Expr) const
        {
            return SmartsParams<NextIndex, PrevIndex, AtomExpr, Expr, RingBonds, Classes, Error>();
        }
        template<typename Bonds>
        constexpr auto setRingBonds(Bonds) const
        {
            return SmartsParams<NextIndex, PrevIndex, AtomExpr, BondExpr, Bonds, Classes, Error>();
        }
        template<typename Cls>
        constexpr auto setClasses(Cls) const
        {
            return SmartsParams<NextIndex, PrevIndex, AtomExpr, BondExpr, RingBonds, Cls, Error>();
        }
        template<typename Tag>
        constexpr auto setError(Tag) const
        {
            return SmartsParams<NextIndex, PrevIndex, AtomExpr, BondExpr, RingBonds, Classes, Tag>();
        }

    };

    template <typename Atoms = ctll::list<>, typename Bonds = ctll::empty_list, typename ParamsT = SmartsParams<0, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, NoErrorTag>>
    struct SmartsContext
    {
        using Params = ParamsT;

        static constexpr inline auto atoms = Atoms();
        static constexpr inline auto bonds = Bonds();
        static constexpr inline auto params = ParamsT();

        static constexpr inline auto valid = !ctll::size(params.ringBonds);

        constexpr SmartsContext() noexcept {}
        constexpr SmartsContext(Atoms, Bonds, ParamsT) noexcept {}
    };

    template<auto V>
    constexpr auto bondPrimitive(ctll::term<V>)
    {
        if constexpr (V == '-') return BondOrder<1>();
        else if constexpr (V == '=') return BondOrder<2>();
        else if constexpr (V == '#') return BondOrder<3>();
        else if constexpr (V == '$') return BondOrder<4>();
        else if constexpr (V == ':') return AromaticBond();
        else if constexpr (V == '~') return AnyBond();
        else if constexpr (V == '@') return RingBond();
        else if constexpr (V == '/') return UpBond();
        else if constexpr (V == '\\') return DownBond();
        else return ImplicitBond();//static_assert(false);
    }

    template<auto V>
    constexpr auto termAtomicNumber(ctll::term<V>)
    {
        switch (V) {
            case 'H': return 1;
            case 'B': case 'b': return 5;
            case 'C': case 'c': return 6;
            case 'N': case 'n': return 7;
            case 'O': case 'o': return 8;
            case 'F': return 9;
            case 'P': case 'p': return 15;
            case 'S': case 's': return 16;
            case 'K': return 19;
            case 'V': return 23;
            case 'Y': return 39;
            case 'I': return 53;
            case 'W': return 74;
            case 'U': return 92;
            default: return -1;
        }
    }

    template<auto C, auto V>
    constexpr auto symbolAtomicNumber(Char<C>, ctll::term<V>)
    {
        switch (C) {
            case 'a':
                return V == 's' ? 33 : -1;
            case 's':
                return V == 'e' ? 34 : -1;
            case 'A':
                switch (V) {
                    case 'c': return 89;
                    case 'g': return 47;
                    case 'l': return 13;
                    case 'm': return 95;
                    case 'r': return 18;
                    case 's': return 33;
                    case 't': return 85;
                    case 'u': return 79;
                    default: break;
                }
            case 'B':
                switch (V) {
                    case 'a': return 56;
                    case 'e': return 4;
                    case 'i': return 83;
                    case 'k': return 97;
                    case 'r': return 35;
                    default: break;
                }
            case 'C':
                switch (V) {
                    case 'a': return 20;
                    case 'd': return 48;
                    case 'e': return 58;
                    case 'f': return 98;
                    case 'l': return 17;
                    case 'm': return 96;
                    case 'o': return 27;
                    case 'r': return 24;
                    case 's': return 55;
                    case 'u': return 29;
                    default: break;
                }
            case 'D':
                switch (V) {
                    case 'y': return 66;
                    default: break;
                }
            case 'E':
                switch (V) {
                    case 'r': return 68;
                    case 's': return 99;
                    case 'u': return 63;
                    default: break;
                }
            case 'F':
                switch (V) {
                    case 'e': return 26;
                    case 'm': return 100;
                    case 'r': return 87;
                    default: break;
                }
            case 'G':
                switch (V) {
                    case 'a': return 31;
                    case 'd': return 64;
                    case 'e': return 32;
                    default: break;
                }
            case 'H':
                switch (V) {
                    case 'e': return 2;
                    case 'f': return 72;
                    case 'g': return 80;
                    case 'o': return 67;
                    default: break;
                }
            case 'I':
                switch (V) {
                    case 'n': return 49;
                    case 'r': return 77;
                    default: break;
                }
            case 'K':
                switch (V) {
                    case 'r': return 36;
                    default: break;
                }
            case 'L':
                switch (V) {
                    case 'a': return 57;
                    case 'i': return 3;
                    case 'r': return 103;
                    case 'u': return 71;
                    default: break;
                }
            case 'M':
                switch (V) {
                    case 'd': return 101;
                    case 'g': return 12;
                    case 'n': return 25;
                    case 'o': return 42;
                    default: break;
                }
            case 'N':
                switch (V) {
                    case 'a': return 11;
                    case 'b': return 41;
                    case 'd': return 60;
                    case 'e': return 10;
                    case 'i': return 28;
                    case 'o': return 102;
                    case 'p': return 93;
                    default: break;
                }
            case 'O':
                switch (V) {
                    case 's': return 76;
                    default: break;
                }
            case 'P':
                switch (V) {
                    case 'a': return 91;
                    case 'b': return 82;
                    case 'd': return 46;
                    case 'm': return 61;
                    case 'o': return 84;
                    case 'r': return 59;
                    case 't': return 78;
                    case 'u': return 94;
                    default: break;
                }
            case 'R':
                switch (V) {
                    case 'a': return 88;
                    case 'b': return 37;
                    case 'e': return 75;
                    case 'h': return 45;
                    case 'n': return 86;
                    case 'u': return 44;
                    default: break;
                }
            case 'S':
                switch (V) {
                    case 'b': return 51;
                    case 'c': return 21;
                    case 'e': return 34;
                    case 'i': return 14;
                    case 'm': return 62;
                    case 'n': return 50;
                    case 'r': return 38;
                    default: break;
                }
            case 'T':
                switch (V) {
                    case 'a': return 73;
                    case 'b': return 65;
                    case 'c': return 43;
                    case 'e': return 52;
                    case 'h': return 90;
                    case 'i': return 22;
                    case 'l': return 81;
                    case 'm': return 69;
                    default: break;
                }
            case 'X':
                switch (V) {
                    case 'e': return 54;
                    default: break;
                }
            case 'Y':
                switch (V) {
                    case 'b': return 70;
                    default: break;
                }
            case 'Z':
                switch (V) {
                    case 'n': return 30;
                    case 'r': return 40;
                    default: break;
                }
        }
    }

    struct SmartsActions
    {

        //
        // Helper actions for parsing strings & numbers
        //

        // push_char
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_char, ctll::term<V>, Context ctx)
        {
            auto atoms = ctll::push_front(Char<V>(), ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params};
        }

        // pop_char
        template <auto V, auto C, typename ... Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::pop_char, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params> ctx)
        {
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params};
        }

        // start_number
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::start_number, ctll::term<V>, Context ctx)
        {
            auto atoms = ctll::push_front(Number<V - '0'>(), ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params};
        }

        // push_number
        template <auto V, auto N, typename ... Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::push_number, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return SmartsContext{ctll::list<Number<10 * N + V - '0'>, Ts...>(), ctx.bonds, ctx.params};
        }

        //
        // Operator helpers to convert infix expression list with primitives
        // and operator tags to the corresponding AST
        //

        template<typename Expr>
        static constexpr auto isOperation(Expr) { return false; }
        static constexpr auto isOperation(NotTag) { return true; }
        static constexpr auto isOperation(AndHighTag) { return true; }
        static constexpr auto isOperation(OrTag) { return true; }
        static constexpr auto isOperation(AndLowTag) { return true; }

        template<typename ...Expr>
        static constexpr auto isImplicitAnd(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::empty(expr))
                return false;
            else
                return !isOperation(ctll::front(expr));
        }

        static constexpr auto pushExpr(auto list, auto expr)
        {
            if constexpr (isImplicitAnd(list))
                return ctll::push_front(expr, ctll::push_front(AndHighTag(), list));
            else
                return ctll::push_front(expr, list);
        }

        template<typename ...Expr>
        static constexpr auto makeNot(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) % 2 == 0)
                return Not(ctll::front(ctll::rotate(expr)));
            else
                return ctll::front(ctll::rotate(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeAndHighHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeNot(ctll::front(expr));
            else
                return And(makeNot(Expr())...);
        }

        static constexpr auto makeAndHigh(auto expr)
        {
            return makeAndHighHelper(split<AndHighTag>(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeOrHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeAndHigh(ctll::front(expr));
            else
                return Or(makeAndHigh(Expr())...);
        }

        static constexpr auto makeOr(auto expr)
        {
            return makeOrHelper(split<OrTag>(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeAndLowHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeOr(ctll::front(expr));
            else
                return And(makeOr(Expr())...);
        }

        static constexpr auto makeAndLow(auto expr)
        {
            return makeAndLowHelper(split<AndLowTag>(expr));
        }

        static constexpr auto makeAtomAST(auto expr)
        {
            return makeAndLow(expr);
        }

        static constexpr auto makeBondAST(auto expr)
        {
            if constexpr (ctll::empty(expr))
                return ImplicitBond();
            else
                return makeAndLow(expr);
        }

        //
        // Transition to next atom:
        // - convert atomExpr to AST and add it to atoms
        // - create bond if needed
        //

        template<typename Atoms, typename Bonds, typename Params>
        static constexpr auto createBond(Atoms atoms, Bonds bonds, Params params)
        {
            if constexpr (ctll::empty(params.prevIndex)) {
                return SmartsContext{atoms, bonds, params.nextAtom()};
            } else {
                constexpr auto prevIndex = ctll::front(params.prevIndex).value;
                auto expr = makeBondAST(params.bondExpr);
                auto bond = Bond<prevIndex, params.nextIndex, decltype(expr)>();
                auto bonds2 = ctll::push_front(bond, bonds);
                return SmartsContext{atoms, bonds2, params.nextAtom()};
            }
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::next_atom, ctll::term<V> term, Context ctx)
        {
            auto tree = makeAtomAST(ctll::rotate(ctx.params.atomExpr));
            auto atoms = ctll::push_front(tree, ctx.atoms);
            return createBond(atoms, ctx.bonds, ctx.params.setAtomExpr(ctll::empty_list())); // FIXME empty atomExpr in nextAtom
        }

        //
        // Atom primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushAtomExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto atomExpr = pushExpr(ctx.params.atomExpr, leaf);
            auto params = ctx.params.setAtomExpr(atomExpr);
            return SmartsContext{atoms, ctx.bonds, params};
        }

        // make_any_atom
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_atom, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAtom());
        }

        // make_any_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aliphatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAliphatic());
        }

        // make_any_aromatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aromatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAromatic());
        }

        // make_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<termAtomicNumber(term)>());
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params> ctx)
        {
            auto expr = AliphaticAtom<symbolAtomicNumber(Char<C>(), term)>();
            return pushAtomExpr(ctx, ctll::list<Ts...>(), expr);
        }

        // make_aromatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_aromatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AromaticAtom<termAtomicNumber(term)>());
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_aromatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params> ctx)
        {
            auto expr = AromaticAtom<symbolAtomicNumber(Char<C>(), term)>();
            return pushAtomExpr(ctx, ctll::list<Ts...>(), expr);
        }

        // make_isotope
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_isotope, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Isotope<N>());
        }

        // make_element
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Element<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Element<N>());
        }

        // make_degree
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Degree<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Degree<N>());
        }

        // make_valence
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Valence<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Valence<N>());
        }

        // make_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Connectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Connectivity<N>());
        }

        // make_cyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_cyclic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Cyclic());
        }

        // make_acyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_acyclic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Acyclic());
        }

        // make_ring_count
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingCount<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingCount<N>());
        }

        // make_ring_size
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingSize<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingSize<N>());
        }

        // make_ring_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingConnectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingConnectivity<N>());
        }

        // make_class
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_class, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            auto cls = ctll::push_front(ClassHelper<N, ctx.params.nextIndex>(), ctx.params.classes);
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params.setClasses(cls)};
        }

        // make_impl_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_impl_h, ctll::term<V> term, Context ctx)
        {
            constexpr auto value = V == 'h' ? 1 : V - '0';
            return pushAtomExpr(ctx, ctx.atoms, ImplicitH<value>());
        }

        // aliphatic, aromatic, '*', 'A' or 'a'
        template<typename Expr>
        static constexpr bool isTotalHExpr(Expr) { return false; }
        static constexpr bool isTotalHExpr(AnyAtom) { return true; }
        static constexpr bool isTotalHExpr(AnyAliphatic) { return true; }
        static constexpr bool isTotalHExpr(AnyAromatic) { return true; }
        template<int Element>
        static constexpr bool isTotalHExpr(AliphaticAtom<Element>) { return true; }
        template<int Element>
        static constexpr bool isTotalHExpr(AromaticAtom<Element>) { return true; }
        // FIXME: operations
        template<typename Expr>
        static constexpr bool isTotalHExpr(Not<Expr>) { return isTotalHExpr(Expr()); }
        template<typename ...Expr>
        static constexpr bool isTotalHExpr(And<Expr...>) { return (isTotalHExpr(Expr()) || ...); }
        template<typename ...Expr>
        static constexpr bool isTotalHExpr(Or<Expr...>) { return (isTotalHExpr(Expr()) || ...); }        

        static constexpr bool isTotalH(auto atomExpr)
        {
            if constexpr (ctll::empty(atomExpr)) {
                return false;
            } else {
                auto [head, tail] = ctll::pop_and_get_front(atomExpr);
                if (isTotalHExpr(head))
                    return true;
                return isTotalH(tail);
            }
        }

        // make_total_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_total_h, ctll::term<V> term, Context ctx)
        {
            if constexpr (!isTotalH(ctx.params.atomExpr)) {
                return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<1>());
            } else {
                constexpr auto value = V == 'H' ? 1 : V - '0';
                return pushAtomExpr(ctx, ctx.atoms, TotalH<value>());
            }
        }

        // make_neg_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_neg_charge, ctll::term<V> term, Context ctx)
        {
            constexpr auto value = V == '-' ? 1 : V - '0';
            return pushAtomExpr(ctx, ctx.atoms, Charge<-value>());
        }

        // make_pos_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_pos_charge, ctll::term<V> term, Context ctx)
        {
            constexpr auto value = V == '+' ? 1 : V - '0';
            return pushAtomExpr(ctx, ctx.atoms, Charge<value>());
        }

        //
        // Operations
        //

        // make_atom_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_not, ctll::term<V> term, Context ctx)
        {            
            auto atomExpr = pushExpr(ctx.params.atomExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr)};
        }

        static constexpr auto makeAtomOp(auto ctx, auto op)
        {
            auto atomExpr = ctll::push_front(op, ctx.params.atomExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr)};
        }

        // make_atom_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_high, ctll::term<V> term, Context ctx)
        {
            return makeAtomOp(ctx, AndHighTag());
        }

        // make_atom_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_or, ctll::term<V> term, Context ctx)
        {            
            return makeAtomOp(ctx, OrTag());
        }

        // set_atom_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_low, ctll::term<V> term, Context ctx)
        {
            return makeAtomOp(ctx, AndLowTag());
        }

        // make_bond_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_not, ctll::term<V> term, Context ctx)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr)};
        }

        static constexpr auto makeBondOp(auto ctx, auto op)
        {
            auto bondExpr = ctll::push_front(op, ctx.params.bondExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr)};
        }

        // make_bond_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_high, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, AndHighTag());
        }

        // make_bond_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_or, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, OrTag());
        }

        // set_bond_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_low, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, AndLowTag());
        }

        //
        // Branches
        //

        // push_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.pushPrevIndex()};
        }

        // pop_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::pop_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex()};
        }

        // reset_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::reset_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex()};
        }

        //
        // Bond primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushBondExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, leaf);
            auto params = ctx.params.setBondExpr(bondExpr);
            return SmartsContext{atoms, ctx.bonds, params};
        }

        // bond_primitive
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V> term, Context ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            return pushBondExpr(ctx, ctx.atoms, bondPrimitive(term));
            //auto expr = leafBondExpr<ctx.params.notSet>(ctx.params.operation, ctx.params.bondExpr, bondPrimitive(term));
            //return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr)};
            //return ctx;

        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params> ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            if constexpr (V == '?') {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, UpOrDownBond());
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr)};
            } else {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, bondPrimitive(ctll::term<C>()));
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr)};
            }
        }

        //
        // Ring bonds
        //

        template<int N, typename ...RingBonds>
        static constexpr auto findRingBond(ctll::list<RingBonds...> ringBonds)
        {
            if constexpr (ctll::empty(ringBonds)) {
                return ctll::_nothing();
            } else {
                auto rb = ctll::front(ringBonds);
                if constexpr (std::is_same_v<decltype(rb), ctll::_nothing>)
                    return findRingBond<N>(ctll::pop_front(ringBonds));
                if constexpr (rb.n == N)
                    return rb;
                else
                    return findRingBond<N>(ctll::pop_front(ringBonds));
            }
        }

        // ring_bond
        template <auto N, typename Context, typename Atoms>
        static constexpr auto handleRingBond(Context ctx, Atoms atoms)
        {
            constexpr auto atomIndex = ctll::front(ctx.params.prevIndex).value;
            auto rb = findRingBond<N>(ctx.params.ringBonds);
            if constexpr (std::is_same_v<decltype(rb), ctll::_nothing>) {
                auto rb2 = RingBondHelper<N, atomIndex, decltype(ctx.params.bondExpr)>();
                auto ringBonds = ctll::push_front(rb2, ctx.params.ringBonds);
                return SmartsContext{atoms, ctx.bonds, ctx.params.setRingBonds(ringBonds)};
            } else {
                constexpr auto prevIndex = rb.atomIndex;

                auto expr = makeBondAST(ctx.params.bondExpr);
                auto bond = Bond<atomIndex, prevIndex, decltype(expr)>();
                auto bonds = ctll::push_front(bond, ctx.bonds);

                //using BondExpr = decltype(ctx.params.bondExpr);
                //using Expr = ctll::conditional<std::is_same_v<BondExpr, ImplicitBond>, BondOrder<1>, BondExpr>;
                //auto bonds = ctll::push_front(Bond<atomIndex, prevIndex, BondExpr>(), ctx.bonds);
                //auto ringBonds = ctll::pop_front(ctx.params.ringBonds); // FIXME: erase rb.... NOT FRONT!!
                auto ringBonds = ctll::remove_item(rb, ctx.params.ringBonds); // FIXME: erase rb.... NOT FRONT!!
                return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list())};
            }
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V> term, Context ctx)
        {
            return handleRingBond<V - '0'>(ctx, ctx.atoms);
        }

        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params> ctx)
        {
            return handleRingBond<10 * N + V - '0'>(ctx, ctll::list<Ts...>());
        }

        //
        // Errors
        //

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::error_empty_bracket, ctll::term<V> term, Context ctx)
        {
            //return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.template setError<EmptyBracketAtomTag>(EmptyBracketAtomTag())};
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setError(EmptyBracketAtomTag())};
        }

    };

} // namespace ctsmarts

#include <ctll/parser.hpp>

#include <ranges>

namespace Kitimar::CTSmarts {

    template <ctll::fixed_string SMARTS, bool IgnoreInvalid = false>
    struct Smarts;

    template<int I = 0>
    constexpr auto rotateAdjacencyList(auto adjList)
    {
        if constexpr (I == ctll::size(adjList))
            return adjList;
        else {
            auto newAdjList = replace<I, decltype(ctll::rotate(get<I>(adjList)))>(adjList);
            return rotateAdjacencyList<I + 1>(newAdjList);

        }
    }

    //
    // Adjacency list
    //

    template<int NumAtoms, int NumBonds, typename ...Bonds, typename AdjList = decltype(resize<NumAtoms, ctll::list<>>())>
    constexpr auto adjacencyList(ctll::list<Bonds...> bonds, AdjList adjList = {})
    {
        if constexpr (ctll::empty(bonds))
            return adjList;
        else {
            constexpr auto source = ctll::front(bonds).source;
            constexpr auto target = ctll::front(bonds).target;
            constexpr auto bondIndex = NumBonds - ctll::size(ctll::list<Bonds...>());
            //using SourceNbrs = decltype(ctll::push_front(Number<target>(), get<source>(adjList)));
            //using TargetNbrs = decltype(ctll::push_front(Number<source>(), get<target>(adjList)));
            using SourceNbrs = decltype(ctll::push_front(Number<bondIndex>(), get<source>(adjList)));
            using TargetNbrs = decltype(ctll::push_front(Number<bondIndex>(), get<target>(adjList)));

            auto newAdjList = replace<source, SourceNbrs>(replace<target, TargetNbrs>(adjList));
            return adjacencyList<NumAtoms, NumBonds>(ctll::pop_front(bonds), newAdjList);
        }
    }

    //
    // Depth-first search bonds
    //

    template<int Source, int Target, bool IsCyclic, bool IsRingClosure, typename SourceExpr, typename TargetExpr, typename BondExpr>
    struct DfsBond
    {
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
        static constexpr inline auto sourceExpr = SourceExpr();
        static constexpr inline auto targetExpr = TargetExpr();
        static constexpr inline auto bondExpr = BondExpr();
    };

    struct DfsSearch
    {

        static constexpr inline auto NoAtomVisitor = [] (auto, auto, auto ctx) { return ctx; };
        static constexpr inline auto NoBondVisitor = [] (auto, auto, auto, auto, auto, auto ctx) { return ctx; };
        static constexpr inline auto NoAtomBacktrack = [] (auto ctx) { return ctx; };
        static constexpr inline auto NoBondBacktrack = [] (auto ctx) { return ctx; };

        static constexpr auto visitAtom(auto smarts, auto atomIdx, auto atomVisitor, auto ctx, auto visitedAtoms) noexcept
        {
            if constexpr (!ctll::exists_in(atomIdx, visitedAtoms))
                return atomVisitor(smarts, atomIdx, ctx);
            else
                return ctx;
        }

        static constexpr auto backtrackAtom(auto cond, auto atomBacktrack, auto ctx) noexcept
        {
            if constexpr (cond.value)
                return atomBacktrack(ctx);
            else
                return ctx;
        }

        template<int sourceIdx, int adjIdx, typename AtomVisitor, typename BondVisitor, typename Context, typename VB, typename VA>
        static constexpr auto visit(auto smarts, AtomVisitor atomVisitor, BondVisitor bondVisitor, auto atomBacktrack, auto bondBacktrack, Context ctx, VB visitedBonds, VA visitedAtoms) noexcept
        {
            auto adjBondIdxs = get<sourceIdx>(smarts.adjList);
            if constexpr (adjIdx >= ctll::size(adjBondIdxs)) {
                // A leaf node has been reached
                return std::make_tuple(visitedBonds, visitedAtoms, ctx);
            } else {
                auto bondIdx = get<adjIdx>(adjBondIdxs);
                if constexpr (ctll::exists_in(bondIdx, visitedBonds)) {
                    // Skip visited bonds
                    return visit<sourceIdx, adjIdx + 1>(smarts, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx, visitedBonds, visitedAtoms);;
                } else {
                    auto bond = get<bondIdx.value>(smarts.bonds);
                    constexpr auto targetIdx = bond.source == sourceIdx ? bond.target : bond.source;
                    auto targetExpr = get<targetIdx>(smarts.atoms);
                    constexpr auto isRingClosure = ctll::exists_in(std::integral_constant<int, targetIdx>(), visitedAtoms);
                    // Visit atoms and bond
                    auto ctx2 = visitAtom(smarts, std::integral_constant<int, sourceIdx>{}, atomVisitor, ctx, visitedAtoms);
                    auto ctx3  = bondVisitor(smarts, std::integral_constant<int, sourceIdx>{}, std::integral_constant<int, targetIdx>{}, bond.expr, std::integral_constant<bool, isRingClosure>{}, ctx2);
                    auto ctx4 = visitAtom(smarts, std::integral_constant<int, targetIdx>{}, atomVisitor, ctx3, visitedAtoms);
                    // Mark bond and atoms as visited
                    auto visitedBonds2 = ctll::push_front(bondIdx, visitedBonds);
                    auto visitedAtoms2 = ctll::push_front(std::integral_constant<int, sourceIdx>{}, ctll::push_front(std::integral_constant<int, targetIdx>{}, visitedAtoms));
                    // Recursive search....
                    auto [visitedBonds3, visitedAtoms3, ctx5] = visit<targetIdx, 0>(smarts, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx4, visitedBonds2, visitedAtoms2);
                    auto ctx6 = backtrackAtom(std::integral_constant<bool, !isRingClosure>{}, atomBacktrack, ctx5);
                    auto ctx7 = bondBacktrack(ctx6);
                    // bond to next nbr
                    return visit<sourceIdx, adjIdx + 1>(smarts, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx7, visitedBonds3, visitedAtoms3);
                }
            }
        }

        template<typename Context = ctll::empty_list>
        static constexpr auto visit(auto smarts, auto atomVisitor, auto bondVisitor, auto atomBacktrack, auto bondBacktrack, Context ctx = {}) noexcept
        {
            auto ctx2 = std::get<2>(visit<0, 0>(smarts, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx, ctll::empty_list{}, ctll::empty_list{}));
            return atomBacktrack(ctx2);
        }

        template<typename Context = ctll::empty_list>
        static constexpr auto visit(auto smarts, auto atomVisitor, auto bondVisitor, Context ctx = {}) noexcept
        {
            return visit(smarts, atomVisitor, bondVisitor, NoAtomBacktrack, NoBondBacktrack, ctx);
        }

    };

    constexpr auto getDfsAtoms(auto smarts) noexcept
    {
        constexpr auto atomVisitor = [] (auto smarts, auto atomIdx, auto ctx) {
            return ctll::push_front(atomIdx, ctx);
        };
        return ctll::rotate(DfsSearch::visit(smarts, atomVisitor, DfsSearch::NoBondVisitor));
    }

    //template<int BondIdx = 0, typename Atoms = ctll::empty_list, typename Bonds = ctll::empty_list, typename Path = ctll::empty_list>
    template<int BondIdx, typename Atoms, typename Bonds, typename Path>
    struct CycleMembershipContext
    {
        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto path = Path{};
        static constexpr inline auto bondIdx = BondIdx;

        constexpr CycleMembershipContext() noexcept = default;
        constexpr CycleMembershipContext(Atoms, Bonds, Path, Number<BondIdx>) noexcept {}

    };

    template<int SourceIdx, int TargetIdx, int BondIdx>
    struct CylceMembershipBond
    {
        static constexpr inline auto source = SourceIdx;
        static constexpr inline auto target = TargetIdx;
        static constexpr inline auto idx = BondIdx;
    };

    constexpr auto getCycleMembership(auto targetIdx, ctll::empty_list path, auto ctx) noexcept
    {
        return ctx;
    }

    template<typename Bond, typename ...Bonds>
    constexpr auto getCycleMembership(auto targetIdx, ctll::list<Bond, Bonds...> path, auto ctx) noexcept
    {
        auto atoms = ctll::add_item(Number<Bond::target>{}, ctll::add_item(Number<Bond::source>{}, ctx.atoms));
        auto bonds = ctll::add_item(Number<Bond::idx>{}, ctx.bonds);
        auto ctx2 = CycleMembershipContext{atoms, bonds, ctx.path, Number<ctx.bondIdx>{}};
        if constexpr (Bond::idx != ctx.bondIdx && (Bond::source == targetIdx.value || Bond::target == targetIdx.value))
            return ctx2;
        else
            return getCycleMembership(targetIdx, ctll::list<Bonds...>{}, ctx2);
    }

    constexpr auto getCycleMembership(auto smarts) noexcept
    {
        constexpr auto bondVisitor = [] (auto smarts, auto sourceIdx, auto targetIdx, auto expr, auto isRingClosure, auto ctx) {
            auto path = ctll::push_front(CylceMembershipBond<sourceIdx.value, targetIdx.value, ctx.bondIdx>{}, ctx.path);

            if constexpr (isRingClosure.value) {
                auto ctx2 = getCycleMembership(targetIdx, path, ctx);
                return CycleMembershipContext(ctx2.atoms, ctx2.bonds, path, Number<ctx2.bondIdx+1>{});
            } else {
                return CycleMembershipContext(ctx.atoms, ctx.bonds, path, Number<ctx.bondIdx+1>{});
            }
        };

        constexpr auto bondBacktrack = [] (auto ctx) {
            auto path = ctll::pop_front(ctx.path);
            return CycleMembershipContext{ctx.atoms, ctx.bonds, path, Number<ctx.bondIdx>{}};
        };

        auto ctx = CycleMembershipContext<0, ctll::empty_list, ctll::empty_list, ctll::empty_list>{};
        return DfsSearch::visit(smarts, DfsSearch::NoAtomVisitor, bondVisitor, DfsSearch::NoAtomBacktrack, bondBacktrack, ctx);
    }

    template<int BondIdx>
    constexpr auto isCyclic(ctll::empty_list)
    {
        return false;
    }

    template<int BondIdx, typename Bond, typename ...Bonds>
    constexpr auto isCyclic(ctll::list<Bond, Bonds...>)
    {
        if (Bond::value == BondIdx)
            return true;
        return isCyclic<BondIdx>(ctll::list<Bonds...>{});
    }

    template<int BondIdx>
    constexpr auto addCycleMembership(ctll::empty_list, auto cyclicBondIdxs) noexcept
    {
        return ctll::empty_list{};
    }

    template<int BondIdx, int S, int T, bool IC, bool IRC, typename SE, typename TE, typename BE, typename ...CycleBonds>
    constexpr auto addCycleMembership(ctll::list<DfsBond<S, T, IC, IRC, SE, TE, BE>, CycleBonds...>, auto cyclicBondIdxs) noexcept
    {
        auto dfsBond = DfsBond<S, T, isCyclic<BondIdx>(cyclicBondIdxs), IRC, SE, TE, BE>{};
        return ctll::push_front(dfsBond, addCycleMembership<BondIdx+1>(ctll::list<CycleBonds...>{}, cyclicBondIdxs));
    }

    constexpr auto getDfsBonds(auto smarts) noexcept
    {
        constexpr auto bondVisitor = [] (auto smarts, auto sourceIdx, auto targetIdx, auto expr, auto isRingClosure, auto ctx) {
            auto sourceExpr = get<sourceIdx.value>(smarts.atoms);
            auto targetExpr = get<targetIdx.value>(smarts.atoms);
            return ctll::push_front(DfsBond<sourceIdx.value, targetIdx.value, false, isRingClosure.value,
                                            decltype(sourceExpr), decltype(targetExpr), decltype(expr)>(), ctx);
        };
        auto dfsBonds = ctll::rotate(DfsSearch::visit(smarts, DfsSearch::NoAtomVisitor, bondVisitor));
        return addCycleMembership<0>(dfsBonds, getCycleMembership(smarts).bonds);

        //return ctll::rotate(DfsSearch::visit(smarts, DfsSearch::NoAtomVisitor, bondVisitor));
    }

    /*
    template<int sourceIndex, int adjIndex, typename SmartsT, typename ...VBs, typename ...VAs, typename ...Ts>
    constexpr auto getDfsBonds(SmartsT smarts,  ctll::list<VBs...> visitedBonds = {}, ctll::list<VAs...> visitedAtoms = {}, ctll::list<Ts...> bonds = {})
    {
        auto sourceExpr = get<sourceIndex>(smarts.atoms);
        auto adj = get<sourceIndex>(smarts.adjList);
        if constexpr (adjIndex >= ctll::size(adj)) {
            return std::make_tuple(bonds, visitedBonds, visitedAtoms);
        } else {
            auto bondIndex = get<adjIndex>(adj);
            if constexpr (ctll::exists_in(bondIndex, visitedBonds)) {
                return getDfsBonds<sourceIndex, adjIndex + 1>(smarts, visitedBonds, visitedAtoms, bonds);;
            } else {
                auto bond = get<bondIndex.value>(smarts.bonds);
                auto bondExpr = bond.expr;
                constexpr auto targetIndex = bond.source == sourceIndex ? bond.target : bond.source;
                auto targetExpr = get<targetIndex>(smarts.atoms);                
                constexpr auto isRingClosure = ctll::exists_in(Number<targetIndex>(), visitedAtoms);
                // mark as visited
                auto visitedBonds2 = ctll::push_front(bondIndex, visitedBonds);
                auto visitedAtoms2 = ctll::push_front(Number<sourceIndex>(), ctll::push_front(Number<targetIndex>(), visitedAtoms));
                // add dfs bond
                auto bonds2 = ctll::push_front(DfsBond<sourceIndex, targetIndex, isRingClosure, decltype(sourceExpr), decltype(targetExpr), decltype(bondExpr)>(), bonds);
                // dfs search....
                auto [bonds3, visitedBonds3, visitedAtoms3] = getDfsBonds<targetIndex, 0>(smarts, visitedBonds2, visitedAtoms2, bonds2);
                // bond to next nbr
                return getDfsBonds<sourceIndex, adjIndex + 1>(smarts, visitedBonds3, visitedAtoms3, bonds3);
            }
        }
    }

    template<typename SmartsT>
    constexpr auto getDfsBonds(SmartsT smarts)
    {
        //return ctll::rotate(std::get<0>(getDfsBonds<0, 0>(smarts)));
        return getDfsBonds2(smarts);
    }
    */

    //
    // Captures
    //

    template<int N, typename Classes>
    constexpr auto captureIndex(Classes classes)
    {
        if constexpr (ctll::empty(classes)) {
            return -1;
        } else {
            auto cls = ctll::front(classes);
            if (cls.n == N)
                return cls.atomIndex;
            return captureIndex<N>(ctll::pop_front(classes));
        }
    }

    template<int I, typename Classes, typename Mapping>
    constexpr auto captureMappingHelper(Classes classes, Mapping mapping)
    {
        constexpr auto index = captureIndex<I>(classes);
        if constexpr (index < 0) {
            return ctll::rotate(mapping);
        } else {
            return captureMappingHelper<I + 1>(classes, ctll::push_front(Number<index>(), mapping));
        }
    }

    template<typename Smarts>
    constexpr auto captureMapping(Smarts smarts)
    {
        constexpr auto mapping = captureMappingHelper<1>(smarts.context.params.classes, ctll::empty_list());
        return to_array(mapping);
    }

    //
    // Errors
    //

    template<int I>
    struct Pos {};

    struct NoError {};

    template<typename>
    struct EmptyBracketAtomError {};

    template<typename>
    constexpr auto UnmatchedRingBondError() { return false; }

    template<typename ...Ns>
    struct RingBondIds {};

    template<typename ...Ids>
    constexpr auto ringBondIds(ctll::empty_list, ctll::list<Ids...>)
    {
        return RingBondIds<Ids...>();
    }

    template<typename Bond, typename ...Bonds, typename Ids = ctll::empty_list>
    constexpr auto ringBondIds(ctll::list<Bond, Bonds...>, Ids = {})
    {
        return ringBondIds(ctll::list<Bonds...>(), ctll::push_front(Number<Bond::n>(), Ids()));
    }

    //template<typename Atoms, typename Bonds>
    template <ctll::fixed_string SMARTS, bool IgnoreInvalid>
    struct Smarts
    {
        using Result = ctll::parser<SmartsGrammar, SMARTS, SmartsActions>::template output<SmartsContext<>>;
        using Context = Result::output_type;

        static constexpr auto input()
        {
            auto str = SMARTS | std::views::transform([] (auto c) { return static_cast<char>(c); });
            return std::string(str.begin(), str.end());
        }

        static constexpr inline auto context = Context();
        static constexpr inline auto valid = Result::is_correct;
        static constexpr inline auto position = Result::position;

        static constexpr inline auto atoms = ctll::rotate(Context::atoms);
        static constexpr inline auto bonds = ctll::rotate(Context::bonds);
        static constexpr inline auto numAtoms = ctll::size(atoms);
        static constexpr inline auto numBonds = ctll::size(bonds);
        static constexpr inline auto adjList = rotateAdjacencyList(adjacencyList<numAtoms, numBonds>(bonds));

        static constexpr inline auto isSingleAtom = numAtoms == 1;
        static constexpr inline auto isSingleBond = numAtoms == 2 && numBonds == 1;

        template<typename ErrorTag>
        static constexpr auto getError(ErrorTag)
        {
            if constexpr (std::is_same_v<ErrorTag, EmptyBracketAtomTag>)
                return EmptyBracketAtomError<Pos<position>>();
            else if constexpr (!ctll::empty(context.params.ringBonds))
                return UnmatchedRingBondError<decltype(ringBondIds(context.params.ringBonds))>();
            else
                return NoError();
        }

        static constexpr inline auto error = getError(context.params.error);

        //static_assert(IgnoreInvalid || valid);
        static_assert(IgnoreInvalid || std::is_same_v<const NoError, decltype(error)>); // FIXME
        //static_assert(IgnoreInvalid || error); // FIXME

    };

} // namespace ctsmarts

namespace Kitimar::CTSmarts {

    //
    // Operators
    //

    template<typename Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Not<Expr>)
    {
        return !matchAtomExpr(mol, atom, Expr());
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Or<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr()) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, And<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr()) && ...);
    }

    template<typename Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Not<Expr>)
    {
        return !matchBondExpr(mol, bond, Expr());
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Or<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr()) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, And<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr()) && ...);
    }

    //
    // Atoms
    //

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAtom)
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
        return get_connectivity(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, TotalH<N>)
    {
        return get_total_hydrogens(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, ImplicitH<N>)
    {
        return get_implicit_hydrogens(mol, atom) == N;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Cyclic)
    {
        return is_cyclic(mol, atom);
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Acyclic)
    {
        return is_acyclic(mol, atom);
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
        return get_ring_connectivity(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Charge<N>)
    {
        return get_charge(mol, atom) == N;
    }

    //
    // Bonds
    //

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, ImplicitBond)
    {
        return get_order(mol, bond) == 1 || is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, AnyBond)
    {
        return true;
    }

    template<int Order>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, BondOrder<Order>)
    {
        return get_order(mol, bond) == Order && !is_aromatic_bond(mol, bond);
        /*
        if constexpr (requires { is_aromatic(mol, bond); })
            return get_order(mol, bond) == Order && !is_aromatic(mol, bond);
        else
            return get_order(mol, bond) == Order;
        */
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, AromaticBond)
    {
        return is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, RingBond)
    {
        return is_cyclic(mol, bond);
    }

} // namespace ctsmarts

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>
#include <variant>
#include <type_traits>
#include <unordered_set>
#include <cassert>

#define DEBUG_ISOMORPHISM 0

namespace Kitimar::CTSmarts {

    using IsomorphismMapping = std::vector<int>;
    using IsomorphismMappings = std::vector<IsomorphismMapping>;

    template<auto N>
    std::ostream& operator<<(std::ostream &os, const std::array<int, N> &map)
    {
        os << "[";
        for (auto i = 0; i < map.size(); ++i)
            os << " " << map[i];
        os << " ]";
        return os;
    }

    // FIXME
    template<int AtomIndex>
    constexpr auto get_degree(ctll::empty_list)
    {
        return 0;
    }

    template<int AtomIndex, typename Bond, typename ...Bonds>
    constexpr auto get_degree(ctll::list<Bond, Bonds...>)
    {
        int d = Bond::source == AtomIndex || Bond::target == AtomIndex;
        return d + get_degree<AtomIndex>(ctll::list<Bonds...>{});
    }

    template<int AtomIndex, int NumAtoms, typename Bonds>
    constexpr auto get_degrees(Bonds)
    {
        if constexpr (AtomIndex == NumAtoms)
            return ctll::empty_list{};
        else
            return ctll::push_front(Number<get_degree<AtomIndex>(Bonds{})>{}, get_degrees<AtomIndex+1, NumAtoms>(Bonds{}));
    }

    template<typename ...T>
    constexpr auto degreesToArray(ctll::list<T...> degrees)
    {
        return std::array<uint8_t, ctll::size(degrees)>({T::value...});
    }

    template<int NumAtoms, typename Bonds>
    constexpr auto get_degrees(Bonds)
    {
        return degreesToArray(get_degrees<0, NumAtoms>(Bonds{}));

    }

    enum class MapType {
        Single,
        Unique,
        All
    };

    template<MapType T>
    using MapTypeTag = std::integral_constant<MapType, T>;

    static constexpr auto Single = MapTypeTag<MapType::Single>{};
    static constexpr auto Unique = MapTypeTag<MapType::Unique>{};
    static constexpr auto All    = MapTypeTag<MapType::All>{};

    template<typename SmartsT, MapType Type>
    class Isomorphism
    {

        public:
            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto dfsBonds = getDfsBonds(smarts);

            static_assert(ctll::size(smarts.bonds) == ctll::size(dfsBonds));

            Isomorphism(SmartsT, MapTypeTag<Type>)
            {
                m_degrees = get_degrees<smarts.numAtoms>(smarts.bonds);
                m_map.fill(-1);
            }

            //
            // Molecule
            //

            bool match(auto &mol)
            {
                //if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                if (num_atoms(mol) < smarts.numAtoms)
                    return false;
                matchComponents(mol, nullptr);
                return isDone();
            }

            auto count(auto &mol)
            {                
                auto n = 0;
                matchComponents(mol, [&n] (const auto &array) { ++n; });
                return n;
            }

            auto single(auto &mol)
            {
                IsomorphismMapping map;
                matchComponents(mol, [&map] (const auto &array) {
                    map.resize(array.size());
                    std::ranges::copy(array, map.begin());
                });
                return map;
            }

            auto all(auto &mol)
            {
                IsomorphismMappings maps;
                matchComponents(mol, [&maps] (const auto &array) {
                    maps.emplace_back(IsomorphismMapping(array.begin(), array.end()));

                });
                return maps;
            }

            //
            // Atom
            //

            bool matchAtom(auto &mol, const auto &atom)
            {
                reset(mol);
                matchComponent(mol, atom,  nullptr);
                return isDone();
            }

            bool matchBond(auto &mol, const auto &bond)
            {
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);

                if (matchAtom(mol, source, 0, get<0>(smarts.atoms))) {
                    auto index = get_index(mol, source);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    auto callback = [this, &mol, target] (const auto &map) {
                        if (m_map[1] == get_index(mol, target))
                            setDone(true);
                    };
                    matchDfs(mol, callback, dfsBonds);
                    if (isDone() && m_map[1] == get_index(mol, target))
                        return true;
                }

                reset(mol); // FIXME: needed?
                if (matchAtom(mol, target, 0, get<0>(smarts.atoms))) {
                    auto index = get_index(mol, target);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    auto callback = [this, &mol, source] (const auto &map) {
                        if (m_map[1] == get_index(mol, source))
                            setDone(true);
                    };
                    matchDfs(mol, callback, dfsBonds);
                    return isDone() && m_map[1] == get_index(mol, source);
                }

                return false;
            }

            auto count(auto &mol, const auto &atom)
            {
                reset(mol);
                auto n = 0;
                matchComponent(mol, atom, [&n] (const auto &array) { ++n; });
                return n;
            }

            auto single(auto &mol, const auto &atom)
            {
                reset(mol);
                IsomorphismMapping map;
                matchComponent(mol, atom, [&map] (const auto &array) {
                    map.resize(array.size());
                    std::ranges::copy(map, map.begin());
                });
                return map;
            }

            auto all(auto &mol, const auto &atom)
            {
                reset(mol);
                IsomorphismMappings maps;
                matchComponent(mol, atom, [&maps] (const auto &array) {
                    maps.emplace_back(IsomorphismMapping(array.begin(), array.end()));

                });
                return maps;
            }

        private:

            bool matchAtom(auto &mol, const auto &molAtom, int qryAtom, auto atomExpr) noexcept
            {
                if (get_degree(mol, molAtom) < m_degrees[qryAtom])
                    return false;

                return matchAtomExpr(mol, molAtom, atomExpr);
            }

            void matchDfs(auto &mol, auto callback, auto bonds)
            {
                static_assert(!ctll::empty(bonds));
                auto queryBond = ctll::front(bonds);
                auto querySource = queryBond.source;
                auto queryTarget = queryBond.target;

                auto source = get_atom(mol, m_map[querySource]);

                //if (DEBUG_ISOMORPHISM)
                //    std::cout << "mapping bond: " << get_index(m_query, querySource) << "-" << get_index(m_query, queryTarget) << std::endl;

                for (auto bond : get_bonds(mol, source)) {

                    if constexpr (queryBond.isCyclic)
                        if (!is_cyclic_bond(mol, bond))
                            continue;

                    if (!matchBondExpr(mol, bond, queryBond.bondExpr))
                        continue;

                    auto target = Molecule::get_nbr(mol, bond, source);
                    auto targetIndex = get_index(mol, target);

                    if constexpr (queryBond.isRingClosure) {
                        if (m_map[queryTarget] != targetIndex)
                            continue;
                    } else {
                        assert(targetIndex < m_mapped.size());
                        if (m_mapped[targetIndex])
                            continue;

                        // match target atom
                        //if (!matchAtomExpr(mol, target, queryBond.targetExpr))
                        if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr))
                            continue;

                        if constexpr (DEBUG_ISOMORPHISM)
                            std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                        // map target atom
                        m_map[queryTarget] = targetIndex;
                        m_mapped[targetIndex] = true;
                    }

                    if constexpr (ctll::size(bonds) == 1) {
                        // Found mapping...

                        //if (stereoMatches())

                        //if (DEBUG_ISOMORPHISM)
                        //    std::cout << "found map: " << m_map << std::endl;

                        if constexpr (Type == MapType::Unique) {
                            // create bit mask of atoms (to ensure uniqueness of mapping)
                            std::vector<bool> atoms(num_atoms(mol));
                            for (auto index : m_map)
                                atoms[index] = true;
                            // add the mapping to the result if it is unique
                            auto hash = std::hash<std::vector<bool>>()(atoms);
                            if (m_maps.find(hash) == m_maps.end()) {
                                m_maps.insert(hash);
                                addMapping(callback);
                            }
                        } else
                            addMapping(callback);
                    } else {
                        matchDfs(mol, callback, ctll::pop_front(bonds));
                    }

                    // exit as soon as possible if only one match is required
                    // (single mapping stored in m_map after returning)
                    if (isDone())
                        return;

                    // bracktrack target atom
                    if constexpr (!queryBond.isRingClosure) {
                        m_map[queryTarget] = -1;
                        m_mapped[targetIndex] = false;
                    }
                }
            }

            void matchComponent(auto &mol, const auto &atom, auto callback)
            {
                if constexpr (!smarts.numAtoms)
                    return;

                if (isDone())
                    return;

                if constexpr (!smarts.numBonds) {
                    constexpr auto queryAtom = 0;

                    //if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms)))
                    if (!matchAtom(mol, atom, 0, get<0>(smarts.atoms)))
                        return;

                    if constexpr (DEBUG_ISOMORPHISM)
                        std::cout << queryAtom << " -> " << get_index(mol, atom) << '\n';

                    // map source atom, recursive dfs, backtrack
                    auto index = get_index(mol, atom);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    addMapping(callback);
                    m_map[0] = -1;
                    m_mapped[index] = false;

                } else {
                    auto queryBond = ctll::front(dfsBonds);
                    auto queryAtom = queryBond.source;

                    //if (!matchAtomExpr(mol, atom, queryBond.sourceExpr))
                    if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr))
                        return;

                    if constexpr (DEBUG_ISOMORPHISM)
                        std::cout << queryAtom << " -> " << get_index(mol, atom) << '\n';

                    // map source atom, recursive dfs, backtrack
                    auto index = get_index(mol, atom);
                    m_map[queryAtom] = index;
                    m_mapped[index] = true;
                    matchDfs(mol, callback, dfsBonds);
                    m_map[queryAtom] = -1;
                    m_mapped[index] = false;
                }
            }

            void matchComponents(auto &mol, auto callback)
            {
                if constexpr (!smarts.numAtoms)
                    return;

                reset(mol);

                if constexpr (!smarts.numBonds) {
                    for (auto atom : get_atoms(mol))
                        matchComponent(mol, atom, callback);
                } else {
                    // try to match each atom in the molecule against the first atom
                    // epxression in the SMARTS
                    for (auto atom : get_atoms(mol))
                        matchComponent(mol, atom, callback);
                }
            }

            constexpr auto reset(auto &mol) noexcept
            {
                setDone(false);
                m_mapped.clear();
                m_mapped.resize(num_atoms(mol));
            }

            constexpr auto isDone() const noexcept
            {
                return m_done;
            }

            constexpr auto setDone(bool done) noexcept
            {
                m_done = done;
            }

            template<typename Callback>
            constexpr auto addMapping(Callback callback) noexcept
            {
                if constexpr (Type == MapType::Single)
                    setDone(true);
                if constexpr (!std::is_same_v<std::nullptr_t, Callback>)
                    callback(m_map);
            }

            /*
            constexpr auto matchAtom(auto &mol, const auto &atom, std::size_t queryAtomIndex) const noexcept
            {
                return with_n<ctll::size(smarts.atoms), bool>(queryAtomIndex, [&mol, &atom] (auto i) {
                    return matchAtomExpr(mol, atom, get<i>(smarts.atoms));
                });
            }

            constexpr auto matchBond(auto &mol, const auto &bond, std::size_t queryBondIndex) const noexcept
            {
                return with_n<ctll::size(smarts.bonds), bool>(queryBondIndex, [&mol, &bond] (auto i) {
                    return matchBondExpr(mol, bond, get<i>(smarts.bonds).expr);
                });
            }
            */

            std::array<int, smarts.numAtoms> m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms

            //std::array<bool, smarts.numAtoms> m_mapped; // current mapping: queried atom index -> true if mapped
            std::vector<bool> m_mapped; // current mapping: queried atom index -> true if mapped
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

} // namespace ctsmarts

#include <cstdint>
#include <vector>
#include <cassert>

namespace Kitimar::Molecule {

    struct MockAtom
    {
        uint8_t element;
        uint8_t isotope;
        int8_t charge;
        uint8_t degree;
        uint8_t implicitHydrogens;
        bool cyclic;
        bool aromatic;
    };

    struct MockBond
    {
        uint32_t source;
        uint32_t target;
        uint8_t order;
        bool cyclic;
        bool aromatic;
    };

    struct MockMolecule
    {
        std::vector<MockAtom> atoms;
        std::vector<MockBond> bonds;
    };

    //
    // Molecule
    //

    inline uint32_t num_atoms(const MockMolecule &mol) noexcept
    {
        return mol.atoms.size();
    }

    inline uint32_t num_bonds(const MockMolecule &mol) noexcept
    {
        return mol.atoms.size();
    }

    inline auto get_atoms(const MockMolecule &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_atoms(mol))>(0), num_atoms(mol));
    }

    inline auto get_bonds(const MockMolecule &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_bonds(mol))>(0), num_bonds(mol));
    }

    inline auto get_atom(const MockMolecule &mol, uint32_t index) noexcept
    {
        assert(index < num_atoms(mol));
        return index;
    }

    inline auto get_bond(const MockMolecule &mol, uint32_t index) noexcept
    {
        assert(index < num_bonds(mol));
        return index;
    }

    inline auto get_index(const MockMolecule &mol, uint32_t index) noexcept
    {
        return index;
    }

    //
    // Atom
    //

    inline auto get_element(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].element;
    }

    inline auto get_isotope(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].isotope;
    }

    inline auto get_charge(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].charge;
    }

    inline auto get_degree(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].degree;
    }

    inline auto get_implicit_hydrogens(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].implicitHydrogens;
    }

    inline auto is_cyclic_atom(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].cyclic;
    }

    inline auto is_aromatic_atom(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].aromatic;
    }

    inline uint32_t get_source(const MockMolecule &mol, uint32_t bond) noexcept;
    inline uint32_t get_target(const MockMolecule &mol, uint32_t bond) noexcept;

    inline auto get_bonds(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_bonds(mol) | std::views::filter([&mol, atom] (auto i) {
            auto bond = get_bond(mol, i);
            return atom == get_source(mol, bond) || atom == get_target(mol, bond);
        });
    }

    inline auto get_nbrs(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
            return Kitimar::Molecule::get_nbr(mol, bond, atom);
        });
    }

    inline uint32_t null_atom(const MockMolecule &mol) noexcept
    {
        return -1;
    }

    //
    // Bond
    //

    inline uint32_t get_source(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].source;
    }

    inline uint32_t get_target(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].target;
    }

    inline auto get_order(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].order;
    }

    inline auto is_cyclic_bond(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].cyclic;
    }

    inline auto is_aromatic_bond(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].aromatic;
    }

    inline auto null_bond(const MockMolecule &mol) noexcept
    {
        return -1;
    }

} // namespace Kitimar::Molecule

namespace Kitimar::CTSmarts {

    namespace detail {

        template<typename Smarts, auto N>
        auto captureAtoms(Molecule::Molecule auto &mol, Smarts, const IsomorphismMapping &map, const std::array<int, N> &cap)
        {
            using Atom = decltype(get_atom(mol, 0));
            if constexpr (N) {
                std::array<Atom, N> atoms = {};
                if (map.empty())
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < N; ++i)
                        atoms[i] = get_atom(mol, map[cap[i]]);
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

        constexpr bool singleAtomMatch(auto smarts, auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms));
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
            if (matchAtomExpr(mol, source, get<0>(smarts.atoms)) && matchAtomExpr(mol, target, get<1>(smarts.atoms)))
                return 1;
            if (matchAtomExpr(mol, source, get<1>(smarts.atoms)) && matchAtomExpr(mol, target, get<0>(smarts.atoms)))
                return 2;
            return 0;
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

    } // namespace detail

    //
    // CTSmarts::contains<"SMARTS">(mol) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool contains(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (detail::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.match(mol);
        }
    }

    //
    // CTSmarts::atom<"SMARTS">(mol, atom) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool atom(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            return detail::singleAtomMatch(smarts, mol, atom);
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms)))
                return false;
            for (auto bond : get_bonds(mol, atom)) {
                if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                    continue;
                if (matchAtomExpr(mol, get_nbr(mol, bond, atom), get<1>(smarts.atoms)))
                    return true;
            }
            return false;
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.matchAtom(mol, atom);
        }
    }

    //
    // CTSmarts::bond<"SMARTS">(mol, bond) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool bond(Molecule::Molecule auto &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        //std::cout << "CTSmarts::bond<" << smarts.input() << ">(mol, " << get_index(mol, bond) << ")" << std::endl;
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return detail::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism{smarts, All};
            return iso.matchBond(mol, bond);
            /*
            for (const auto &map : iso.all(mol)) {
                if (map[queryBond.source] == sourceIndex && map[queryBond.target] == targetIndex)
                    return true;
                if (map[queryBond.source] == targetIndex && map[queryBond.target] == sourceIndex)
                    return true;
            }
            return false;
            */
        }
    }

    //
    // CTSmarts::count<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::integeral
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    constexpr auto count(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::contains<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            auto n = 0;
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol))
                if (detail::singleBondMatch(smarts, mol, bond))
                    ++n;
            return n;
        } else {
            auto iso = Isomorphism{smarts, mapType};
            return iso.count(mol);
        }
    }

    //
    // CTSmarts::single<"SMARTS">(mol) -> std::vector<int>
    //

    template<ctll::fixed_string SMARTS>
    constexpr auto single(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return IsomorphismMapping{1, atom};
            return IsomorphismMapping{};
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                switch (detail::singleBondMatch(smarts, mol, bond)) {
                    case 1:
                        return IsomorphismMapping{1, get_source(mol, bond), get_target(mol, bond)};
                    case 2:
                        return IsomorphismMapping{1, get_target(mol, bond), get_source(mol, bond)};
                    default:
                        break;
                }
            }
            return IsomorphismMapping{};
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.single(mol);
        }
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
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            constexpr auto cap = captureMapping(smarts);
            for (auto bond : get_bonds(mol)) {
                auto matchType = detail::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return detail::singleBondCapture(smarts, mol, bond, matchType);
            }
            return detail::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism{smarts, Single};
            constexpr auto cap = captureMapping(smarts);
            auto map = iso.single(mol);
            return detail::captureMatchAtoms(mol, smarts, map, cap);
        }
    }

    //
    // CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS>
    auto capture(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};

        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(smarts);
        auto map = iso.single(mol, atom);
        return detail::captureMatchAtoms(mol, smarts, map, cap);
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, map, cap);
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
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol, atom) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol, atom));
    }

} // namespace Kitimar::CTSmarts
