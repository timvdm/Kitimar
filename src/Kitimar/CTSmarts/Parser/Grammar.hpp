#pragma once

#include <ctll/list.hpp>
#include <ctll/grammars.hpp>

namespace Kitimar::CTSmarts {




    struct SmartsGrammar
    {
        template<bool Component = false, typename Recursive = ctll::empty_list, int Branch = 0>
        struct Parenthesis
        {
            static constexpr inline auto component = Component;
            static constexpr inline auto recursive = Recursive{};
            static constexpr inline auto branch = Branch;
        };

        template<typename P>
        using ParenthesisPushComponent = Parenthesis<true, decltype(P::recursive), P::branch>;

        template<typename P>
        using ParenthesisPopComponent = Parenthesis<false, decltype(P::recursive), P::branch>;

        template<typename P>
        using ParenthesisPushRecursive = Parenthesis<P::component, decltype(ctll::push_front(Number<P::branch>{}, P::recursive)), 0>;

        template<typename P>
        using ParenthesisPopRecursive = Parenthesis<P::component, decltype(ctll::pop_front(P::recursive)), ctll::front(P::recursive).value>;

        template<typename P>
        using ParenthesisPushBranch = Parenthesis<P::component, decltype(P::recursive), P::branch + 1>;

        template<typename P>
        using ParenthesisPopBranch = Parenthesis<P::component, decltype(P::recursive), P::branch - 1>;


        //
        // Symbols
        //

        template<typename P> struct atom {};
        template<typename P> struct atom_B {};
        template<typename P> struct atom_C {};
        template<typename P> struct atom_expr {}; // FIXME: rename to bracket_atom
        template<typename P> struct atom_expr2 {}; // FIXME: rename to atom_expr
        template<typename P> struct atom_exprA {}; // atom_expr_A
        template<typename P> struct atom_exprB {};
        template<typename P> struct atom_exprC {};
        template<typename P> struct atom_exprD {};
        template<typename P> struct atom_exprE {};
        template<typename P> struct atom_exprF {};
        template<typename P> struct atom_exprG {};
        template<typename P> struct atom_exprH {};
        template<typename P> struct atom_exprI {};
        template<typename P> struct atom_exprK {};
        template<typename P> struct atom_exprL {};
        template<typename P> struct atom_exprM {};
        template<typename P> struct atom_exprN {};
        template<typename P> struct atom_exprO {};
        template<typename P> struct atom_exprP {};
        template<typename P> struct atom_exprR {};
        template<typename P> struct atom_exprS {};
        template<typename P> struct atom_exprT {};
        template<typename P> struct atom_exprX {};
        template<typename P> struct atom_exprY {};
        template<typename P> struct atom_exprZ {};
        template<typename P> struct atom_expr_a {};
        template<typename P> struct atom_expr_s {};
        template<typename P> struct atom_expr_isotope {};
        template<typename P> struct atom_expr_isotope2 {};
        template<typename P> struct atom_expr_element {};
        template<typename P> struct atom_expr_element2 {};
        template<typename P> struct atom_expr_degree {};
        template<typename P> struct atom_expr_valence {};
        template<typename P> struct atom_expr_valence2 {};
        template<typename P> struct atom_expr_connectivity {};
        template<typename P> struct atom_expr_total_h {};
        template<typename P> struct atom_expr_impl_h {};
        template<typename P> struct atom_expr_ring_count {};
        template<typename P> struct atom_expr_ring_size {};
        template<typename P> struct atom_expr_ring_size2 {};
        template<typename P> struct atom_expr_ring_connectivity {};
        template<typename P> struct atom_expr_ring_connectivity2 {};
        template<typename P> struct atom_expr_neg_charge {};
        template<typename P> struct atom_expr_neg_charge2 {};
        template<typename P> struct atom_expr_pos_charge {};
        template<typename P> struct atom_expr_pos_charge2 {};
        template<typename P> struct atom_expr_chiral {};
        template<typename P> struct atom_expr_class {};
        template<typename P> struct atom_expr_class2 {};
        template<typename P> struct atom_expr_recursive {};
        template<typename P> struct bond_expr {};
        template<typename P> struct bond_expr2 {};
        template<typename P> struct chain {};
        template<typename P> struct chain_up_down {};
        template<typename P> struct ring_bond {};
        template<typename P> struct ring_bond2 {};

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
        struct start_charge : ctll::action {};
        struct increment_charge : ctll::action {};
        struct decrement_charge : ctll::action {};
        struct make_charge : ctll::action {};
        struct make_chiral : ctll::action {};
        struct make_class : ctll::action {};

        struct make_bond_primitive : ctll::action {};
        //struct make_up_or_down_bond : ctll::action {};

        // Create bond AST elements
        struct next_atom : ctll::action {};


        // Components
        struct push_component : ctll::action {};
        struct pop_component : ctll::action {};

        // Branches
        struct push_prev : ctll::action {};
        struct pop_prev : ctll::action {};
        struct reset_prev : ctll::action {};
        struct set_bond_type : ctll::action {};
        struct handle_ring_bond : ctll::action {};



        // Recursive
        struct push_recursive : ctll::action {};
        struct pop_recursive : ctll::action {};

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


        using _start = atom<Parenthesis<>>;








        //
        // Organic atoms
        //

        // aliphatic atom: 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'N', 'O', 'P', 'S', 'F', 'I'>) -> ctll::push<ctll::anything, make_aliphatic, next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_B<P>>;
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_C<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::term<'r'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::neg_set<'r'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::term<'l'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::neg_set<'l'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        // aromatic atom: 'b' | 'c' | 'n' | 'o' | 's' | 'p'
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'b', 'c','n','o','p','s'>) -> ctll::push<ctll::anything, make_aromatic, next_atom, chain<P>>;
        // any atom: '*'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'*'>) -> ctll::push<ctll::anything, make_any_atom, next_atom, chain<P>>;
        // any aromatic: 'a'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'a'>) -> ctll::push<ctll::anything, make_any_aromatic, next_atom, chain<P>>;
        // any aliphatic: 'A'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'A'>) -> ctll::push<ctll::anything, make_any_aliphatic, next_atom, chain<P>>;
        // bracket atom: '[' atom_expression+ ']'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'['>) -> ctll::push<ctll::anything, atom_expr<P>>;


        //
        // Atom expressions (i.e. [ ... ] )
        //

        using digit_chars = ctll::range<'0', '9'>;
        using not_digit_chars = ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;

        // check for [], [!] or [...!]
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::set<']'>) -> ctll::push<error_empty_bracket, ctll::reject>;
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::neg_set<']', '!'>) -> ctll::push<atom_expr2<P>>;

        // operations + end ]
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'&'>) -> ctll::push<ctll::anything, make_atom_and_high, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<','>) -> ctll::push<ctll::anything, make_atom_or, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<';'>) -> ctll::push<ctll::anything, make_atom_and_low, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<']'>) -> ctll::push<ctll::anything, next_atom, chain<P>>;

        // isotope: NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_isotope<P>>;
        template<typename P> static constexpr auto rule(atom_expr_isotope<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_isotope<P>>;
        template<typename P> static constexpr auto rule(atom_expr_isotope<P>, not_digit_chars) -> ctll::push<make_isotope, atom_expr2<P>>;

        // element: '#' | '#' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'#'>) -> ctll::push<ctll::anything, atom_expr_element<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_element2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element<P>, not_digit_chars) -> ctll::push<make_element, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_element2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element2<P>, not_digit_chars) -> ctll::push<make_element, atom_expr2<P>>;

        // valence: 'v' | 'v' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'v'>) -> ctll::push<ctll::anything, atom_expr_valence<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_valence2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence<P>, not_digit_chars) -> ctll::push<make_valence, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_valence2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence2<P>, not_digit_chars) -> ctll::push<make_valence, atom_expr2<P>>;

        // implicit hydrogens: 'h' | 'h' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'h'>) -> ctll::push<ctll::anything, atom_expr_impl_h<P>>;
        template<typename P> static constexpr auto rule(atom_expr_impl_h<P>, digit_chars) -> ctll::push<ctll::anything, make_impl_h, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_impl_h<P>, not_digit_chars) -> ctll::push<make_impl_h, atom_expr2<P>>;

        // cyclic: 'r'
        // acyclic: 'r0'
        // ring size: 'r' NUMBER
        // FIXME: r1 & r2  = error??
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'r'>) -> ctll::push<ctll::anything, atom_expr_ring_size<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_size2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_size2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size2<P>, not_digit_chars) -> ctll::push<make_ring_size, atom_expr2<P>>;

        // cyclic: 'x'
        // acyclic: 'x0'
        // ring connectivity: 'x' NUMBER
        // FIXME: x1 = error?
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'x'>) -> ctll::push<ctll::anything, atom_expr_ring_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_connectivity2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_connectivity2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity2<P>, not_digit_chars) -> ctll::push<make_ring_connectivity, atom_expr2<P>>;

        // charge: '-' | '--' | '---' | ... | '-' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, start_charge, atom_expr_neg_charge<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, digit_chars) -> ctll::push<ctll::anything, make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, decrement_charge, atom_expr_neg_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, ctll::neg_set<'-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge2<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, decrement_charge, atom_expr_neg_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge2<P>, ctll::neg_set<'-'>) -> ctll::push<make_charge, atom_expr2<P>>;

        // charge: '+' | '++' | '+++' | ... | '+' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, start_charge, atom_expr_pos_charge<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, digit_chars) -> ctll::push<ctll::anything, make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, increment_charge, atom_expr_pos_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, ctll::neg_set<'+', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge2<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, increment_charge, atom_expr_pos_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge2<P>, ctll::neg_set<'+'>) -> ctll::push<make_charge, atom_expr2<P>>;

        // atom class: ':' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<':'>) -> ctll::push<ctll::anything, atom_expr_class<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_class2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class<P>, not_digit_chars) -> ctll::reject; //ctll::push<make_class, set_and_high, atom_expr2>;
        template<typename P> static constexpr auto rule(atom_expr_class2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_class2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class2<P>, not_digit_chars) -> ctll::push<make_class, atom_expr2<P>>;

        // symbol: 'U' | 'V' | 'W'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'U', 'V', 'W'>) -> ctll::push<ctll::anything, make_aliphatic, atom_expr2<P>>;
        // symbol: 'b' | 'c' | 'n' | 'o' | 'p'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'b', 'c', 'n', 'o', 'p'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        // any atom: '*'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'*'>) -> ctll::push<ctll::anything, make_any_atom, atom_expr2<P>>;

        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'A'>) -> ctll::push<ctll::anything, push_char, atom_exprA<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_exprB<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_exprC<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'D'>) -> ctll::push<ctll::anything, push_char, atom_exprD<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'E'>) -> ctll::push<ctll::anything, push_char, atom_exprE<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'F'>) -> ctll::push<ctll::anything, push_char, atom_exprF<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'G'>) -> ctll::push<ctll::anything, push_char, atom_exprG<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'H'>) -> ctll::push<ctll::anything, push_char, atom_exprH<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'I'>) -> ctll::push<ctll::anything, push_char, atom_exprI<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'K'>) -> ctll::push<ctll::anything, push_char, atom_exprK<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'L'>) -> ctll::push<ctll::anything, push_char, atom_exprL<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'M'>) -> ctll::push<ctll::anything, push_char, atom_exprM<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'N'>) -> ctll::push<ctll::anything, push_char, atom_exprN<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'O'>) -> ctll::push<ctll::anything, push_char, atom_exprO<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'P'>) -> ctll::push<ctll::anything, push_char, atom_exprP<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'R'>) -> ctll::push<ctll::anything, push_char, atom_exprR<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'S'>) -> ctll::push<ctll::anything, push_char, atom_exprS<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'T'>) -> ctll::push<ctll::anything, push_char, atom_exprT<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'X'>) -> ctll::push<ctll::anything, push_char, atom_exprX<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'Y'>) -> ctll::push<ctll::anything, push_char, atom_exprY<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'Z'>) -> ctll::push<ctll::anything, push_char, atom_exprZ<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'a'>) -> ctll::push<ctll::anything, push_char, atom_expr_a<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'s'>) -> ctll::push<ctll::anything, push_char, atom_expr_s<P>>;

        template<typename P> using Symbol1 = ctll::push<pop_char, make_aliphatic, atom_expr2<P>>;
        template<typename P> using Symbol2 = ctll::push<ctll::anything, make_aliphatic, atom_expr2<P>>;

        // symbol: 'Ac' | 'Ag' | 'Al' | 'Am' | 'Ar' | 'As' | 'At' | 'Au'
        // any aliphatic: 'A'
        template<typename P> static constexpr auto rule(atom_exprA<P>,     ctll::set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprA<P>, ctll::neg_set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> ctll::push<pop_char, make_any_aliphatic, atom_expr2<P>>;
        // symbol: 'B' | 'Ba' | 'Be' | 'Bh' | 'Bi' | 'Bk' | 'Br'
        template<typename P> static constexpr auto rule(atom_exprB<P>,     ctll::set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprB<P>, ctll::neg_set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol1<P>;
        // C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu
        template<typename P> static constexpr auto rule(atom_exprC<P>,     ctll::set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprC<P>, ctll::neg_set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol1<P>;
        // symbol: 'Db' | 'Ds' | 'Dy'
        // degree: 'D' | 'D' NUMBER
        template<typename P> static constexpr auto rule(atom_exprD<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_degree<P>>;
        template<typename P> static constexpr auto rule(atom_exprD<P>,     ctll::set<'b', 's', 'y'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprD<P>, ctll::neg_set<'b', 's', 'y', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_degree, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_degree<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_degree<P>>;
        template<typename P> static constexpr auto rule(atom_expr_degree<P>, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_degree, atom_expr2<P>>;
        // symbol: 'Er' | 'Es' | 'Eu'
        template<typename P> static constexpr auto rule(atom_exprE<P>,     ctll::set<'r', 's', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprE<P>, ctll::neg_set<'r', 's', 'u'>) -> ctll::reject;
        // symbol: 'F' | 'Fe' | 'Fm' | 'Fr'
        template<typename P> static constexpr auto rule(atom_exprF<P>,     ctll::set<'e', 'm', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprF<P>, ctll::neg_set<'e', 'm', 'r'>) -> Symbol1<P>;
        // symbol: 'Ga' | 'Gd' | 'Ge'
        template<typename P> static constexpr auto rule(atom_exprG<P>,     ctll::set<'a', 'd', 'e'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprG<P>, ctll::neg_set<'a', 'd', 'e'>) -> ctll::reject;
        // symbol: 'H' | 'He' | 'Hf' | 'Hg' | 'Ho' | 'Hs'
        // total hydrogens: 'H' | 'H' DIGIT
        template<typename P> static constexpr auto rule(atom_exprH<P>,     ctll::set<'e', 'f', 'g', 'o', 's'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprH<P>,   ctll::range<'0', '9'>) -> ctll::list<ctll::anything, pop_char, make_total_h, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_exprH<P>, ctll::neg_set<'e', 'f', 'g', 'o', 's', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::list<pop_char, make_total_h, atom_expr2<P>>;
        // symbol: 'I' | 'In' | 'Ir'
        template<typename P> static constexpr auto rule(atom_exprI<P>,     ctll::set<'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprI<P>, ctll::neg_set<'n', 'r'>) -> Symbol1<P>;
        // symbol: 'K' | 'Kr'
        template<typename P> static constexpr auto rule(atom_exprK<P>,    ctll::term<'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprK<P>, ctll::neg_set<'r'>) -> Symbol1<P>;
        // La Li Lr Lu
        template<typename P> static constexpr auto rule(atom_exprL<P>,     ctll::set<'a', 'i', 'r', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprL<P>, ctll::neg_set<'a', 'i', 'r', 'u'>) -> ctll::reject;
        // Md Mg Mn Mo Mt
        template<typename P> static constexpr auto rule(atom_exprM<P>,     ctll::set<'d', 'g', 'n', 'o', 't'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprM<P>, ctll::neg_set<'d', 'g', 'n', 'o', 't'>) -> ctll::reject;
        // N Na Nb Nd Ne Ni No Np
        template<typename P> static constexpr auto rule(atom_exprN<P>,     ctll::set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprN<P>, ctll::neg_set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol1<P>;
        // O Os
        template<typename P> static constexpr auto rule(atom_exprO<P>,    ctll::term<'s'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprO<P>, ctll::neg_set<'s'>) -> Symbol1<P>;
        // P Pa Pb Pd Pm Po Pr Pt Pu
        template<typename P> static constexpr auto rule(atom_exprP<P>,     ctll::set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprP<P>, ctll::neg_set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol1<P>;
        // symbol: 'Ra' | 'Rb' | 'Re' | 'Rf' | 'Rg' | 'Rh' | 'Rn' | 'Ru'
        // cyclic: 'R'
        // acyclic: 'R0'
        // ring count: 'R' NUMBER
        template<typename P> static constexpr auto rule(atom_exprR<P>,   ctll::term<'0'>) -> ctll::push<ctll::anything, pop_char, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_exprR<P>,   ctll::range<'1', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_ring_count<P>>;
        template<typename P> static constexpr auto rule(atom_exprR<P>,     ctll::set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprR<P>, ctll::neg_set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_count<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_ring_count<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_count<P>, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_ring_count, atom_expr2<P>>;
        // S Sb Sc Se Sg Si Sm Sn Sr
        template<typename P> static constexpr auto rule(atom_exprS<P>,     ctll::set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprS<P>, ctll::neg_set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol1<P>;
        // Ta Tb Tc Te Th Ti Tl Tm
        template<typename P> static constexpr auto rule(atom_exprT<P>,     ctll::set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprT<P>, ctll::neg_set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> ctll::reject;
        // symbol: 'Xe'
        // connectivity: 'X' | 'X' NUMBER
        template<typename P> static constexpr auto rule(atom_exprX<P>, digit_chars) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_exprX<P>,    ctll::term<'e'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprX<P>, ctll::neg_set<'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_connectivity, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_connectivity<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_expr_connectivity<P>, not_digit_chars) -> ctll::push<make_connectivity, atom_expr2<P>>;
        // symbol: 'Y' | 'Yb'
        template<typename P> static constexpr auto rule(atom_exprY<P>,    ctll::term<'b'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprY<P>, ctll::neg_set<'b'>) -> Symbol1<P>;
        // symbol: 'Zn' | 'Zr'
        template<typename P> static constexpr auto rule(atom_exprZ<P>,     ctll::set<'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprZ<P>, ctll::neg_set<'n', 'r'>) -> ctll::reject;

        // 'a' | 'as'
        template<typename P> static constexpr auto rule(atom_expr_a<P>,     ctll::set<'s'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_a<P>, ctll::neg_set<'s'>) -> ctll::push<pop_char, make_any_aromatic, atom_expr2<P>>;
        // 's' | 'se'
        template<typename P> static constexpr auto rule(atom_expr_s<P>,     ctll::set<'e'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_s<P>, ctll::neg_set<'e'>) -> ctll::push<pop_char, make_aromatic, atom_expr2<P>>;

        //
        // Recursive SMARTS
        //

        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'$'>) -> ctll::push<ctll::anything, atom_expr_recursive<P>>;
        template<typename P> static constexpr auto rule(atom_expr_recursive<P>, ctll::term<'('>) -> ctll::push<ctll::anything, push_recursive, atom<ParenthesisPushRecursive<P>>>;
        template<typename P> static constexpr auto rule(atom_expr_recursive<P>, ctll::neg_set<'('>) -> ctll::reject;

        //
        // Bond expressions
        //

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::set<'-', '=', '#', '$', ':', '~', '@'>) -> ctll::push<ctll::anything, make_bond_primitive, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, push_char, bond_expr2<P>>;

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<'&'>) -> ctll::push<ctll::anything, make_bond_and_high, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<','>) -> ctll::push<ctll::anything, make_bond_or, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<';'>) -> ctll::push<ctll::anything, make_bond_and_low, bond_expr<P>>;

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '!', '&', ',', ';'>) -> ctll::push<chain<P>>; // FIXME: ring bonds?

        template<typename P> static constexpr auto rule(bond_expr2<P>, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, pop_char, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr2<P>, ctll::neg_set<'?'>) -> ctll::push<make_bond_primitive, pop_char, bond_expr<P>>;




        //
        // Chain expressions
        //



        template<typename P> static constexpr auto rule(chain<P>, ctll::epsilon) -> ctll::epsilon;
        template<typename P> static constexpr auto rule(chain<P>, ctll::set<'-', '=', '#', '$', ':', '~', '@', '/', '\\'>) -> ctll::push<bond_expr<P>>;
        //template<typename P> static constexpr auto rule(chain, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, make_bond_primitive, chain_up_down<P>>;
        //template<typename P> static constexpr auto rule(chain_up_down, ctll::neg_set<'?'>) -> ctll::push<chain<P>>;
        //template<typename P> static constexpr auto rule(chain_up_down, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, set_and_high, reset_not, atom<P>>;

        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr<P>>;

        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'.'>) -> ctll::push<ctll::anything, reset_prev, chain<P>>;
        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'('>) -> ctll::push<ctll::anything, push_prev, chain<P>>;

        //
        // Close parenthesis
        //

        // pop_recusrive
        template<bool Component, typename R, typename ...Rs>
        static constexpr auto rule(chain<Parenthesis<Component, ctll::list<R, Rs...>, 0>>, ctll::term<')'>) -> ctll::push<ctll::anything, pop_recursive, atom_expr2<Parenthesis<Component, ctll::list<Rs...>, R::value>>>;

        // pop_prev
        template<bool Component, typename Recursive, int Branch>
        static constexpr auto rule(chain<Parenthesis<Component, Recursive, Branch>>, ctll::term<')'>) -> ctll::push<ctll::anything, pop_prev, chain<Parenthesis<Component, Recursive, Branch - 1>>>;

        // pop_component
        static constexpr auto rule(chain<Parenthesis<true, ctll::empty_list, 0>>, ctll::term<')'>) -> ctll::push<ctll::anything, pop_prev, chain<Parenthesis<>>>;






        // ring bond: DIGIT | '%' DIGIT DIGIT
        template<typename P> static constexpr auto rule(chain<P>, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain<P>>;
        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'%'>) -> ctll::push<ctll::anything, ring_bond<P>>;
        template<typename P> static constexpr auto rule(ring_bond<P>, digit_chars) -> ctll::push<ctll::anything, start_number, ring_bond2<P>>;
        template<typename P> static constexpr auto rule(ring_bond2<P>, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain<P>>;



        // FIXME: add ops  etc
        using not_chain_chars = ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;
        template<typename P> static constexpr auto rule(chain<P>, not_chain_chars) -> ctll::push<atom<P>>;




        // chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom



    };

} // namespace ctsmarts
