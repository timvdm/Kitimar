#pragma once

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
        static constexpr auto rule(atom_exprR,   ctll::term<'0'>) -> ctll::push<ctll::anything, pop_char, make_acyclic, atom_expr2>;
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
