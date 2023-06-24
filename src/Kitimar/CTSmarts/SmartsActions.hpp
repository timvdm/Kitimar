#pragma once

#include "SmartsGrammar.hpp"
#include "SmartsAST.hpp"
#include "Util.hpp"

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

        template<int Source, int Target, typename Expr1, typename Expr2>
        static constexpr auto makeRingBond(Expr1, Expr2)
        {
            if constexpr (std::is_same_v<Expr1, Expr2>)
                return std::make_tuple(Bond<Source, Target, Expr1>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr1, ImplicitBond>)
                return std::make_tuple(Bond<Source, Target, Expr2>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr2, ImplicitBond>)
                return std::make_tuple(Bond<Source, Target, Expr1>{}, NoErrorTag{});
            else
                return std::make_tuple(Bond<Source, Target, Expr1>{}, ConflicingRingBondTag{});
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
                auto [bond, error] = makeRingBond<atomIndex, prevIndex>(makeBondAST(rb.bondExpr), makeBondAST(ctx.params.bondExpr));
                auto bonds = ctll::push_front(bond, ctx.bonds);
                auto ringBonds = ctll::remove_item(rb, ctx.params.ringBonds); // FIXME: erase rb.... NOT FRONT!!
                if constexpr (std::is_same_v<NoErrorTag, decltype(error)>)
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list())};
                else
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()).setError(error)};
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
