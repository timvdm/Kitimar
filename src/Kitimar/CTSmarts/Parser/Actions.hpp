#pragma once

#include "Grammar.hpp"
#include "../AST/AST.hpp"
#include "../Util/Ctll.hpp"

#include <type_traits>
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


    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    struct SmartsParams
    {
        using NextIndex = T1;
        using PrevIndex = T2;
        using AtomExpr = T3;
        using BondExpr = T4;
        using RingBonds = T5;
        using Classes = T6;
        using Error = T7;

        static constexpr inline auto nextIndex = NextIndex{}; // Number
        static constexpr inline auto prevIndex = PrevIndex{}; // ctll::list<Number>
        static constexpr inline auto atomExpr = AtomExpr{}; // ctll::list<T>
        static constexpr inline auto bondExpr = BondExpr{};  // ctll::list<T>
        static constexpr inline auto ringBonds = RingBonds{}; // ctll::list<RingBondHelper>
        static constexpr inline auto classes = Classes{}; // ctll::list<ClassHelper>
        static constexpr inline auto error = Error{};

        template<typename T> static consteval auto setNextIndex(T) noexcept { return SmartsParams<T, T2, T3, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setPrevIndex(T) noexcept { return SmartsParams<T1, T, T3, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setAtomExpr(T) noexcept  { return SmartsParams<T1, T2, T, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setBondExpr(T) noexcept  { return SmartsParams<T1, T2, T3, T, T5, T6, T7>{}; }
        template<typename T> static consteval auto setRingBonds(T) noexcept { return SmartsParams<T1, T2, T3, T4, T, T6, T7>{}; }
        template<typename T> static consteval auto setClasses(T) noexcept   { return SmartsParams<T1, T2, T3, T4, T5, T, T7>{}; }
        template<typename T> static consteval auto setError(T) noexcept     { return SmartsParams<T1, T2, T3, T4, T5, T6, T>{}; }

        consteval SmartsParams() noexcept = default;

        constexpr auto pushPrevIndex() const
        {
            return setPrevIndex(ctll::push_front(ctll::front(prevIndex), prevIndex));
        }
        constexpr auto popPrevIndex() const
        {
            static_assert(!ctll::empty(prevIndex));
            return setPrevIndex(ctll::pop_front(prevIndex));
        }

        template<typename PrevIndexT>
        constexpr auto nextAtomHelper(PrevIndexT) const
        {
            return SmartsParams<Number<nextIndex.value + 1>, PrevIndexT, AtomExpr, ctll::empty_list, RingBonds, Classes, Error>{};
        }

        constexpr auto nextAtom() const
        {
            if constexpr (ctll::empty(prevIndex))
                return nextAtomHelper(ctll::push_front(nextIndex, prevIndex));
            else
                return nextAtomHelper(ctll::pop_front_and_push_front(nextIndex, prevIndex));
        }

    };


    template <typename Atoms = ctll::empty_list, typename Bonds = ctll::empty_list,
              typename ParamsT = SmartsParams<Number<0>, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, NoErrorTag>,
              typename ParentT = ctll::_nothing>
    struct SmartsContext
    {
        using Params = ParamsT;
        using Parent = ParentT;

        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto params = ParamsT{};
        static constexpr inline auto parent = ParentT{};

        static constexpr inline auto valid = !ctll::size(params.ringBonds);

        constexpr SmartsContext() noexcept {}
        constexpr SmartsContext(Atoms, Bonds, Params, Parent) noexcept {}
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
                break;
            case 'B':
                switch (V) {
                    case 'a': return 56;
                    case 'e': return 4;
                    case 'i': return 83;
                    case 'k': return 97;
                    case 'r': return 35;
                    default: break;
                }
                break;
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
                break;
            case 'D':
                switch (V) {
                    case 'y': return 66;
                    default: break;
                }
                break;
            case 'E':
                switch (V) {
                    case 'r': return 68;
                    case 's': return 99;
                    case 'u': return 63;
                    default: break;
                }
                break;
            case 'F':
                switch (V) {
                    case 'e': return 26;
                    case 'm': return 100;
                    case 'r': return 87;
                    default: break;
                }
                break;
            case 'G':
                switch (V) {
                    case 'a': return 31;
                    case 'd': return 64;
                    case 'e': return 32;
                    default: break;
                }
                break;
            case 'H':
                switch (V) {
                    case 'e': return 2;
                    case 'f': return 72;
                    case 'g': return 80;
                    case 'o': return 67;
                    default: break;
                }
                break;
            case 'I':
                switch (V) {
                    case 'n': return 49;
                    case 'r': return 77;
                    default: break;
                }
                break;
            case 'K':
                switch (V) {
                    case 'r': return 36;
                    default: break;
                }
                break;
            case 'L':
                switch (V) {
                    case 'a': return 57;
                    case 'i': return 3;
                    case 'r': return 103;
                    case 'u': return 71;
                    default: break;
                }
                break;
            case 'M':
                switch (V) {
                    case 'd': return 101;
                    case 'g': return 12;
                    case 'n': return 25;
                    case 'o': return 42;
                    default: break;
                }
                break;
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
                break;
            case 'O':
                switch (V) {
                    case 's': return 76;
                    default: break;
                }
                break;
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
                break;
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
                break;
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
                break;
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
                break;
            case 'X':
                switch (V) {
                    case 'e': return 54;
                    default: break;
                }
                break;
            case 'Y':
                switch (V) {
                    case 'b': return 70;
                    default: break;
                }
                break;
            case 'Z':
                switch (V) {
                    case 'n': return 30;
                    case 'r': return 40;
                    default: break;
                }
                break;
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
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // pop_char
        template <auto V, auto C, typename ... Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::pop_char, ctll::term<V>, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params, ctx.parent};
        }

        // start_number
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::start_number, ctll::term<V>, Context ctx)
        {
            auto atoms = ctll::push_front(Number<V - '0'>(), ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // push_number
        template <auto V, auto N, typename ... Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::push_number, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return SmartsContext{ctll::list<Number<10 * N + V - '0'>, Ts...>(), ctx.bonds, ctx.params, ctx.parent};
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
            constexpr auto expr2 = makeAndLow(expr);
            if constexpr (SmartsActions::isTotalHExpr(expr2))
                return SmartsActions::toTotalHExpr(expr2);
            else
                return expr2;
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



        template<typename Atoms, typename Bonds, typename Params, typename Parent>
        static constexpr auto makeBond(Atoms atoms, Bonds bonds, Params params, Parent parent)
        {
            if constexpr (ctll::empty(params.prevIndex)) {
                return SmartsContext{atoms, bonds, params.nextAtom(), parent};
            } else {
                constexpr auto prevIndex = ctll::front(params.prevIndex).value;
                auto expr = makeBondAST(ctll::rotate(params.bondExpr));
                auto bond = Bond<ctll::size(bonds), prevIndex, params.nextIndex(), decltype(expr)>();
                auto bonds2 = ctll::push_front(bond, bonds);
                return SmartsContext{atoms, bonds2, params.nextAtom(), parent};
            }
        }


        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::next_atom, ctll::term<V>, Context ctx)
        {
            auto expr = makeAtomAST(ctll::rotate(ctx.params.atomExpr));
            auto atom = Atom<ctll::size(ctx.atoms), decltype(expr)>{};
            auto atoms = ctll::push_front(atom, ctx.atoms);
            return makeBond(atoms, ctx.bonds, ctx.params.setAtomExpr(ctll::empty_list()), ctx.parent); // FIXME empty atomExpr in nextAtom
        }

        //
        // Atom primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushAtomExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto atomExpr = pushExpr(ctx.params.atomExpr, leaf);
            auto params = ctx.params.setAtomExpr(atomExpr);
            return SmartsContext{atoms, ctx.bonds, params, ctx.parent};
        }

        // make_any_atom
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_atom, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAtom());
        }

        // make_any_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aliphatic, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAliphatic());
        }

        // make_any_aromatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aromatic, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAromatic());
        }

        // make_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<termAtomicNumber(term)>());
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
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
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_aromatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            auto expr = AromaticAtom<symbolAtomicNumber(Char<C>(), term)>();
            return pushAtomExpr(ctx, ctll::list<Ts...>(), expr);
        }

        // make_isotope
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_isotope, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Isotope<N>());
        }


        // make_element
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Element<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Element<N>());
        }

        // make_degree
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Degree<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Degree<N>());
        }

        // make_valence
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Valence<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Valence<N>());
        }

        // make_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Connectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Connectivity<N>());
        }

        // make_cyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_cyclic, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Cyclic());
        }

        // make_acyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_acyclic, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Acyclic());
        }

        // make_ring_count
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingCount<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingCount<N>());
        }

        // make_ring_size
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingSize<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingSize<N>());
        }

        // make_ring_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingConnectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingConnectivity<N>());
        }

        // make_class
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_class, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            auto cls = ctll::push_front(ClassHelper<N, ctx.params.nextIndex()>(), ctx.params.classes);
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params.setClasses(cls), ctx.parent};
        }

        // make_has_impl_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_has_impl_h, ctll::term<V>, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, HasImplicitH{});
        }

        // make_impl_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_impl_h, ctll::term<V>, Context ctx)
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
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(AliphaticAtom<AtomicNumber>) { return true; }
        static constexpr bool isTotalHExpr(AliphaticAtom<1>) { return false; }
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(AromaticAtom<AtomicNumber>) { return true; }
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(Element<AtomicNumber>) { return true; }
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

        template<typename Expr>
        static constexpr auto toTotalHExpr(Expr) { return Expr{}; }
        static constexpr auto toTotalHExpr(AliphaticAtom<1>) { return TotalH<1>{}; }
        template<typename Expr>
        static constexpr auto toTotalHExpr(Not<Expr>) { return Not<decltype(toTotalHExpr(Expr{}))>{}; }



        template<typename ...Expr>
        static constexpr auto toTotalHExpr(And<Expr...> op) { return And(toTotalHExpr(op.expr)); }
        template<typename ...Expr>
        static constexpr auto toTotalHExpr(Or<Expr...> op) { return Or(toTotalHExpr(op.expr)); }

        template<typename ...Expr>
        static constexpr auto toTotalHExpr(ctll::list<Expr...> expr)
        {
            return transform(expr, [] <typename Expr2> (Expr2 expr2) {
                if constexpr (std::is_same_v<Expr2, AliphaticAtom<1>>)
                    return TotalH<1>{};
                else
                    return toTotalHExpr(expr2);
            });
        }


        // make_total_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_total_h, ctll::term<V>, Context ctx)
        {
            if constexpr (V == 'H')
                // Will be replaced by TotalH later if needed
                return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<1>());
            else
                return pushAtomExpr(ctx, ctx.atoms, TotalH<V - '0'>());
        }

        // start_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::start_charge, ctll::term<V>, Context ctx)
        {
            constexpr auto value = V == '+' ? 1 : -1;
            auto atoms = ctll::push_front(Charge<value>{}, ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        template <int Value, typename Context>
        static constexpr auto add_charge(Context ctx)
        {
            auto [charge, tail] = ctll::pop_and_get_front(ctx.atoms);
            auto atoms = ctll::push_front(Charge<charge.value + Value>{}, tail);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // increment_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::increment_charge, ctll::term<V>, Context ctx)
        {
            return add_charge<1>(ctx);
        }

        // decrement_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::decrement_charge, ctll::term<V>, Context ctx)
        {
            return add_charge<-1>(ctx);
        }

        // make_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_charge, ctll::term<V>, Context ctx)
        {
            auto [charge, atoms] = ctll::pop_and_get_front(ctx.atoms);
            constexpr bool isDigit = V >= '0' && V <= '9';
            if constexpr (isDigit)
                return pushAtomExpr(ctx, atoms, Charge<charge.value * (V - '0')>{});
            else
                return pushAtomExpr(ctx, atoms, charge);
        }







        //
        // Operations
        //



        // make_atom_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_not, ctll::term<V>, Context ctx)
        {
            auto atomExpr = pushExpr(ctx.params.atomExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr), ctx.parent};
        }

        static constexpr auto makeAtomOp(auto ctx, auto op)
        {
            auto atomExpr = ctll::push_front(op, ctx.params.atomExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr), ctx.parent};
        }

        // make_atom_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_high, ctll::term<V>, Context ctx)
        {
            return makeAtomOp(ctx, AndHighTag());
        }

        // make_atom_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_or, ctll::term<V>, Context ctx)
        {
            return makeAtomOp(ctx, OrTag());
        }

        // set_atom_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_low, ctll::term<V>, Context ctx)
        {
            return makeAtomOp(ctx, AndLowTag());
        }

        // make_bond_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_not, ctll::term<V>, Context ctx)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr), ctx.parent};
        }

        static constexpr auto makeBondOp(auto ctx, auto op)
        {
            auto bondExpr = ctll::push_front(op, ctx.params.bondExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr), ctx.parent};
        }

        // make_bond_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_high, ctll::term<V>, Context ctx)
        {
            return makeBondOp(ctx, AndHighTag());
        }

        // make_bond_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_or, ctll::term<V>, Context ctx)
        {
            return makeBondOp(ctx, OrTag());
        }

        // set_bond_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_low, ctll::term<V>, Context ctx)
        {
            return makeBondOp(ctx, AndLowTag());
        }

        //
        // Branches
        //

        // push_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_prev, ctll::term<V>, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.pushPrevIndex(), ctx.parent};
        }

        // pop_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::pop_prev, ctll::term<V>, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex(), ctx.parent};
        }

        // reset_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::reset_prev, ctll::term<V>, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex(), ctx.parent};
        }

        //
        // Recursive
        //

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_recursive, ctll::term<V>, Context ctx)
        {
            return SmartsContext{{}, {}, {}, ctx};
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::pop_recursive, ctll::term<V>, Context ctx)
        {
            auto expr = BasicSmarts{ctll::rotate(Context::atoms), ctll::rotate(Context::bonds), Context::params.classes};
            return pushAtomExpr(ctx.parent, ctx.parent.atoms, expr);
        }

        //
        // Bond primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushBondExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, leaf);
            auto params = ctx.params.setBondExpr(bondExpr);
            return SmartsContext{atoms, ctx.bonds, params, ctx.parent};
        }

        // bond_primitive
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V> term, Context ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            return pushBondExpr(ctx, ctx.atoms, bondPrimitive(term));
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V>, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            if constexpr (V == '?') {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, UpOrDownBond());
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr), ctx.parent};
            } else {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, bondPrimitive(ctll::term<C>()));
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr), ctx.parent};
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

        template<int Index, int Source, int Target, typename Expr1, typename Expr2>
        static constexpr auto makeRingBond(Expr1, Expr2)
        {
            if constexpr (std::is_same_v<Expr1, Expr2>)
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr1, ImplicitBond>)
                return std::make_tuple(Bond<Index, Source, Target, Expr2>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr2, ImplicitBond>)
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, NoErrorTag{});
            else
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, ConflicingRingBondTag{});
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
                return SmartsContext{atoms, ctx.bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()), ctx.parent};
            } else {
                constexpr auto prevIndex = rb.atomIndex;
                constexpr auto ringBond = makeRingBond<ctll::size(ctx.bonds), atomIndex, prevIndex>(makeBondAST(ctll::rotate(rb.bondExpr)), makeBondAST(ctll::rotate(ctx.params.bondExpr)));
                constexpr auto bond = std::get<0>(ringBond);
                constexpr auto error = std::get<1>(ringBond);
                auto bonds = ctll::push_front(bond, ctx.bonds);
                auto ringBonds = ctll::remove_item(rb, ctx.params.ringBonds);
                if constexpr (std::is_same_v<NoErrorTag, decltype(error)>)
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()), ctx.parent};
                else
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()).setError(error), ctx.parent};
            }
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V>, Context ctx)
        {
            return handleRingBond<V - '0'>(ctx, ctx.atoms);
        }

        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V>, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return handleRingBond<10 * N + V - '0'>(ctx, ctll::list<Ts...>());
        }

        //
        // Errors
        //

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::error_empty_bracket, ctll::term<V>, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setError(EmptyBracketAtomTag()), ctx.parent};
        }




    };

} // namespace ctsmarts
