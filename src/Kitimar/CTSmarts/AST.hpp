#pragma once

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    //
    // Operators
    //

    template<typename Expr>
    struct Not
    {
        static constexpr inline auto expr = Expr{};

        constexpr Not() noexcept {}
        constexpr Not(Expr) noexcept {}
    };

    template<typename ...Expr>
    struct And
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr And() noexcept {}
        constexpr And(Expr...) noexcept {}
        constexpr And(ctll::list<Expr...>) noexcept {}
    };

    template<typename ...Expr>
    struct Or
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr Or() noexcept {}
        constexpr Or(Expr...) noexcept {}
        constexpr Or(ctll::list<Expr...>) noexcept {}
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

    // 'h'
    struct HasImplicitH {};

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
    struct Charge
    {
        static constexpr inline int value = N;
    };

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

    template<int Index, typename Expr>
    struct Atom
    {
        static constexpr inline auto index = Index;
        static constexpr inline auto expr = Expr{};
    };

    template<int Index, int Source, int Target, typename Expr>
    struct Bond
    {
        static constexpr inline auto index = Index;
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto expr = Expr{};
    };

    //template<typename Atoms, typename Bonds = ctll::empty_list, typename Classes = ctll::empty_list>
    template<typename Atoms, typename Bonds, typename Classes>
    struct BasicSmarts
    {
        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto classes = Classes{};
        static constexpr inline auto numAtoms = ctll::size(atoms);
        static constexpr inline auto numBonds = ctll::size(bonds);
        static constexpr inline auto numClasses = ctll::size(classes);

        static constexpr inline auto isSingleAtom = numAtoms == 1;
        static constexpr inline auto isSingleBond = numAtoms == 2 && numBonds == 1;

        consteval BasicSmarts() noexcept = default;
        consteval BasicSmarts(Atoms, Bonds, Classes) noexcept {}
    };

} // namespace ctsmarts
