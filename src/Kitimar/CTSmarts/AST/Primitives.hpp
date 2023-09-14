#pragma once

namespace Kitimar::CTSmarts {

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

} // namespace Kitimar::CTSmarts
