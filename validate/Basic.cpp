#include "Validate.hpp"

template<Molecule::Molecule Mol>
void validate_basic(Mol &mol)
{
    //
    // Atom primitives
    //

    // Any...
    validate<"*">(mol);
    validate<"[*]">(mol);
    validate<"A">(mol);
    validate<"[A]">(mol);
    validate<"a">(mol);
    validate<"[a]">(mol);

    // AliphaticAtom
    validate<"B">(mol);
    validate<"C">(mol);
    validate<"N">(mol);
    validate<"O">(mol);
    validate<"S">(mol);
    validate<"P">(mol);
    validate<"F">(mol);
    validate<"Cl">(mol);
    validate<"Br">(mol);
    validate<"I">(mol);

    // AromaticAtom
    validate<"b">(mol);
    validate<"c">(mol);
    validate<"n">(mol);
    validate<"o">(mol);
    validate<"s">(mol);
    validate<"p">(mol);
    validate<"[se]">(mol);
    validate<"[as]">(mol);

    // Isotope
    // Element

    // Degree
    validate<"[D]">(mol);
    validate<"[D0]">(mol);
    validate<"[D1]">(mol);
    validate<"[D2]">(mol);
    validate<"[D3]">(mol);
    validate<"[D4]">(mol);
    validate<"[D5]">(mol);
    validate<"[D6]">(mol);
    validate<"[D7]">(mol);
    validate<"[D8]">(mol);
    validate<"[D9]">(mol);

    // Valence
    validate<"[v]">(mol);
    validate<"[v0]">(mol);
    validate<"[v1]">(mol);
    validate<"[v2]">(mol);
    validate<"[v3]">(mol);
    validate<"[v4]">(mol);
    validate<"[v5]">(mol);
    validate<"[v6]">(mol);
    validate<"[v7]">(mol);
    validate<"[v8]">(mol);
    validate<"[v9]">(mol);

    // Connectivity
    validate<"[X]">(mol);
    if constexpr (!std::is_same_v<Mol, OpenBabel::OBMol>) // OB: 'X0' == 'X1' ? (parsmarts.cpp:1061)
        validate<"[X0]">(mol);
    validate<"[X1]">(mol);
    validate<"[X2]">(mol);
    validate<"[X3]">(mol);
    validate<"[X4]">(mol);
    validate<"[X5]">(mol);
    validate<"[X6]">(mol);
    validate<"[X7]">(mol);
    validate<"[X8]">(mol);
    validate<"[X9]">(mol);

    // TotalH
    validate<"[H]">(mol);
    validate<"[*H]">(mol);
    validate<"[H0]">(mol);
    validate<"[H1]">(mol);
    validate<"[H2]">(mol);
    validate<"[H3]">(mol);
    validate<"[*H0]">(mol);
    validate<"[*H1]">(mol);
    validate<"[*H2]">(mol);
    validate<"[*H3]">(mol);
    validate<"[*H4]">(mol);
    validate<"[*H5]">(mol);
    validate<"[*H6]">(mol);
    validate<"[*H7]">(mol);
    validate<"[*H8]">(mol);
    validate<"[*H9]">(mol);

    // ImplicitH
    validate<"[h]">(mol);
    validate<"[h0]">(mol);
    validate<"[h1]">(mol);
    validate<"[h2]">(mol);
    validate<"[h3]">(mol);
    validate<"[h4]">(mol);
    validate<"[h5]">(mol);
    validate<"[h6]">(mol);
    validate<"[h7]">(mol);
    validate<"[h8]">(mol);
    validate<"[h9]">(mol);

    // Cyclic
    validate<"[R]">(mol);
    validate<"[r]">(mol);
    validate<"[x]">(mol);

    // Acyclic
    validate<"[R0]">(mol);
    validate<"[r0]">(mol);
    validate<"[x0]">(mol);

    // RingCount
    validate<"[R]">(mol);
    validate<"[R0]">(mol);
    validate<"[R1]">(mol);
    validate<"[R2]">(mol);
    validate<"[R3]">(mol);
    validate<"[R4]">(mol);
    validate<"[R5]">(mol);
    validate<"[R6]">(mol);
    validate<"[R7]">(mol);
    validate<"[R8]">(mol);
    validate<"[R9]">(mol);

    // RingSize
    validate<"[r]">(mol);
    validate<"[r0]">(mol);
    validate<"[r1]">(mol);
    validate<"[r2]">(mol);
    validate<"[r3]">(mol);
    validate<"[r4]">(mol);
    validate<"[r5]">(mol);
    validate<"[r6]">(mol);
    validate<"[r7]">(mol);
    validate<"[r8]">(mol);
    validate<"[r9]">(mol);

    // RingConnectivity
    validate<"[x]">(mol);
    validate<"[x0]">(mol);
    validate<"[x1]">(mol);
    validate<"[x2]">(mol);
    validate<"[x3]">(mol);
    validate<"[x4]">(mol);
    validate<"[x5]">(mol);
    validate<"[x6]">(mol);
    validate<"[x7]">(mol);
    validate<"[x8]">(mol);
    validate<"[x9]">(mol);

    // Charge
    validate<"[-]">(mol);
    validate<"[--]">(mol);
    validate<"[+]">(mol);
    validate<"[++]">(mol);
    validate<"[-9]">(mol);
    validate<"[-8]">(mol);
    validate<"[-7]">(mol);
    validate<"[-6]">(mol);
    validate<"[-5]">(mol);
    validate<"[-4]">(mol);
    validate<"[-3]">(mol);
    validate<"[-2]">(mol);
    validate<"[-1]">(mol);
    validate<"[-0]">(mol);
    validate<"[+0]">(mol);
    validate<"[+1]">(mol);
    validate<"[+2]">(mol);
    validate<"[+3]">(mol);
    validate<"[+4]">(mol);
    validate<"[+5]">(mol);
    validate<"[+6]">(mol);
    validate<"[+7]">(mol);
    validate<"[+8]">(mol);
    validate<"[+9]">(mol);


    //validate<"[+,++,+++]">(mol); RDKit does not support '+++'?
    validate<"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]">(mol);
    validate<"[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]">(mol);
    validate<"[Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]">(mol);
    validate<"[!H]">(mol);
    validate<"[R](-*(-*))~*~*~*~[a]">(mol);
/*
FAIL: [+,++,+++]
FAIL: [$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]
FAIL: [AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]
FAIL: [Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]
FAIL: [!H]
FAIL: [R](-*(-*))~*~*~*~[a]
  */

}

template void validate_basic<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void validate_basic<RDKit::ROMol>(RDKit::ROMol &mol);
