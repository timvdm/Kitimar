#include "Validate.hpp"

void RDKit_smarts_part_6(OpenBabel::OBMol &mol)
{
    // SMARTS 251 - 300
    validate_contains<"C#C-C(=O)-O-C-C">(mol);
    validate_contains<"C#C-[C;!H0](=O)">(mol);
    validate_contains<"C#C-C(=O)[N,!H0]">(mol);
    validate_contains<"C=C-[C;H1](=O)">(mol);
    validate_contains<"[O;D2]-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"[C;!H3]-O-S(=O)(=O)O">(mol);
    validate_contains<"[C;!a]-S(=O)(=O)[O-]">(mol);
    validate_contains<"[N;!a]-S(=O)(=O)[O-]">(mol);
    validate_contains<"C(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"C(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"C(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"C(=O)-[C;H2]([F,Br,I,Cl])">(mol);
    validate_contains<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"S(=O)(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate_contains<"S(=O)(=O)-[C;H2]([F,Br,I,Cl])">(mol);
    //validate_contains<"[$(C-[Cl,Br,I]),$(O=C-[Cl,Br,I,F]),$(O=C([CH,CH2][Cl,Br,I,F])[O,C]),$(C~O~[Cl,Br,I,F][CH,CH2]),$(n1c([Cl,Br,I,F])nccc1);!$(C=C-[Cl,Br,I]);!$(ClC-[Cl,Br,I,F])]">(mol); // FIXME: recursive
    validate_contains<"[O,N,S]1CC1">(mol);
    validate_contains<"O=COC=O">(mol);
    validate_contains<"N=C=O">(mol);
    validate_contains<"S=C=N">(mol);
    validate_contains<"O~O">(mol);
    validate_contains<"[Si]~N">(mol);
    validate_contains<"P[S,N]">(mol);
    validate_contains<"[N;R0]=[N;R0]=[C;R0]">(mol);
    validate_contains<"[N+]#N-*">(mol);
    validate_contains<"[C,c]-[S;D2]-[O,N]">(mol);
    validate_contains<"[Cl,Br,I,F][S,P,Si,N]">(mol);
    validate_contains<"[Cl,Br,I]CC[S,N]">(mol);
    validate_contains<"[Si]-O-*">(mol);
    validate_contains<"S(O[C,c,N,n])(~O)[C,c,N,n]">(mol);
    validate_contains<"[N;R0](~N)~O">(mol);
    validate_contains<"N(~O)(~O)(~O)-*">(mol);
    validate_contains<"[N+]([O-])(=C)-*">(mol);
    //validate_contains<"[!$([C,c]-N(=O)~O);$([!O]~[N;R0]=O)]">(mol); // FIXME: recursive
    validate_contains<"N([H,C;X4])([H,C;X4])-[C;R0;X4]-N[C;X4]([H,C;X4])([H,C;X4])">(mol);
    //validate_contains<"N#C[C;R0;X4]O[!$(O=[C,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;R0](=[C;R0])-[S,O,N;R0][!$(O=[C,S])]">(mol); // FIXME: recursive
    validate_contains<"[O-]-[C;R0]=[C;R0]">(mol);
    validate_contains<"[C;R0](=[C;R0])-[OH]">(mol);
    validate_contains<"[Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]">(mol);
    validate_contains<"O=CON1C(=O)CCC1=O">(mol);
    validate_contains<"O=COn1cncc1">(mol);
    validate_contains<"Fc1c(OC=O)c(F)c(F)c(F)c1F">(mol);
    validate_contains<"[S;D2][C;R0](C)(C)C">(mol);
    validate_contains<"[C;R0;X4]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH][CH2]">(mol);
    validate_contains<"[CH]=O">(mol);
    validate_contains<"[C;R0;X4]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[C;R0;X4]">(mol);
    validate_contains<"[S;D2;R0]-[S;D2]">(mol);
    validate_contains<"[SH]">(mol);
}
