#include "Validate.hpp"

void RDKit_smarts_part_7(OpenBabel::OBMol &mol)
{
    // SMARTS 301 - 350
    validate_contains<"[S,s;D2]C[S,s;D2]">(mol);
    //validate_contains<"[$([S,s]~[S,s]~[C,c]=S),$([S,s]~[C,c](=S)~[S,s,N]),$([S;D2;R0]-S~O)]">(mol); // FIXME: recursive
    validate_contains<"[n+]-C">(mol);
    validate_contains<"[NH]=C([NH2])c">(mol);
    validate_contains<"O=CN[OH]">(mol);
    validate_contains<"[NH;R0][NH;R0]">(mol);
    //validate_contains<"[$(O=C[CH](C=O)C=O),$(N#C[CH](-C=O)-C=O)]">(mol); // FIXME: recursive
    validate_contains<"P(=O)(O[H,C])O[H,C]">(mol);
    //validate_contains<"[$(N#C-C=[CH][C,c]),$([CH](=[C;R0]-[CH]=O)),$([CH](=[C,R]-C(=O)-C));!$([CH]1=CC(=O)C=CC1=*);!$([CH]1=CC(=O)C(=[N,O])C=C1);!$([CH](=C-C=O)-C=O)]">(mol); // FIXME: recursive
    //validate_contains<"[$(N#C-C#C[C,c]),$(C#C-[CH]=O),$(C(#C-C(=O)-[C,c]))]">(mol); // FIXME: recursive
    //validate_contains<"[$(N#CSc1sc(nc1)N),$([S,Se]1C(N)C(=O)[#6][#6]1)]">(mol); // FIXME: recursive
    //validate_contains<"[!$(O=[C,S])][N;R0]=[C;R0]([C,c])[C,c]">(mol); // FIXME: recursive
    validate_contains<"O=CC([Cl,Br,I,F])([Cl,Br,I,F])[Cl,Br,I,F]">(mol);
    validate_contains<"[#6]SC(=O)">(mol);
    //validate_contains<"[$(O([CX4,c])!@[CH,CH2]!@O[CX4,c])]">(mol); // FIXME: recursive
    validate_contains<"C(=O)Cl">(mol);
    //validate_contains<"[$(C#N),$([C,N,S]=O)][CH2,CH][$([C,N,S]=O),$(C#N)]">(mol); // FIXME: recursive
    validate_contains<"C(=O)[Cl,Br]">(mol);
    validate_contains<"[CX4][OH]">(mol);
    validate_contains<"[CX4;H2][OH]">(mol);
    validate_contains<"[CX4;H][OH]">(mol);
    validate_contains<"[CX4;H0][OH]">(mol);
    validate_contains<"[#6][CH]=O">(mol);
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]">(mol); // FIXME: recursive
    validate_contains<"[CX4][Cl,Br,I]">(mol);
    validate_contains<"C#C">(mol);
    validate_contains<"[Cl,Br,I][CH]C=C">(mol);
    validate_contains<"cCC(=O)[OH]">(mol);
    validate_contains<"C([Cl,Br,I])([Cl,Br,I])([Cl,Br,I])C(=O)[OH]">(mol);
    validate_contains<"[#6]C(=O)C(=O)[OH]">(mol);
    validate_contains<"C=CC(=O)[OH]">(mol);
    validate_contains<"N#C-C(=O)[OH]">(mol);
    validate_contains<"N(~O)(~O)-C-C(=O)[OH]">(mol);
    //validate_contains<"[ND3]([CX4,c,H])([CX4,c,H])[CX4][$([CH]),$(C([CX4,c]))]=O">(mol); // FIXME: recursive
    //validate_contains<"[OH][CX4][$([CH]),$(C([CX4,c]))]=O">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-[$([CH]),$(C([CX4,c]))]=O">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-C#N">(mol); // FIXME: recursive
    validate_contains<"[OH][CX4][Cl,Br,I]">(mol);
    //validate_contains<"[Cl,Br,I][$([CX4][CH]=O),$([CX4]C(=O)[CX4,c])]">(mol); // FIXME: recursive
    //validate_contains<"[OH][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol); // FIXME: recursive
    //validate_contains<"O=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2][CX4]),$([$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]),$([ND3]([CX4])([CX4])[CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4]);!$(NC=O)][CX4]C(=O)[OH]">(mol); // FIXME: recursive
    validate_contains<"O(C(=O)[#6])(C(=O)[#6])">(mol);
    validate_contains<"c[NH,NH2]">(mol);
    validate_contains<"c-[ND3]([#6])[#6]">(mol);
    validate_contains<"c">(mol);
    validate_contains<"c-[F,Cl,Br,I]">(mol);
    validate_contains<"c-[Br,I]">(mol);
    //validate_contains<"Fcaa[F,Cl,Br,I,$([C,N,S]=O)]">(mol); // FIXME: recursive
}
