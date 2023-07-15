#include "Validate.hpp"

void RDKit_smarts_part_8(OpenBabel::OBMol &mol)
{
    // SMARTS 351 - 400
    validate_contains<"cB([OH])[OH]">(mol);
    validate_contains<"[N;R0;D2]~[N;R0]~[N;R0;D1]">(mol);
    //validate_contains<"[N;D2]([C,c;!$(C=[O,S,N])])=[N;D2]-[C,c;!$(C=[O,S,N])]">(mol); // FIXME: recursive
    validate_contains<"[OH][CX4][CX4][Cl,Br,I]">(mol);
    //validate_contains<"[OH][CX4][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C[CX4]C(=O)[OH]">(mol); // FIXME: recursive
    validate_contains<"[CX4]([OH])[CX4]C(=O)[OH]">(mol);
    validate_contains<"[CX4]([Cl,Br,I])[CX4]C(=O)[OH]">(mol);
    validate_contains<"[#6]C(=O)[CX4]C(=O)[OH]">(mol);
    //validate_contains<"[#6]C(=O)[CX4]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)[OH]">(mol); // FIXME: recursive
    validate_contains<"C(=O)[OH]">(mol);
    validate_contains<"O([CX4,c])C(=O)O[CX4,c]">(mol);
    //validate_contains<"[$([CX4,c][CH]=O),$([CX4,c]C(=O)[CX4,c])]">(mol); // FIXME: recursive
    validate_contains<"ClC(=O)O[CX4,c]">(mol);
    validate_contains<"[CX4,c]-C#N">(mol);
    validate_contains<"N[CX4]C(=O)N[CX4]C(=O)">(mol);
    validate_contains<"[S;D2]-[S;D2]">(mol);
    //validate_contains<"[$([S;D2]([CX4,c])!@[CH,CH2]!@[S;D2][CX4,c])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][$([NH2]),$([NH][CX4]),$([N;R0]([CX4])[CX4])]);!$(CC=[O,S,N]);!$(C[O,S])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[O,S,N])]=[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[N,S])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][OH]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]O[CX4]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol); // FIXME: recursive
    validate_contains<"O1[CX4][CX4]1">(mol);
    //validate_contains<"[$([#6]);!$(C=[O,S,N])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$(O([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$(O1[CX4][CX4]1)]">(mol); // FIXME: recursive
    validate_contains<"[OH][CH,CH2]O[CX4,c]">(mol);
    validate_contains<"O([#6])-C([#6])([#6])-[OH]">(mol);
    //validate_contains<"[OH][$([NX3]([C;!$(C=[O,S,N])])[C;!$(C=[O,S,N])]),$([NH][CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH;R0]([C;!$(C=[O,S,N])]))][$([NH;R0][C;!$(C=[O,S,N])])]">(mol); // FIXME: recursive
    validate_contains<"C=N[NH2]">(mol);
    validate_contains<"C=NC=O">(mol);
    //validate_contains<"[$([C;R0]=[N;R0]);!$(C(~[N,O,S])(~[N,O,S]));!$([C;R0]=[N;R0]~[N,O,n])]">(mol); // FIXME: recursive
    validate_contains<"C#N-[#6]">(mol);
    validate_contains<"O=C=N-[#6]">(mol);
    validate_contains<"S=C=N-[#6]">(mol);
    validate_contains<"O([CX4,c])-C([CX4,c])([CX4,c])-O([CX4,c])">(mol);
    validate_contains<"[CX4,c]C(=O)[CX4,c]">(mol);
    validate_contains<"[C;R0;X4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[C;R0;X4]">(mol);
    //validate_contains<"[$([C;R1]);!$(C(N)N)](=O)@[$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([C;R1]);!$(C(O)N);!$(C(O)O)](=O)@[$(O);!$(O(C=O))]">(mol); // FIXME: recursive
    validate_contains<"OC(=O)CC(=O)[OH]">(mol);
    validate_contains<"[r8,r9,r10,r11,r12,r13,r14]">(mol);
    validate_contains<"*-C#N">(mol);
    validate_contains<"*-N(=O)(~O)">(mol);
    validate_contains<"[N;D4]">(mol);
    validate_contains<"OC(=O)C(=O)O">(mol);
    validate_contains<"C=[N;R0]-[OH]">(mol);
    validate_contains<"O~O">(mol);
}
