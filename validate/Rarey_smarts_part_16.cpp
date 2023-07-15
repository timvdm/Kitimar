#include "Validate.hpp"

void Rarey_smarts_part_16(OpenBabel::OBMol &mol)
{
    // SMARTS 751 - 800
    validate_contains<"[CX3](=[OX1])C">(mol);
    validate_contains<"[CX3](=[OX1])O">(mol);
    validate_contains<"[CX3](=[OX1])[F,Cl,Br,I]">(mol);
    validate_contains<"[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])">(mol);
    validate_contains<"[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])">(mol);
    validate_contains<"[CX3](=[OX1])[NX3H][CX3](=[OX1])">(mol);
    validate_contains<"[CX3](=[OX1])[OX2][CX3](=[OX1])">(mol);
    validate_contains<"[CX3]=[CD2H1]-[CD2H1]=[O,N,S]">(mol);
    validate_contains<"[CX3]=[OX1]">(mol);
    validate_contains<"[CX4]">(mol);
    validate_contains<"[CX4][CH]=[O,SX2]">(mol);
    validate_contains<"[CX4v4]">(mol);
    validate_contains<"[C](=O)([C,c,O,S])[C,c,O,S]">(mol);
    validate_contains<"[Cl,Br,I][CH2]">(mol);
    validate_contains<"[Cl]-c:2:c:c:1:n:o:n:c:1:c:c:2">(mol);
    validate_contains<"[Cl]C([C&R0])=N">(mol);
    validate_contains<"[F,Cl,Br,I]">(mol);
    validate_contains<"[H+]">(mol);
    validate_contains<"[H,#1]">(mol);
    validate_contains<"[H1,H0-]">(mol);
    validate_contains<"[HC]=O">(mol);
    validate_contains<"[H]">(mol);
    validate_contains<"[H][H]">(mol);
    validate_contains<"[N&D2](=O)">(mol);
    validate_contains<"[N&X4&+,N&X3&+0]">(mol);
    validate_contains<"[N+]#[C-]">(mol);
    validate_contains<"[N,n;+1]">(mol);
    //validate_contains<"[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!D1][!D1])]-;!@[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[N;!D1][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[N;!H0;!$([N;!H0]#C);!$([NH](C(F)(F)(F))S(=O)=O);!$([nH]1nnnc1);!$([nH]1nncn1);!$([nH]1ncnn1)]">(mol); // FIXME: recursive
    //validate_contains<"[N;+0,+1;$(N(=O)~[O;H0;-0,-1])]">(mol); // FIXME: recursive
    validate_contains<"[N;R0]=[N;R0]C#N">(mol);
    validate_contains<"[N;R0]=[N;R0]CC=O">(mol);
    validate_contains<"[N;R0][N;R0]C(=O)">(mol);
    //validate_contains<"[N;X4]-;!@[C;!D4;!D1;!R;!$(C(=O));$(C([N;X4;!D1])[#6;!D1][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[N;X4]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[N;X4]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[N;X4]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[N;X4]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    validate_contains<"[ND1H2]a">(mol);
    validate_contains<"[ND1H2]aaaa[ND1H2]">(mol);
    validate_contains<"[ND2H0]#[CD1H0]">(mol);
    validate_contains<"[ND2H0]=[CD2H0]=[ND2H0]">(mol);
    validate_contains<"[ND2H0]=[CD2H0]=[OD1H0,SD1H0]">(mol);
    validate_contains<"[ND2H0]=[ND1H0]">(mol);
    validate_contains<"[ND2H0]=[ND2H0]">(mol);
    validate_contains<"[ND2H0][OD1H0]">(mol);
    validate_contains<"[ND2H1,ND3H0](a1aaa2aaaaa2a1)a1aaaaa1">(mol);
    validate_contains<"[ND3H0](~[OD1H0])(~[OD1H0])-A">(mol);
    validate_contains<"[ND3H0](~[OD1H0])(~[OD1H0])-a">(mol);
    validate_contains<"[NH,O,S;r3]">(mol);
}
