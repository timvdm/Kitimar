#include "Validate.hpp"

void Rarey_smarts_part_19(OpenBabel::OBMol &mol)
{
    // SMARTS 901 - 950
    validate_contains<"[SD2H0,OD2H0]-[SD2H0,OD2H0]">(mol);
    validate_contains<"[SD2H0]-[Cl,Br,F]">(mol);
    validate_contains<"[SD4H0](=[OD1H0])(=[OD1H0])(-[OD2H0]-C)(-[OD2H0]-C)">(mol);
    validate_contains<"[SD4H0](=[OD1H0])(=[OD1H0])(-[OD2H0]-C)-C">(mol);
    validate_contains<"[SD4H0](=[OD1H0])(=[OD1H0])(-a)-a">(mol);
    validate_contains<"[SD4H0](=[OD1H0])(=[OD1H0])[Cl,Br,F]">(mol);
    validate_contains<"[SH]">(mol);
    validate_contains<"[SX2]">(mol);
    validate_contains<"[SX4](C)(C)(=O)=N">(mol);
    validate_contains<"[X4;R2;r4,r5,r6](@[r4,r5,r6])(@[r4,r5,r6])(@[r4,r5,r6])@[r4,r5,r6]">(mol);
    validate_contains<"[a][CH]=O">(mol);
    validate_contains<"[a]~*~*-[CH3]">(mol);
    validate_contains<"[c,n;R][O;!R][C;!R](=O)[#6,#7,#8;!R]">(mol);
    //validate_contains<"[c;$([*Cl]),$([*H1])]1ccc(O)c(C)c1">(mol); // FIXME: recursive
    //validate_contains<"[c;$(c([cH])([cH,n&H0,O]))][F,Cl,Br,I]">(mol); // FIXME: recursive
    //validate_contains<"[c;$(c([n&H0,cH,S,O]))]-;!@[N;X3;!$(N[N,O,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;$(c[n&H0,cH,S,O])]-;!@[c;$(c[n&H0,cH,S,O])]">(mol); // FIXME: recursive
    //validate_contains<"[c;$(c[n&H0,cH,S,O])]-;!@[n]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    validate_contains<"[cD3H0;r5,r6;R2](:[cD3H0;r5,r6;R2])(:[cD3H0;r5,r6;R2])">(mol);
    validate_contains<"[cD3H0;r6;R2][*;r5,r6;R1][cD3H0;r6;R2]">(mol);
    validate_contains<"[cR1]1[cR1][cR1][cR1][cR1][cR1]1">(mol);
    validate_contains<"[nD3H0,R](~[OD1H0])(a)a">(mol);
    //validate_contains<"[nH0;!$(n-C);!$(n(:c)(:c):a)]1ccccc1">(mol); // FIXME: recursive
    validate_contains<"[nH1](n)nn">(mol);
    validate_contains<"[nH1]1C(=O)CC(=O)O1">(mol);
    validate_contains<"[nH1]1ccc(=O)o1">(mol);
    validate_contains<"[nH1]1cnc(n1)C(F)(F)F">(mol);
    validate_contains<"[nH1]1cnnc1C(F)(F)F">(mol);
    validate_contains<"[nH1]1occc1=O">(mol);
    validate_contains<"[nH1]nnn">(mol);
}
