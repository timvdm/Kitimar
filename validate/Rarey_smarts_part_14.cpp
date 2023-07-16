#include "Validate.hpp"

void Rarey_smarts_part_14(OpenBabel::OBMol &mol)
{
    // SMARTS 651 - 700
    //validate_contains<"[C;!D4;!D1;!R;$(C(=O));$(C([#7;R;D3])[!S;!D1][!D1])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!N;!D1][!D1])]-;!@[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!N;!D1][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!S;!D1][!D1])]-;!@[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!S;!D1][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;$(C(=O));$(C([O;D2])[#6;!D1]!=[!D1])]-;!@[O;D2;$(O(C(=O))[#6;!D1][#6;!D1]);!$(OCO);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;$(C(=O));$(C([n])[!S;!D1][!D1])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;$(C(S=O)[#6;!D1][!D1])]-;!@[S;$(S(=O)([C;!D1;!R])[#6;!D1][!D1])]">(mol); // FIXME: recursive
    validate_contains<"[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]">(mol);
    validate_contains<"[C;!r5]([C;!r5])=[C;!r5](C)[C;!r5]=[O,SX2;!r5]">(mol);
    //validate_contains<"[C;$(C(=O))]-;!@[N;!$(N(C=O)(C=O));!$(N[N,O,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;$(C(=O))]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol); // FIXME: recursive
    //validate_contains<"[C;$(C(=O))]-;!@[n]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1])[#6;!D1;R]!=,@[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    validate_contains<"[C;H1](=[O,S])[C,c]">(mol);
    //validate_contains<"[C;H2;$([CH2]([CH2])[CH2])][H]">(mol); // FIXME: recursive
    //validate_contains<"[C;H3;$(Cc([cH])([cH,n&H0,O]))][H]">(mol); // FIXME: recursive
    validate_contains<"[C;X4;!D4]-;!@[C;D3]">(mol);
    validate_contains<"[C;X4;!D4]-;!@[C;X4;!D4;!R]">(mol);
    validate_contains<"[C;X4;!D4]-;!@[N;X3;!D3]">(mol);
    //validate_contains<"[C;X4;!D4]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol); // FIXME: recursive
    //validate_contains<"[C;X4;!D4]-;!@[c;$(c[n&H0,cH,S,O])]">(mol); // FIXME: recursive
    validate_contains<"[C;X4;!D4]-;!@[n]">(mol);
    //validate_contains<"[CD1H2,CD2H1,CD3H0;!$(C=,#*)]">(mol); // FIXME: recursive
    validate_contains<"[CD1H2]=[CD2H1]-[OD2H0,SD2H0]-C">(mol);
    validate_contains<"[CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0]">(mol);
    validate_contains<"[CD2H1]=[CD3H0]-[ND1H2]">(mol);
}