#include "Validate.hpp"

void Rarey_smarts_part_13(OpenBabel::OBMol &mol)
{
    // SMARTS 601 - 650
    //validate_contains<"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3]C=N)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3]N=C)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3]N=N)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX4+]),$([NX4]=*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([OH]-*=[!#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([OX1]=[CX3])]">(mol); // FIXME: recursive
    //validate_contains<"[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([OX2])]">(mol); // FIXME: recursive
    //validate_contains<"[$([OX2]C=N)]">(mol); // FIXME: recursive
    //validate_contains<"[$([SX1]=[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([SX1]~P)]">(mol); // FIXME: recursive
    //validate_contains<"[$([SX3]=N)]">(mol); // FIXME: recursive
    //validate_contains<"[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]">(mol); // FIXME: recursive
    //validate_contains<"[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]">(mol); // FIXME: recursive
    //validate_contains<"[$([cX2+](:*):*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([cX3](:*):*),$([cX2+](:*):*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([nX2r5]:[a-]),$([nX2r5]:[a]:[a-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]">(mol); // FIXME: recursive
    //validate_contains<"[$(c1(=O)ccn([C,c])cc1),$(c1(=O)n([C,c])cccc1)]">(mol); // FIXME: recursive
    //validate_contains<"[$(c:cCl),$(c:c:cCl),$(c:c:c:cCl)]-[$(c:cCl),$(c:c:cCl),$(c:c:c:cCl)]">(mol); // FIXME: recursive
    validate_contains<"[*!H0,#1]">(mol);
    validate_contains<"[+,++,+++]">(mol);
    validate_contains<"[+1]~*~*~[-1]">(mol);
    //validate_contains<"[+;!$([+]~[-])]">(mol); // FIXME: recursive
    validate_contains<"[+H]">(mol);
    validate_contains<"[+]">(mol);
    validate_contains<"[-,--,---]">(mol);
    //validate_contains<"[-;!$([-]~[+])]">(mol); // FIXME: recursive
    validate_contains<"[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]">(mol);
    validate_contains<"[Br,Cl,I][CX4,CH,CH2]">(mol);
    validate_contains<"[C+,Cl+,I+,P+,S+]">(mol);
    validate_contains<"[C,N,O,S,a]c1[n,o,s]c2ccccc2[n,o,s]1">(mol);
    validate_contains<"[C,O]=[#6]1[#7,#8,#16][#6](=[O,N,SX1])c2ccccc12">(mol);
    validate_contains<"[C,P;H1](=[O,S])[O,S]">(mol);
    validate_contains<"[C,c](=N)N">(mol);
    validate_contains<"[C,c](=N)N[C,S](=O)">(mol);
    //validate_contains<"[C;!D4;!D1;!R;!$(C(=O));$(C([N;X3;!D1])[#6;!D1][!D1])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([C;!D1;!R])[#6;!D1]!=[!D1])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C(=O));$(C([O;D2])[#6;!D1]~[!D1])]-;!@[O;D2;$(O([C;!D1;!R])[#6;!D1][!D1]);!$(OCO);!$(OC=*);!$(O[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C(=S));$(C([S;D2])[#6;!D1]~[!D1])]-;!@[S;D2;$(S([C;!D1;!R])[#6;!D1][!D1]);!$(SCS);!$(SC=*);!$(S[P,S])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([#7;R;D3])!=[#6;!D1]~[!D1])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([n])!=[#6;!D1]~[!D1])]-;!@[n;!$(n@[#6](=O))]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
    //validate_contains<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol); // FIXME: recursive
}
