#include "Validate.hpp"

void Rarey_smarts_part_18(OpenBabel::OBMol &mol)
{
    // SMARTS 851 - 900
    validate_contains<"[OH1]C1=NC(=O)CC1=O">(mol);
    validate_contains<"[OH1]C1=NC(=O)NO1">(mol);
    validate_contains<"[OH1]C1=NC(=O)ON1">(mol);
    validate_contains<"[OH1]C1=NC(=O)c2ccccc21">(mol);
    validate_contains<"[OH1]C1=NC=NO1">(mol);
    validate_contains<"[OH1]C1=NN=CO1">(mol);
    validate_contains<"[OH1]C1=NOC=C1">(mol);
    validate_contains<"[OH1]C1=NS(=O)(=O)c2ccccc21">(mol);
    validate_contains<"[OH1]C1=NSN=C1">(mol);
    validate_contains<"[OH1]C1=N[NH1]C=N1">(mol);
    validate_contains<"[OH1]C1NC(=O)C(=O)C1">(mol);
    validate_contains<"[OH1]C=C">(mol);
    validate_contains<"[OH1]NC(=O)">(mol);
    validate_contains<"[OH1]NC=[O,S]">(mol);
    //validate_contains<"[OH1][C,c,N;!$(C=O)]">(mol); // FIXME: recursive
    validate_contains<"[OH1][P,C,S](=O)">(mol);
    validate_contains<"[OH1]c1c[c,n]ccc1">(mol);
    validate_contains<"[OH1]c1oncc1">(mol);
    validate_contains<"[OH]c1cccc2ccccc12">(mol);
    validate_contains<"[OX1]=CN">(mol);
    validate_contains<"[OX2,OX1-][OX2,OX1-]">(mol);
    validate_contains<"[OX2H+]=*">(mol);
    validate_contains<"[OX2H,OX1H0-]">(mol);
    validate_contains<"[OX2H0]">(mol);
    validate_contains<"[OX2H]">(mol);
    validate_contains<"[OX2H]P">(mol);
    validate_contains<"[OX2H][#6X3]=[#6]">(mol);
    //validate_contains<"[OX2H][$(C=C),$(cc)]">(mol); // FIXME: recursive
    validate_contains<"[OX2H][CX3]=[OX1]">(mol);
    validate_contains<"[OX2H][cX3]:[c]">(mol);
    validate_contains<"[OX3H2+]">(mol);
    //validate_contains<"[P,S]-;!@[O;D2;$(O([P,S])[#6;!D1][!D1]);!$(OCO);!$(OC=*)]">(mol); // FIXME: recursive
    //validate_contains<"[P,S]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol); // FIXME: recursive
    validate_contains<"[P,S][Cl,Br,F,I]">(mol);
    validate_contains<"[R0;D2][R0;D2][R0;D2][R0;D2]">(mol);
    validate_contains<"[R0;D2]~[R0;D2]~[R0;D2]~[R0;D2]">(mol);
    validate_contains<"[R0]">(mol);
    validate_contains<"[R]">(mol);
    validate_contains<"[R](-*(-*))~*~*~*~[a]">(mol);
    validate_contains<"[R](-*(-*))~*~*~[a]">(mol);
    validate_contains<"[R](-*(-*))~*~[a]">(mol);
    validate_contains<"[R]~*-[CH3]">(mol);
    validate_contains<"[R]~*~*-[CH3]">(mol);
    validate_contains<"[R]~[D3]">(mol);
    validate_contains<"[S,C](=[O,S])[F,Br,Cl,I]">(mol);
    validate_contains<"[S-][CX3](=S)[#6]">(mol);
    validate_contains<"[S;!H0]">(mol);
    //validate_contains<"[S;X4;$(S(=O)(=O))]-;!@[#7;X3;!D1;$([#7](S(=O)(=O))[#6;!D1][!D1]);!$([#7]C(=O))]">(mol); // FIXME: recursive
    //validate_contains<"[S;X4;$(S(=O)(=O))]-;!@[N;X3;!$(N[N,O,S])]">(mol); // FIXME: recursive
    validate_contains<"[SD1H1]">(mol);
}
