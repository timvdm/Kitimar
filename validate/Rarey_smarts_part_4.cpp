#include "Validate.hpp"

void Rarey_smarts_part_4(OpenBabel::OBMol &mol)
{
    // SMARTS 151 - 200
    validate_contains<"O-@S">(mol);
    validate_contains<"O1CCCCC1C2CCCO2">(mol);
    validate_contains<"O1CCCCC1OC2CCC3CCCCC3C2">(mol);
    validate_contains<"O=!@P">(mol);
    validate_contains<"O=!@S">(mol);
    validate_contains<"O=C1CCCC(N1)=O">(mol);
    validate_contains<"O=C1NCC2CCCCC21">(mol);
    validate_contains<"O=C1[#6]~[#6]C(=O)[#6]~[#6]1">(mol);
    validate_contains<"O=CN=[N+]=[N-]">(mol);
    validate_contains<"O=[C,N]">(mol);
    validate_contains<"OO">(mol);
    validate_contains<"OS(=O)(=O)C(F)(F)F">(mol);
    validate_contains<"P(=O)([OH])OP(=O)[OH]">(mol);
    validate_contains<"P(=O)[O,S]">(mol);
    validate_contains<"P(=S)(S)S">(mol);
    validate_contains<"P(OCC)(OCC)(=O)C#N">(mol);
    validate_contains<"P-!@O">(mol);
    validate_contains<"P-!@P">(mol);
    validate_contains<"P-!@[Br]">(mol);
    validate_contains<"P-!@[Cl]">(mol);
    validate_contains<"P-!@[F]">(mol);
    validate_contains<"P-!@[I]">(mol);
    validate_contains<"P-@O">(mol);
    validate_contains<"P-@P">(mol);
    validate_contains<"P=!@P">(mol);
    validate_contains<"P=@P">(mol);
    validate_contains<"P=P@[Br]">(mol);
    validate_contains<"P=P@[Cl]">(mol);
    validate_contains<"P=P@[F]">(mol);
    validate_contains<"P=P@[I]">(mol);
    validate_contains<"P[OH1]">(mol);
    validate_contains<"S(=O)(=O)C#N">(mol);
    validate_contains<"S(=O)(=O)[F,Cl,Br,I]">(mol);
    validate_contains<"S([#6])[CX3](=O)[#6]">(mol);
    validate_contains<"S-!@P">(mol);
    validate_contains<"S-!@S">(mol);
    validate_contains<"S-!@[Br]">(mol);
    validate_contains<"S-!@[Cl]">(mol);
    validate_contains<"S-!@[F]">(mol);
    validate_contains<"S-!@[I]">(mol);
    validate_contains<"S-@P">(mol);
    validate_contains<"S-@S">(mol);
    validate_contains<"S=!@P">(mol);
    validate_contains<"S=!@S">(mol);
    validate_contains<"S=@S">(mol);
    validate_contains<"SC#N">(mol);
    validate_contains<"SS">(mol);
    validate_contains<"[!#1&!#6&!#7&!#8&!#9&!#15&!#16&!#17&!#35&!#53]">(mol);
    validate_contains<"[!#1&!#6&!#7&!#8&!#9&!#16&!#17]">(mol);
    validate_contains<"[!#1;!#2;!#3;!#5;!#6;!#7;!#8;!#9;!#11;!#12;!#15;!#16;!#17;!#19;!#20;!#35;!#53]">(mol);
}