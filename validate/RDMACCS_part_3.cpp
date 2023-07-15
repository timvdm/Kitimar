#include "Validate.hpp"

void RDMACCS_part_3(OpenBabel::OBMol &mol)
{
    // SMARTS 101 - 150
    validate_contains<"Cl">(mol);
    validate_contains<"[!#6;!#1;!H0]~*~[CH2]~[!#1]">(mol);
    validate_contains<"*@*(@*)@*">(mol);
    validate_contains<"[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]">(mol);
    validate_contains<"[F,Cl,Br,I]~*(~[!#1])~[!#1]">(mol);
    validate_contains<"[CH3]~*~*~*~[CH2]~[!#1]">(mol);
    validate_contains<"[!#1]~[CH2]~[#8]">(mol);
    validate_contains<"[#7]~[#6]~[#8]">(mol);
    validate_contains<"[#7]~*~[CH2]~[!#1]">(mol);
    validate_contains<"[!#1]~*(~[!#1])(~[!#1])~[!#1]">(mol);
    validate_contains<"[#8]!:*:*">(mol);
    validate_contains<"[CH3]~[CH2]~[!#1]">(mol);
    validate_contains<"[CH3]~*~[CH2]~[!#1]">(mol);
    //validate_contains<"[$([CH3]~*~*~[CH2]~[!#1]),$([CH3]~*1~*~[CH2]1)]">(mol); // FIXME: recursive
    validate_contains<"[#7]~*~[#8]">(mol);
    //validate_contains<"[$([!#1]~[CH2]~[CH2]~[!#1]),$(*1~[CH2]~[CH2]1)]">(mol); // FIXME: recursive
    validate_contains<"[#7]=*">(mol);
    validate_contains<"[!#6;R]">(mol);
    validate_contains<"[#7;R]">(mol);
    validate_contains<"[!#1]~[#7](~[!#1])~[!#1]">(mol);
    validate_contains<"[#8]~[#6]~[#8]">(mol);
    validate_contains<"[!#6;!#1]~[!#6;!#1]">(mol);
    validate_contains<"[!#1]!@[#8]!@[!#1]">(mol);
    validate_contains<"*@*!@[#8]">(mol);
    //validate_contains<"[$([!#1]~[CH2]~*~*~*~[CH2]~[!#1]),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~*~[R]1@[R]@[CH2;R]1)]">(mol); // FIXME: recursive
    //validate_contains<"[$([!#1]~[CH2]~*~*~[CH2]~[!#1]),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~[R]1@[R]@[CH2;R]1)]">(mol); // FIXME: recursive
    validate_contains<"[!#6;!#1]~[!#6;!#1]">(mol);
    validate_contains<"[!#6;!#1;!H0]">(mol);
    validate_contains<"[#8]~*~[CH2]~[!#1]">(mol);
    validate_contains<"*@*!@[#7]">(mol);
    validate_contains<"[F,Cl,Br,I]">(mol);
    validate_contains<"[#7]!:*:*">(mol);
    validate_contains<"[#8]=*">(mol);
    validate_contains<"[!#6;R]">(mol);
    validate_contains<"[!#6;!#1]~[CH2]~[!#1]">(mol);
    validate_contains<"[O;!H0]">(mol);
    validate_contains<"[#8]">(mol);
    validate_contains<"[CH3]">(mol);
    validate_contains<"[#7]">(mol);
    validate_contains<"*@*!@[#8]">(mol);
    validate_contains<"[!#1]!:*:*!:[!#1]">(mol);
    validate_contains<"*~1~*~*~*~*~*~1">(mol);
    validate_contains<"[#8]">(mol);
    //validate_contains<"[$([!#1]~[CH2]~[CH2]~[!#1]),$([R]1@[CH2;R]@[CH2;R]1)]">(mol); // FIXME: recursive
    validate_contains<"[!#1]~[!#6;!#1](~[!#1])~[!#1]">(mol);
    validate_contains<"[C;H3,H4]">(mol);
    validate_contains<"[!#1]!@*@*!@[!#1]">(mol);
    validate_contains<"[#7;!H0]">(mol);
    validate_contains<"[#8]~[#6](~[#6])~[#6]">(mol);
    validate_contains<"[!#6;!#1]~[CH2]~[!#1]">(mol);
}
