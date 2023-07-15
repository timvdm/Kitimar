#include "Validate.hpp"

void RDMACCS_part_2(OpenBabel::OBMol &mol)
{
    // SMARTS 51 - 100
    validate_contains<"[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]">(mol);
    validate_contains<"[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]">(mol);
    validate_contains<"[#8]~[#16]~[#8]">(mol);
    validate_contains<"[#8]~[#7](~[#8])~[#6]">(mol);
    validate_contains<"[#8R]">(mol);
    validate_contains<"[!#6;!#1]~[#16]~[!#6;!#1]">(mol);
    validate_contains<"[#16]!:*:*">(mol);
    validate_contains<"[#16]=[#8]">(mol);
    validate_contains<"[!#1]~[#16](~[!#1])~[!#1]">(mol);
    validate_contains<"*@*!@*@*">(mol);
    validate_contains<"[#7]=[#8]">(mol);
    validate_contains<"*@*!@[#16]">(mol);
    validate_contains<"c:n">(mol);
    validate_contains<"[#6]~[#6](~[#6])(~[#6])~[!#1]">(mol);
    validate_contains<"[!#6;!#1]~[#16]">(mol);
    validate_contains<"[!#6;!#1;!H0]~[!#6;!#1;!H0]">(mol);
    validate_contains<"[!#6;!#1]~[!#6;!#1;!H0]">(mol);
    validate_contains<"[!#6;!#1]~[#7]~[!#6;!#1]">(mol);
    validate_contains<"[#7]~[#8]">(mol);
    validate_contains<"[#8]~*~*~[#8]">(mol);
    validate_contains<"[#16]=*">(mol);
    validate_contains<"[CH3]~*~[CH3]">(mol);
    validate_contains<"[!#1]!@[#7]@[!#1]">(mol);
    validate_contains<"[#6]=[#6](~[!#1])~[!#1]">(mol);
    validate_contains<"[#7]~*~[#7]">(mol);
    validate_contains<"[#6]=[#7]">(mol);
    validate_contains<"[#7]~*~*~[#7]">(mol);
    validate_contains<"[#7]~*~*~*~[#7]">(mol);
    validate_contains<"[#16]~*(~[!#1])~[!#1]">(mol);
    validate_contains<"[!#1]~[CH2]~[!#6;!#1;!H0]">(mol);
    validate_contains<"[!#6;!#1]~1~*~*~*~*~1">(mol);
    validate_contains<"[NH2]">(mol);
    validate_contains<"[#6]~[#7](~[#6])~[#6]">(mol);
    validate_contains<"[C;H2,H3][!#6;!#1][C;H2,H3]">(mol);
    validate_contains<"[F,Cl,Br,I]!@*@*">(mol);
    validate_contains<"[#16]">(mol);
    validate_contains<"[#8]~*~*~*~[#8]">(mol);
    //validate_contains<"[$([!#6;!#1;!H0]~*~*~[CH2]~[!#1]),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]">(mol); // FIXME: recursive
    //validate_contains<"[$([!#6;!#1;!H0]~*~*~*~[CH2]~[!#1]),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]">(mol); // FIXME: recursive
    validate_contains<"[#8]~[#6](~[#7])~[#6]">(mol);
    validate_contains<"[!#6;!#1]~[CH3]">(mol);
    validate_contains<"[!#6;!#1]~[#7]">(mol);
    validate_contains<"[#7]~*~*~[#8]">(mol);
    validate_contains<"*~1~*~*~*~*~1">(mol);
    validate_contains<"[#7]~*~*~*~[#8]">(mol);
    validate_contains<"[!#6;!#1]~1~*~*~*~*~*~1">(mol);
    validate_contains<"[#6]=[#6]">(mol);
    validate_contains<"[!#1]~[CH2]~[#7]">(mol);
    //validate_contains<"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]">(mol); // FIXME: recursive
    validate_contains<"[!#6;!#1]~[#8]">(mol);
}
