#include "Validate.hpp"

void RDMACCS_part_4(OpenBabel::OBMol &mol)
{
    // SMARTS 151 - 162
    validate_contains<"[#6]=[#8]">(mol);
    validate_contains<"[!#1]!@[CH2]!@[!#1]">(mol);
    validate_contains<"[#7]~[!#1](~[!#1])~[!#1]">(mol);
    validate_contains<"[#6]-[#8]">(mol);
    validate_contains<"[#6]-[#7]">(mol);
    validate_contains<"[#8]">(mol);
    validate_contains<"[C;H3,H4]">(mol);
    validate_contains<"[#7]">(mol);
    validate_contains<"a">(mol);
    validate_contains<"*~1~*~*~*~*~*~1">(mol);
    validate_contains<"[#8]">(mol);
    validate_contains<"[R]">(mol);
}
