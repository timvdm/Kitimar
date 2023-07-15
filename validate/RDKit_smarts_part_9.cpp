#include "Validate.hpp"

void RDKit_smarts_part_9(OpenBabel::OBMol &mol)
{
    // SMARTS 401 - 428
    validate_contains<"c[OH]">(mol);
    validate_contains<"[PX3](=O)(~O)~[OH]">(mol);
    validate_contains<"[PX3](=O)(~O)~O-[#6]">(mol);
    validate_contains<"[PX4](=O)(~O)(~O)~[OH]">(mol);
    validate_contains<"[PX4](=O)(~O)(~O)~O-[#6]">(mol);
    validate_contains<"[PX4](=O)(~O)~O">(mol);
    validate_contains<"O=P(~O)(~O)(~O)">(mol);
    validate_contains<"[CX4][NH2]">(mol);
    validate_contains<"[c,C]1(~[O;D1])~*!-*~[c,C](~[O;D1])~*!-*~1">(mol);
    //validate_contains<"[$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]">(mol); // FIXME: recursive
    validate_contains<"C14~*~*~*~*~C~1~*~*~C2~C3~*~*~*~C~3~*~*~C~2~4">(mol);
    validate_contains<"[#6][SD3](~O)[#6]">(mol);
    validate_contains<"[#6][SD4](~O)(~O)[#6]">(mol);
    validate_contains<"[#6][SD4](~O)(~O)N">(mol);
    validate_contains<"[#6]S(~O)(~O)[OH]">(mol);
    validate_contains<"[#6]S(~O)(~O)O[#6]">(mol);
    validate_contains<"[#6][SD4](~O)(~O)[Cl,Br]">(mol);
    validate_contains<"[ND3]([CX4])([CX4])[CX4]">(mol);
    validate_contains<"s1c(=N)nnc1[S,N]">(mol);
    //validate_contains<"S=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#6]);!$(C=[O,S,N])]C(=S)O[$([#6]);!$(C=[O,S,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16;D2]([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$([#16;D2]1[CX4][CX4]1);!$(s1aaaa1)]">(mol); // FIXME: recursive
    validate_contains<"[CX4][SH]">(mol);
    validate_contains<"c[SH]">(mol);
    validate_contains<"NC(=S)N">(mol);
    //validate_contains<"[$([CX4,c][CH]=S),$([CX4,c]C(=S)[CX4,c])]">(mol); // FIXME: recursive
    validate_contains<"[D2R0]-[D2R0]-[D2R0]-[D2R0]-[D2R0]">(mol);
    validate_contains<"NC(=O)N">(mol);
}
