#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_109(Mol &mol)
{
    // SMARTS 541 - 545
    validate<"[#8]=[#6]-c2c1nc(-[#6](-[#1])-[#1])cc(-[#8]-[#1])n1nc2">(mol);
    validate<"[#8]=[C,N]">(mol);
    validate<"[#8]=[C,N]-aaa[F,Cl,Br,I]">(mol);
    validate<"[#9,#17,#35,#53]">(mol);
    validate<"[#9,#17,#35,#53]!@*@*">(mol);
}

template void Rarey_smarts_part_109<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_109<RDKit::ROMol>(RDKit::ROMol &mol);
