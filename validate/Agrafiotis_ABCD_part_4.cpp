#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_4(Mol &mol)
{
    // SMARTS 16 - 20
    validate<"CCC1CCCC1">(mol);
    validate<"CCCC1CCCC1">(mol);
    validate<"[*D1]=[#6D2]-[#6D2]-[#6D2]=[#8D1]">(mol);
    validate<"[*D1]=[#6]-[#6]-[#6]=[#8D1]">(mol);
    validate<"[#6]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D3](=[#8D1])-[#8]">(mol);
}

template void Agrafiotis_ABCD_part_4<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_4<RDKit::ROMol>(RDKit::ROMol &mol);
