#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Hicks_and_Jochum_part_6(Mol &mol)
{
    // SMARTS 26 - 29
    validate<"[*;D1,D2,D3,D4]-[#6D3]-1=[#6D3](-[*;D1,D2,D3,D4])-[#34D2]-[#6D2]=[#6D2]-1">(mol);
    validate<"[#6]-1=[#6]-[#34D2]-[#6]=[#6]-1">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][*D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
    validate<"[#6]=,:1[#6]=,:[#6][*D2][#6]=,:1">(mol);
}

template void Hicks_and_Jochum_part_6<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Hicks_and_Jochum_part_6<RDKit::ROMol>(RDKit::ROMol &mol);
