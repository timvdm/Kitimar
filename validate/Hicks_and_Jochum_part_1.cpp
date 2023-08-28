#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Hicks_and_Jochum_part_1(Mol &mol)
{
    // SMARTS 1 - 5
    validate<"[*D1]=[#6D2]-[#6D2]-[#6D2]=[#8D1]">(mol);
    validate<"[*D1]=[#6]-[#6]-[#6]=[#8D1]">(mol);
    validate<"[#6]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D2]-[#6D3](=[#8D1])-[#8]">(mol);
    validate<"[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8D1])-[#8]">(mol);
    validate<"[#6D1]-[#6]-1-2-[#6]-[#6D2]-[#6D2]-[#6D3]-1-[#6D3]-3-[#6D2]-[#6D2]-[#6D3]-4-[#6D2]-[#6]-[#6D2]-[#6D2]-[#6]-4(-[#6D1])-[#6D3]-3-[#6D2]-[#6D2]-2">(mol);
}

template void Hicks_and_Jochum_part_1<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Hicks_and_Jochum_part_1<RDKit::ROMol>(RDKit::ROMol &mol);
