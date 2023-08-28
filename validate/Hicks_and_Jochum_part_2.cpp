#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Hicks_and_Jochum_part_2(Mol &mol)
{
    // SMARTS 6 - 10
    validate<"[#6]-[#6]-1-2-[#6]-[#6]-[#6]-[#6]-1-[#6]-3-[#6]-[#6]-[#6]-4-[#6]-[#6]-[#6]-[#6]-[#6]-4(-[#6])-[#6]-3-[#6]-[#6]-2">(mol);
    validate<"[#6D2]-1=[#6D2]-[#6D3]-2-[#6]-[#6]-[#6D3]-1-[#6D2]-2">(mol);
    validate<"[#6]-1=[#6]-[#6]-2-[#6]-[#6]-[#6]-1-[#6]-2">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]1=,:[#6D2][#6D2]=,:[#6D2][#6D2]=,:[#6D3]1-[*;D1,D2,D3,D4]">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]1=,:[#6D2][#6D2]=,:[#6D2][#6D3](-[*;D1,D2,D3,D4])=,:[#6D2]1">(mol);
}

template void Hicks_and_Jochum_part_2<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Hicks_and_Jochum_part_2<RDKit::ROMol>(RDKit::ROMol &mol);
