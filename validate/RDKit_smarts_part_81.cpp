#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_81(Mol &mol)
{
    // SMARTS 401 - 405
    validate<"c[OH]">(mol);
    validate<"[PX3](=O)(~O)~[OH]">(mol);
    validate<"[PX3](=O)(~O)~O-[#6]">(mol);
    validate<"[PX4](=O)(~O)(~O)~[OH]">(mol);
    validate<"[PX4](=O)(~O)(~O)~O-[#6]">(mol);
}

template void RDKit_smarts_part_81<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_81<RDKit::ROMol>(RDKit::ROMol &mol);
