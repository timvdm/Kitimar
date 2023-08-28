#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_3(Mol &mol)
{
    // SMARTS 11 - 15
    validate<"N2CCC13CCCCC1C2Cc4c3cccc4">(mol);
    validate<"s1cncc1">(mol);
    validate<"C34CCC1C(CCC2CC(=O)CCC12)C3CCC4">(mol);
    validate<"CCCCCCCCCCCP(O)(O)=O">(mol);
    validate<"CC1CCCC1">(mol);
}

template void Agrafiotis_ABCD_part_3<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_3<RDKit::ROMol>(RDKit::ROMol &mol);
