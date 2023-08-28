#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_1(Mol &mol)
{
    // SMARTS 1 - 5
    validate<"ONC1CC(C(O)C1O)[n]2cnc3c(NC4CC4)ncnc23">(mol);
    validate<"Nc1ncnc2[n]cnc12">(mol);
    validate<"CNc1ncnc2[n](C)cnc12">(mol);
    validate<"Nc1ncnc2[n](cnc12)C3CCCC3">(mol);
    validate<"CC12CCC3C(CCC4=CC(O)CCC34C)C1CCC2">(mol);
}

template void Agrafiotis_ABCD_part_1<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_1<RDKit::ROMol>(RDKit::ROMol &mol);
