#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_2(Mol &mol)
{
    // SMARTS 6 - 10
    validate<"OC2=CC(=O)c1c(cccc1)O2">(mol);
    validate<"Nc1nnc(S)s1">(mol);
    validate<"C1C2SCCN2C1">(mol);
    validate<"CP(O)(O)=O">(mol);
    validate<"CCCCCP(O)(O)=O">(mol);
}

template void Agrafiotis_ABCD_part_2<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_2<RDKit::ROMol>(RDKit::ROMol &mol);
