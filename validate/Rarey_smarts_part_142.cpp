#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_142(Mol &mol)
{
    // SMARTS 706 - 710
    validate<"[CD3H0](-[NX3])(-[NX3])=[SD1H0]">(mol);
    validate<"[CD3H0](-[NX3])(=[NX2])-[SD1H1]">(mol);
    validate<"[CD3H0](=[OD1H0])(-[Cl,Br,F])">(mol);
    validate<"[CD3H0](=[OD1H0])-[OD2H0]-[CD3H0]=[OD1H0]">(mol);
    validate<"[CD3H0]1(=[OD1H0])[NX3][CD3H0](=[OD1H0])[NX3][CD3H0](=[OD1H0])[CX4]1">(mol);
}

template void Rarey_smarts_part_142<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_142<RDKit::ROMol>(RDKit::ROMol &mol);
