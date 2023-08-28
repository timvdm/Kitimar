#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_163(Mol &mol)
{
    // SMARTS 811 - 815
    validate<"[NX2]=N">(mol);
    validate<"[NX2]=[NX2]">(mol);
    validate<"[NX2]=[OX1]">(mol);
    validate<"[NX3+]=[CX3]">(mol);
    validate<"[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]">(mol);
}

template void Rarey_smarts_part_163<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_163<RDKit::ROMol>(RDKit::ROMol &mol);
