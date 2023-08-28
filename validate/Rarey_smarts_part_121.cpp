#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_121(Mol &mol)
{
    // SMARTS 601 - 605
    validate<"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]">(mol);
    validate<"[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]">(mol);
    validate<"[$([NX3]C=N)]">(mol);
    validate<"[$([NX3]N=C)]">(mol);
    validate<"[$([NX3]N=N)]">(mol);
}

template void Rarey_smarts_part_121<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_121<RDKit::ROMol>(RDKit::ROMol &mol);
