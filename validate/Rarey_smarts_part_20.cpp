#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_20(Mol &mol)
{
    // SMARTS 96 - 100
    validate<"COS(=O)(=O)[C,c]">(mol);
    validate<"COS(=O)O[C,c]">(mol);
    //validate<"C[C@?H](Cl)Br">(mol); // FIXME: stereo
    //validate<"C[C@?](F)(Cl)Br">(mol); // FIXME: stereo
    validate<"C[O,S;R0][C;R0](=S)">(mol);
}

template void Rarey_smarts_part_20<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_20<RDKit::ROMol>(RDKit::ROMol &mol);
