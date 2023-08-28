#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_6(Mol &mol)
{
    // SMARTS 26 - 30
    validate<"C!@[N$(C@*)]">(mol);
    validate<"C!@[N$(N@*)]">(mol);
    validate<"C!@[N$(O@*)]">(mol);
    validate<"C!@[N$(S@*)]">(mol);
    validate<"C#!@C">(mol);
}

template void Rarey_smarts_part_6<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_6<RDKit::ROMol>(RDKit::ROMol &mol);
