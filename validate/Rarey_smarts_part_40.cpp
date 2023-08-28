#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_40(Mol &mol)
{
    // SMARTS 196 - 200
    validate<"SC#N">(mol);
    validate<"SS">(mol);
    validate<"[!#1&!#6&!#7&!#8&!#9&!#15&!#16&!#17&!#35&!#53]">(mol);
    validate<"[!#1&!#6&!#7&!#8&!#9&!#16&!#17]">(mol);
    validate<"[!#1;!#2;!#3;!#5;!#6;!#7;!#8;!#9;!#11;!#12;!#15;!#16;!#17;!#19;!#20;!#35;!#53]">(mol);
}

template void Rarey_smarts_part_40<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_40<RDKit::ROMol>(RDKit::ROMol &mol);
