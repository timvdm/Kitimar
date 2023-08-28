#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_234(Mol &mol)
{
    // SMARTS 1166 - 1170
    validate<"c:2:c:c:1:n:c(:c(:n:c:1:c:c:2)-c:3:c:c:c:c:c:3)-c:4:c:c:c:c:c:4-[#8]-[#1]">(mol);
    validate<"c:2:c:c:1:n:c:3:c(:n:c:1:c:c:2):c:c:c:4:c:3:c:c:c:c:4">(mol);
    validate<"c:2:c:c:1:n:n:c(:n:c:1:c:c:2)-[#6](-[#1])(-[#1])-[#6]=[#8]">(mol);
    validate<"c:2:c:c:c:1:c(:c:c:c:1):c:c:2">(mol);
    validate<"cC[N+]">(mol);
}

template void Rarey_smarts_part_234<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_234<RDKit::ROMol>(RDKit::ROMol &mol);
