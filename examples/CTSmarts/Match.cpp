#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecules.hpp>

#include <iostream>

using namespace Kitimar;

int main()
{
    // SMILES: CC(=O)[O-]
    auto molecule = Molecule::mockAcetateAnion();

    bool matches = CTSmarts::match<"C[O-]">(molecule);

    std::cout << "Matching SMARTS \"C[O-]\" in SMILES \"CC(=O)[O-]\": "
              << matches << std::endl;
}
