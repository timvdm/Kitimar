#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecule.hpp>

#include <iostream>

// SMILES: CC(=O)[O-]
// Generated by SmiToMock
auto molecule = Kitimar::Molecule::MockMolecule{
    {{6,0,0,1,3,false,false},{6,0,0,3,0,false,false},{8,0,0,1,0,false,false},{8,0,-1,1,0,false,false}},
    {{0,1,1,false,false},{1,2,2,false,false},{1,3,1,false,false}}
};

using namespace Kitimar;

int main()
{
    bool matches = CTSmarts::contains<"C[O-]">(molecule);

    std::cout << "Matching SMARTS \"C[O-]\" in SMILES \"CC(=O)[O-]\": "
              << matches << std::endl;
}
