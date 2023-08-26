#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecules.hpp>

#include <iostream>

using namespace Kitimar;

bool isCarbonOfDegree3(auto &mol, auto atom)
{
    // Optimized case: single atom (requires no runtime memory)
    return CTSmarts::match_atom<"[CD3]">(mol, atom);
}

bool isCarboxylCarbon(auto &mol, auto atom)
{
    return CTSmarts::match_atom<"C(=O)O">(mol, atom); // atom matches first SMARTS atom
}

bool isCarboxylDoubleBond(auto &mol, auto bond)
{
    return CTSmarts::match_bond<"O=CO">(mol, bond); // bond matches first SMARTS bond
}

int main()
{
    // SMILES: CC(=O)[O-]
    auto molecule = Molecule::mockAcetateAnion();

    std::cout << "Is carbon of degree 3: "
              << isCarbonOfDegree3(molecule, get_atom(molecule, 1)) << std::endl;

    std::cout << "Is carboxyl carbon: "
              << isCarboxylCarbon(molecule, get_atom(molecule, 1)) << std::endl;

    // Optimized case: single bond (requires no runtime memory)
    std::cout << "Bond 1 matches \"C=O\": "
              << CTSmarts::match_bond<"C=O">(molecule, get_bond(molecule, 1)) << std::endl;

    std::cout << "Bond 0 matches \"C=O\": "
              << CTSmarts::match_bond<"C=O">(molecule, get_bond(molecule, 0)) << std::endl;

    std::cout << "Bond 1 matches \"O=C\": "
              << CTSmarts::match_bond<"C=O">(molecule, get_bond(molecule, 1)) << std::endl;

    std::cout << "Is carboxyl double bond: "
              << isCarboxylDoubleBond(molecule, get_bond(molecule, 1)) << std::endl;
}
