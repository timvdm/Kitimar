#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>

#include <iostream>

using namespace Kitimar;

int main()
{
    auto molecule = Toolkit::readSmiles<Toolkit::openbabel>("CC(=O)[O-]");

    bool matches = CTSmarts::match<"C[O-]">(*molecule);

    std::cout << "Matching SMARTS \"C[O-]\" in SMILES \"CC(=O)[O-]\": "
              << matches << std::endl;
}
