#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecule.hpp>

#include <iostream>

// SMILES: Nc1ccc(N)cc1
// Generated by SmiToMock
auto molecule = Kitimar::Molecule::MockMolecule{
    {{7,0,0,1,2,false,false},{6,0,0,3,0,true,true},{6,0,0,2,1,true,true},{6,0,0,2,1,true,true},{6,0,0,3,0,true,true},{7,0,0,1,2,false,false},{6,0,0,2,1,true,true},{6,0,0,2,1,true,true}},
    {{0,1,1,false,false},{1,2,1,true,true},{2,3,2,true,true},{3,4,1,true,true},{4,5,1,false,false},{4,6,2,true,true},{6,7,1,true,true},{1,7,2,true,true}}
};

using namespace Kitimar;

int main()
{
    // Find a single mapping.

    auto [match, single] = CTSmarts::single<"c1ccccc1">(molecule);

    std::cout << "Mapping from SMARTS \"c1ccccc1\" in SMILES \"Nc1ccc(N)cc1\":" << std::endl;
    if (match)
        for (auto i = 0; i < single.size(); ++i)
            std::cout << i << " -> " << single[i] << std::endl;

    // Find all unique mappings.
    // A mapping is considered unique if it's atom set is unique.

    auto unique = CTSmarts::multi<"Nc1ccccc1">(molecule);

    std::cout << "Unique mappings from SMARTS \"Nc1ccccc1\" to SMILES \"Nc1ccc(N)cc1\":" << std::endl;
    auto i = 0;
    for (const auto &map : unique) {
        std::cout << "mapping " << i++ << ":" << std::endl;
        for (auto j = 0; j < map.size(); ++j)
            std::cout << "    " << j << " -> " << map[j] << std::endl;
    }

    // Find all mappings.

    auto all = CTSmarts::multi<"Nc1ccccc1">(molecule, CTSmarts::All);

    std::cout << "All mappings from SMARTS \"Nc1ccccc1\" to SMILES \"Nc1ccc(N)cc1\":" << std::endl;
    i = 0;
    for (const auto &map : unique) {
        std::cout << "mapping " << i++ << ":" << std::endl;
        for (auto j = 0; j < map.size(); ++j)
            std::cout << "    " << j << " -> " << map[j] << std::endl;
    }
}
