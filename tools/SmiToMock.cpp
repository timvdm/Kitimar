#include <Kitimar/Molecule/MockMolecule.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/Util/Util.hpp>

using namespace Kitimar;

int main(int argc, char **argv)
{
    /*
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <SMILES>" << std::endl;
        return -1;
    }

    auto SMILES = argv[1];
    */
    auto SMILES = "OCC(O)CO";
    auto mol = readSmilesOpenBabel(SMILES);

    fmt::println("// SMILES: {}", SMILES);
    fmt::println("// Generated by SmiToMock", SMILES);
    fmt::println("auto molecule = Kitimar::Molecule::MockMolecule{{");
    // Atoms
    fmt::print("    {{");
    for (auto atom : get_atoms(mol)) {
        if (get_index(mol, atom))
            fmt::print(",");
        fmt::print("{{{},{},{},{},{},{},{}}}",
                   get_element(mol, atom),
                   get_isotope(mol, atom),
                   get_charge(mol, atom),
                   get_degree(mol, atom),
                   get_implicit_hydrogens(mol, atom),
                   is_cyclic_atom(mol, atom),
                   is_aromatic_atom(mol, atom));
    }
    fmt::println("}},");
    // Bond
    fmt::print("    {{");
    for (auto bond : get_bonds(mol)) {
        if (get_index(mol, bond))
            fmt::print(",");
        fmt::print("{{{},{},{},{},{}}}",
                   get_index(mol, get_source(mol, bond)),
                   get_index(mol, get_target(mol, bond)),
                   get_order(mol, bond),
                   is_cyclic_bond(mol, bond),
                   is_aromatic_bond(mol, bond));
    }
    fmt::println("}}");
    fmt::println("}};");
}
