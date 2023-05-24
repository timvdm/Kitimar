
#include "TestData.hpp"

#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/RDKit/RDKit.hpp>



using namespace Kitimar;
using namespace Kitimar::CTLayout;



template<typename Layout>
void serializeOpenBabelSmilesToKitimar()
{
    OpenBabelSmilesMolSource source{chembl_smi_filename("1K")};
    serializeMolSource<Layout>(source, chembl_serialized_filename(Layout{}));
}


template<typename Layout>
void serializeOpenBabelSmilesToKitimarTypeIndex()
{
    OpenBabelSmilesMolSource source{chembl_smi_filename()};

    std::vector<AtomType> atomTypes;
    std::vector<BondType> bondTypes;

    // Pass 1: Determine types
    for (auto mol : source.molecules()) {
        for (auto atom : get_atoms(mol)) {
            auto type = AtomType{
                static_cast<uint8_t>(get_element(mol, atom)),
                static_cast<uint8_t>(get_isotope(mol, atom)),
                static_cast<int8_t>(get_charge(mol, atom)),
                static_cast<uint8_t>(get_degree(mol, atom)),
                static_cast<uint8_t>(get_implicit_hydrogens(mol, atom)),
                is_aromatic_atom(mol, atom)
            };

            if (std::ranges::find(atomTypes, type) == atomTypes.end())
                atomTypes.push_back(type);
        }

        for (auto bond : get_bonds(mol)) {
            auto type = BondType{
                static_cast<uint8_t>(get_order(mol, bond)),
                is_cyclic_bond(mol, bond),
                is_aromatic_bond(mol, bond)
            };

            if (std::ranges::find(bondTypes, type) == bondTypes.end())
                bondTypes.push_back(type);
        }

    }

    std::cout << "# atom types: " << atomTypes.size() << std::endl;
    std::cout << "# bond types: " << bondTypes.size() << std::endl;

}




void serializeRDKit(auto &source)
{
    std::ofstream ofs(chembl_rdkit_filename(), std::ios_base::binary);
    std::string pickle;
    int n = 0;
    for (auto mol : source.molecules()) {
        RDKit::MolPickler::pickleMol(*mol, ofs);
        n++;
        if (n % 1000 == 0) std::cout << n << std::endl;
    }
}

void serializeRDKitSmilesToRDKit()
{
    RDKitSmilesMolSource source{chembl_smi_filename()};
    serializeRDKit(source);
}




int main()
{
    // OpenBabel SMILES -> Kitimar
    serializeOpenBabelSmilesToKitimar<Vector<StructMoleculeIncident>>();
    //serializeOpenBabelSmilesToKitimar<AtomTypeMolecule>();

    //serializeOpenBabelSmilesToKitimarTypeIndex<StructMoleculeIncident>();

    // RDKit SMILES -> RDKit Pickle
    //serializeRDKitSmilesToRDKit();


}
