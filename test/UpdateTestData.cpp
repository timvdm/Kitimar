
#include "TestData.hpp"

#ifdef KITIMAR_WITH_OPENBABEL
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#endif

#ifdef KITIMAR_WITH_RDKIT
#include <Kitimar/RDKit/RDKit.hpp>
#endif


using namespace Kitimar;
using namespace Kitimar::CTLayout;

#ifdef KITIMAR_WITH_OPENBABEL

template<typename Layout>
void serializeOpenBabelSmilesToKitimar(const std::string &suffix = {})
{
    OpenBabelSmilesMolSource source{chembl_smi_filename(suffix)};
    FileStreamSink sink{chembl_serialized_filename(Layout{}, suffix)};
    serializeMolSource<Layout>(source, sink);
}


template<typename Layout>
void serializeOpenBabelSmilesToKitimarTypeIndex()
{
    OpenBabelSmilesMolSource source{chembl_smi_filename()};

    std::vector<AtomType> atomTypes;
    std::vector<BondType> bondTypes;

    // Pass 1: Determine types
    for (auto mol : source.molecules()) {
        for (auto atom : get_atoms(*mol)) {
            auto type = AtomType{
                static_cast<uint8_t>(get_element(*mol, atom)),
                static_cast<uint8_t>(get_isotope(*mol, atom)),
                static_cast<int8_t>(get_charge(*mol, atom)),
                static_cast<uint8_t>(get_degree(*mol, atom)),
                static_cast<uint8_t>(get_implicit_hydrogens(*mol, atom)),
                is_aromatic_atom(*mol, atom)
            };

            if (std::ranges::find(atomTypes, type) == atomTypes.end())
                atomTypes.push_back(type);
        }

        for (auto bond : get_bonds(*mol)) {
            auto type = BondType{
                static_cast<uint8_t>(get_order(*mol, bond)),
                is_ring_bond(*mol, bond),
                is_aromatic_bond(*mol, bond)
            };

            if (std::ranges::find(bondTypes, type) == bondTypes.end())
                bondTypes.push_back(type);
        }

    }

    std::cout << "# atom types: " << atomTypes.size() << std::endl;
    std::cout << "# bond types: " << bondTypes.size() << std::endl;

}

#endif // KITIMAR_WITH_OPENBABEL

#ifdef KITIMAR_WITH_RDKIT

void serializeRDKit(auto &source, const std::string &suffix)
{
    std::ofstream ofs(chembl_rdkit_filename(suffix), std::ios_base::binary);
    std::string pickle;
    int n = 0;
    for (auto mol : source.molecules()) {
        RDKit::MolPickler::pickleMol(*mol, ofs);
        n++;
        if (n % 1000 == 0) std::cout << n << std::endl;
    }
}

void serializeRDKitSmilesToRDKit(const std::string &suffix = {})
{
    RDKitSmilesMolSource source{chembl_smi_filename(suffix)};
    serializeRDKit(source, suffix);
}

#endif // #ifdef KITIMAR_WITH_RDKIT


int main()
{
#ifdef KITIMAR_WITH_OPENBABEL
    // OpenBabel SMILES -> Kitimar
    //serializeOpenBabelSmilesToKitimar<Vector<StructMolecule>>();
    //serializeOpenBabelSmilesToKitimar<Vector<StructMoleculeIncident>>();
    //serializeOpenBabelSmilesToKitimar<Vector<ListMoleculeIncident>>();
    //serializeOpenBabelSmilesToKitimar<TypeMolecules>();

    //serializeOpenBabelSmilesToKitimarTypeIndex<StructMoleculeIncident>();
#endif

#ifdef KITIMAR_WITH_RDKIT
    // RDKit SMILES -> RDKit Pickle
    serializeRDKitSmilesToRDKit("1K");
#endif


}
