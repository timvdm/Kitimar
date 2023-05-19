
#include "TestData.hpp"

#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/RDKit/RDKit.hpp>



using namespace Kitimar;
using namespace Kitimar::CTLayout;



template<typename Layout>
void serializeKitimar(auto &source)
{
    CTLayout::MoleculeSink<Layout> sink{chembl_serialized_filename(Layout{})};
    for (auto mol : source.molecules()) {
        sink.write(mol);
        if (sink.size() % 1000 == 0) std::cout << sink.size() << std::endl;
    }
}

template<typename Layout>
void serializeOpenBabelSmilesToKitimar()
{
    OpenBabelSmilesMolSource source{chembl_smi_filename()};
    serializeKitimar<Layout>(source);
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
    serializeOpenBabelSmilesToKitimar<StructMoleculeIncident>();

    // RDKit SMILES -> RDKit Pickle
    //serializeRDKitSmilesToRDKit();


}
