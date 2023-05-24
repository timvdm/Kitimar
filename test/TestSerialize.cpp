#include <Kitimar/CTLayout/Molecule.hpp>
#include <Kitimar/CTLayout/Sink.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/CTSmarts/Isomorphism.hpp>

#include "TestData.hpp"
//#include "Test.hpp"

#include <openbabel/parsmart.h>

#include <gtest/gtest.h>

using namespace Kitimar;
using namespace Kitimar::CTLayout;
using namespace Kitimar::CTSmarts;

template<typename Layout, typename Sink>
void serialize(Molecule auto &mol, Sink &sink)
{    
    serializeMolecule<Layout>(mol, sink);
}



template<typename MolObj>
void test_molecule()
{
    static_assert(IsAtomList<MolObj>);
    static_assert(IsBondList<MolObj>);
    static_assert(MoleculeGraph<MolObj>);
    static_assert(IncidentBondList<MolObj>);
    static_assert(AdjacentAtomList<MolObj>);
    static_assert(ElementLayer<MolObj>);
    static_assert(IsotopeLayer<MolObj>);
    static_assert(ChargeLayer<MolObj>);
    static_assert(BondOrderLayer<MolObj>);
    static_assert(ImplicitHydrogensLayer<MolObj>);
    //static_assert(AromaticLayer<MolObj>);
}

TEST(TestSerialize, Molecule)
{
    test_molecule<StructObject<StructMolecule>>();

    test_molecule<StructObject<StructMoleculeIncident>>();
    test_molecule<StructObject<StructMoleculeAdjacent>>();
    test_molecule<StructObject<StructMoleculeIncidentAdjacent>>();

    /*
    test_molecule<StructObject<ArrayMolecule>>();
    test_molecule<StructObject<ArrayMoleculeIncident>>();
    test_molecule<StructObject<ArrayMoleculeAdjacent>>();
    test_molecule<StructObject<ArrayMoleculeIncidentAdjacent>>();
    */
}


void compare(auto &mol, auto &ref)
{
    ASSERT_EQ(num_atoms(mol), num_atoms(ref));
    ASSERT_EQ(num_bonds(mol), num_bonds(ref));

    for (auto i = 0; i < num_atoms(ref); ++i) {
        auto refAtom = get_atom(ref, i);
        auto atom = get_atom(mol, i);
        ASSERT_EQ(get_index(mol, atom), get_index(ref, refAtom));
        ASSERT_EQ(get_element(mol, atom), get_element(ref, refAtom));
        ASSERT_EQ(get_isotope(mol, atom), get_isotope(ref, refAtom));
        ASSERT_EQ(get_charge(mol, atom), get_charge(ref, refAtom));
        ASSERT_EQ(get_implicit_hydrogens(mol, atom), get_implicit_hydrogens(ref, refAtom));
        ASSERT_EQ(get_degree(mol, atom), get_degree(ref, refAtom));
        ASSERT_EQ(is_aromatic_atom(mol, atom), is_aromatic_atom(ref, refAtom));
    }

    for (auto i = 0; i < num_bonds(ref); ++i) {
        auto refBond = get_bond(ref, i);
        auto bond = get_bond(mol, i);
        ASSERT_EQ(get_index(mol, bond), get_index(ref, refBond));
        ASSERT_EQ(get_index(mol, get_source(mol, bond)),
                  get_index(ref, get_source(ref, refBond)));
        ASSERT_EQ(get_index(mol, get_target(mol, bond)),
                  get_index(ref, get_target(ref, refBond)));
        ASSERT_EQ(get_order(mol, bond), get_order(ref, refBond));
        ASSERT_EQ(is_aromatic_bond(mol, bond), is_aromatic_bond(ref, refBond));
    }

    for (auto i = 0; i < num_atoms(ref); ++i) {
        auto refAtom = get_atom(ref, i);
        auto refBonds = get_bonds(ref, refAtom);
        auto refNbrs = get_nbrs(ref, refAtom);

        auto atom = get_atom(mol, i);
        auto bonds = get_bonds(mol, atom);
        auto nbrs = get_nbrs(mol, atom);

        ASSERT_EQ(std::ranges::distance(bonds), std::ranges::distance(refBonds));
        ASSERT_EQ(std::ranges::distance(nbrs), std::ranges::distance(refNbrs));

        // not sorted...
        /*
        for (auto b : bonds)
            std::cout << get_index(mol, b) << " ";
        std::cout << std::endl;

        for (auto b : obbonds)
            std::cout << get_index(obmol, b) << " ";
        std::cout << std::endl;


        assert(std::ranges::equal(bonds, obbonds, {},
            [&mol] (auto bond) { return get_index(mol, bond); },
            [&obmol] (auto obbond) { return get_index(obmol, obbond); }));

        assert(std::ranges::equal(nbrs, obnbrs, {},
            [&mol] (auto nbr) { return get_index(mol, nbr); },
            [&obmol] (auto obnbr) { return get_index(obmol, obnbr); }));
        */
    }
}

template<typename Layout>
void test_serialize(const std::string &smiles)
{
    auto obmol = readSmilesOpenBabel(smiles);

    std::vector<std::byte> data;
    StlVectorSink sink{data};
    serialize<Layout>(obmol, sink);

    auto source = BytePtrSource{data.data()};
    auto mol = toObject(Layout{}, source);

    compare(mol, obmol);
}

TEST(TestSerialize, Serialize)
{
    test_serialize<StructMolecule>("CC(=O)[O-]");
    test_serialize<StructMoleculeIncident>("CC(=O)[O-]");
    test_serialize<StructMoleculeAdjacent>("CC(=O)[O-]");

    test_serialize<StructMoleculeIncidentAdjacent>("CC(=O)[O-]");

    test_serialize<StructMolecule>("c1ccccc1");
}

template<typename Layout>
auto serializeSmiles(const std::string &SMILES, const std::string &filename)
{
    // Read SMILES
    auto obmol = readSmilesOpenBabel(SMILES);

    // Serialize SMILES
    std::vector<std::byte> data;
    StlVectorSink sink{data};
    serialize<Layout>(obmol, sink);

    auto source = BytePtrSource{data.data()};
    auto ref = toObject(Layout{}, source);

    // Write to file
    std::ofstream ofs(filename, std::ios_base::binary | std::ios_base::out);
    ofs.write(reinterpret_cast<char*>(data.data()), data.size());
    ofs.close();

    return obmol;
}


TEST(TestSerialize, InMemorySource)
{
    using Layout = StructMoleculeIncident;
    // Compare
    auto filename = "test_mmap.bin";
    auto ref = serializeSmiles<Layout>("CC(=O)[O-]", filename);
    // Deserialize
    InMemorySource source{filename};
    auto mol = toObject(Layout{}, source);
    // Compare
    compare(mol, ref);
}

TEST(TestSerialize, MemoryMappedSource)
{
    using Layout = StructMoleculeIncident;
    // Compare
    auto filename = "test_mmap.bin";
    auto ref = serializeSmiles<Layout>("CC(=O)[O-]", filename);
    // Deserialize
    MemoryMappedSource source{filename};
    auto mol = toObject(Layout{}, source);
    // Compare
    compare(mol, ref);
}

/*
TEST(TestSerialize, ValidateMemory)
{
    using Layout = StructMoleculeIncident;

    MemoryMappedSource<Layout> mmap{chembl_serialized_filename(Layout{})};
    InMemorySource<Layout> mem{chembl_serialized_filename(Layout{})};

    ASSERT_EQ(mem.size(), mmap.size());

    for (auto i = 0; i < mmap.size(); ++i) {
        auto memMol = mem.read();
        auto mmapMol = mmap.read();
        compare(memMol, mmapMol);
    }
}
*/

//TEST(TestSerialize, CompareMemoryMappedInMemory)
//{
//    using Layout = StructMoleculeIncident;

//    MemoryMappedSource<Layout> mmap{chembl_serialized_filename(Layout{})};
//    InMemorySource<Layout> mem{chembl_serialized_filename(Layout{})};

//    ASSERT_EQ(mem.size(), mmap.size());

//    for (auto i = 0; i < mmap.size(); ++i) {
//        auto memMol = mem.read();
//        auto mmapMol = mmap.read();
//        compare(memMol, mmapMol);
//    }
//}




//auto findOBMatches(const std::string &SMARTS, const std::string &filename = chembl_smi_filename())
//{
//    OpenBabel::OBSmartsPattern smarts;
//    smarts.Init(SMARTS);

//    std::vector<bool> matches;
//    OpenBabelSmilesMolSource source{filename};
//    for (auto mol : source.molecules())
//        matches.push_back(smarts.Match(mol));

//    return matches;
//}

//template<ctll::fixed_string SMARTS, typename Layout>
//auto findCTMatches()
//{
//    FileStreamSource<Layout> source{chembl_serialized_filename(Layout{})};
//    SingleIsomorphism<SMARTS> smarts{};

//    std::vector<bool> matches;
//    for (auto mol : source.objects()) {
//        matches.push_back(smarts.match(mol));
//    }

//    return matches;
//}

//template<ctll::fixed_string SMARTS, typename Layout>
//auto validateOpenBabel()
//{
//    auto ctMatches = findCTMatches<SMARTS, Layout>();
//    auto obMatches = findOBMatches(std::string{SMARTS.begin(), SMARTS.end()});


//    ASSERT_EQ(ctMatches.size(), obMatches.size());

//    ASSERT_EQ(std::ranges::count(ctMatches, true),
//              std::ranges::count(obMatches, true));

//}


//#define VALIDATE_OPENBABEL 0


//TEST(TestSerialize, ValidateOpenBabel)
//{
//    if (!VALIDATE_OPENBABEL)
//        return;

//    validateOpenBabel<"C", StructMoleculeIncident>();

//    validateOpenBabel<"*1~*~*~*~*~*~1", StructMoleculeIncident>();
//    validateOpenBabel<"c1ccccc1-c2ccccc2", StructMoleculeIncident>();
//    validateOpenBabel<"c1ccccc1CCCCc1ccccc1", StructMoleculeIncident>();
//    validateOpenBabel<"c1cc2c(cc1)cccc2", StructMoleculeIncident>();
//    validateOpenBabel<"Clc1ccccc1", StructMoleculeIncident>();
//    validateOpenBabel<"BrCCCCCCCCC", StructMoleculeIncident>();
//    //validateOpenBabel<"[nH]1ccc2c1cccc2", StructMoleculeIncident>();
//    //validateOpenBabel<"[nH]", StructMoleculeIncident>();
//    //validateOpenBabel<"c1cc(=O)cc[nH]1", StructMoleculeIncident>();
//    validateOpenBabel<"O1CCOC12CCNCC2", StructMoleculeIncident>();

//}





















//void test_phenol()
//{
//    auto obmol = readSmilesOpenBabel("c1ccccc1O");
//    std::vector<std::byte> data;
//    serialize<StructMolecule>(obmol, data);

//    auto mol = Object<StructMolecule>(data.data());

//    auto iso = SingleIsomorphism<"c1ccccc1O">{};

//    //std::cout << "phenol: " << iso.match(mol, get_atom(mol, 0)) << std::endl;
//    std::cout << "phenol: " << iso.match(mol) << std::endl;

//}




