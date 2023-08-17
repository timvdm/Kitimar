#include <Kitimar/CTLayout/Molecule.hpp>
#include <Kitimar/CTLayout/Sink.hpp>

#ifdef KITIMAR_WITH_OPENBABEL
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <openbabel/parsmart.h>
#endif

#include "TestData.hpp"

#include <catch2/catch_test_macros.hpp>


using namespace Kitimar;
using namespace Kitimar::CTLayout;

template<typename Layout, typename Sink>
void serialize(Molecule::Molecule auto &mol, Sink &sink)
{    
    serializeMolecule<Layout>(mol, sink);
}



template<typename MolObj>
void test_molecule()
{
    static_assert(Molecule::AtomList<MolObj>);
    static_assert(Molecule::BondList<MolObj>);
    static_assert(Molecule::MoleculeGraph<MolObj>);
    static_assert(Molecule::IncidentBondList<MolObj>);
    static_assert(Molecule::AdjacentAtomList<MolObj>);
    static_assert(Molecule::ElementLayer<MolObj>);
    static_assert(Molecule::IsotopeLayer<MolObj>);
    static_assert(Molecule::ChargeLayer<MolObj>);
    static_assert(Molecule::BondOrderLayer<MolObj>);
    //static_assert(Molecule::ImplicitHydrogensLayer<MolObj>);
    //static_assert(AromaticLayer<MolObj>);
}

TEST_CASE("Molecule")
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
    REQUIRE(num_atoms(mol) == num_atoms(ref));
    REQUIRE(num_bonds(mol) == num_bonds(ref));

    for (auto i = 0; i < num_atoms(ref); ++i) {
        auto refAtom = get_atom(ref, i);
        auto atom = get_atom(mol, i);
        REQUIRE(get_index(mol, atom) == get_index(ref, refAtom));
        REQUIRE(get_element(mol, atom) == get_element(ref, refAtom));
        REQUIRE(get_isotope(mol, atom) == get_isotope(ref, refAtom));
        REQUIRE(get_charge(mol, atom) == get_charge(ref, refAtom));
        REQUIRE(get_implicit_hydrogens(mol, atom) == get_implicit_hydrogens(ref, refAtom));
        REQUIRE(get_degree(mol, atom) == get_degree(ref, refAtom));
        REQUIRE(is_aromatic_atom(mol, atom) == is_aromatic_atom(ref, refAtom));
    }

    for (auto i = 0; i < num_bonds(ref); ++i) {
        auto refBond = get_bond(ref, i);
        auto bond = get_bond(mol, i);
        REQUIRE(get_index(mol, bond) == get_index(ref, refBond));
        REQUIRE(get_index(mol, get_source(mol, bond)) == get_index(ref, get_source(ref, refBond)));
        REQUIRE(get_index(mol, get_target(mol, bond)) == get_index(ref, get_target(ref, refBond)));
        REQUIRE(get_order(mol, bond) == get_order(ref, refBond));
        REQUIRE(is_aromatic_bond(mol, bond) == is_aromatic_bond(ref, refBond));
    }

    for (auto i = 0; i < num_atoms(ref); ++i) {
        auto refAtom = get_atom(ref, i);
        auto refBonds = get_bonds(ref, refAtom);
        auto refNbrs = get_nbrs(ref, refAtom);

        auto atom = get_atom(mol, i);
        auto bonds = get_bonds(mol, atom);
        auto nbrs = get_nbrs(mol, atom);

        REQUIRE(std::ranges::distance(bonds) == std::ranges::distance(refBonds));
        REQUIRE(std::ranges::distance(nbrs) == std::ranges::distance(refNbrs));

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

#ifdef KITIMAR_WITH_OPENBABEL

template<typename Layout>
void test_serialize(const std::string &smiles)
{
    auto obmol = readSmilesOpenBabel(smiles);

    std::vector<std::byte> data;
    StlVectorSink sink{data};
    serialize<Layout>(obmol, sink);

    auto source = PtrSource{data.data()};
    auto mol = toObject(Layout{}, source);

    compare(mol, obmol);
}

TEST_CASE("Serialize")
{
    test_serialize<StructMolecule>("CC(=O)[O-]");
    test_serialize<StructMoleculeIncident>("CC(=O)[O-]");
    test_serialize<StructMoleculeAdjacent>("CC(=O)[O-]");

    test_serialize<StructMoleculeIncidentAdjacent>("CC(=O)[O-]");

    test_serialize<StructMolecule>("c1ccccc1");
}

template<typename Layout>
auto serializeSmilesStlVectorSink(const std::string &SMILES, const std::string &filename)
{
    // Read SMILES
    auto obmol = readSmilesOpenBabel(SMILES);
    // Serialize SMILES
    std::vector<std::byte> data;
    StlVectorSink sink{data};
    serialize<Layout>(obmol, sink);
    // Write to file
    Util::writeFileData(filename, data);

    return obmol;
}

template<typename Layout>
auto serializeSmilesFileStreamSink(const std::string &SMILES, const std::string &filename)
{
    // Read SMILES
    auto obmol = readSmilesOpenBabel(SMILES);
    // Serialize SMILES
    FileStreamSink sink{filename};
    serialize<Layout>(obmol, sink);

    return obmol;
}


TEST_CASE("StlVectorSinkInMemorySource")
{
    using Layout = StructMoleculeIncident;
    // Compare
    auto filename = "Kitimar_StlVectorSinkInMemorySource.bin";
    auto ref = serializeSmilesStlVectorSink<Layout>("CC(=O)[O-]", filename);
    // Deserialize
    InMemorySource source{filename};
    auto mol = toObject(Layout{}, source);
    // Compare
    compare(mol, ref);
}

TEST_CASE("StlVectorSinkMemoryMappedSource")
{
    using Layout = StructMoleculeIncident;
    // Compare
    auto filename = "Kitimar_StlVectorSinkMemoryMappedSource.bin";
    auto ref = serializeSmilesStlVectorSink<Layout>("CC(=O)[O-]", filename);
    // Deserialize
    MemoryMappedSource source{filename};
    auto mol = toObject(Layout{}, source);
    // Compare
    compare(mol, ref);
}

TEST_CASE("FileStreamSinkMemoryMappedSource")
{
    using Layout = StructMoleculeIncident;
    // Compare
    auto filename = "Kitimar_FileStreamSinkMemoryMappedSource.bin";
    auto ref = serializeSmilesFileStreamSink<Layout>("CC(=O)[O-]", filename);
    // Deserialize
    InMemorySource source{filename};
    auto mol = toObject(Layout{}, source);
    // Compare
    compare(mol, ref);
}

template<typename Layout>
void test_validate()
{
    auto suffix = "1K";
    //auto suffix = "";
    auto path = chembl_serialized_filename(Layout{}, suffix);
    FileStreamSink sink{path};

    OpenBabelSmilesMolSource obMolSource{chembl_smi_filename(suffix)};
    std::cout << "counting molecules..." << std::endl;
    auto numMolecules = obMolSource.numMolecules();
    std::cout << "# molecules: " << numMolecules << std::endl;
    std::cout << "serializing molecules..." << std::endl;
    serializeMolSource<Layout>(obMolSource, sink);
    sink.close();
    obMolSource.reset();

    std::cout << "validating molecules..." << std::endl;
    auto source = InMemorySource{path};
    BytePtrMolSource<Layout> molSource{source.toPtrSource()};

    REQUIRE(molSource.numMolecules() == obMolSource.numMolecules());

    for (auto i = 0; i < numMolecules; ++i) {
        auto obmol = obMolSource.read();
        auto mol = molSource.read();
        compare(mol, obmol);
    }
}

TEST_CASE("Validate")
{
    test_validate<Vector<StructMoleculeIncident>>();
    test_validate<Vector<ListMoleculeIncident>>();
    //test_validate<TypeMolecules>();
}


#endif


/*
TEST_CASE("ValidateMemory)
{
    using Layout = StructMoleculeIncident;

    MemoryMappedSource<Layout> mmap{chembl_serialized_filename(Layout{})};
    InMemorySource<Layout> mem{chembl_serialized_filename(Layout{})};

    REQUIRE(mem.size(), mmap.size());

    for (auto i = 0; i < mmap.size(); ++i) {
        auto memMol = mem.read();
        auto mmapMol = mmap.read();
        compare(memMol, mmapMol);
    }
}
*/

//TEST_CASE("CompareMemoryMappedInMemory)
//{
//    using Layout = StructMoleculeIncident;

//    MemoryMappedSource<Layout> mmap{chembl_serialized_filename(Layout{})};
//    InMemorySource<Layout> mem{chembl_serialized_filename(Layout{})};

//    REQUIRE(mem.size(), mmap.size());

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
//    Isomorphism<SMARTS, MapType::Single> smarts{};

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


//    REQUIRE(ctMatches.size(), obMatches.size());

//    REQUIRE(std::ranges::count(ctMatches, true),
//              std::ranges::count(obMatches, true));

//}


//#define VALIDATE_OPENBABEL 0


//TEST_CASE("ValidateOpenBabel)
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

//    auto iso = Isomorphism<"c1ccccc1O", MapType::Single>{};

//    //std::cout << "phenol: " << iso.match(mol, get_atom(mol, 0)) << std::endl;
//    std::cout << "phenol: " << iso.match(mol) << std::endl;

//}




