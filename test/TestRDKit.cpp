#include <Kitimar/RDKit/RDKit.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>

#include "TestData.hpp"

#include <catch2/catch_test_macros.hpp>

static_assert(Kitimar::Molecule::AtomList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::BondList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::MoleculeGraph<RDKit::ROMol>);
static_assert(Kitimar::Molecule::IncidentBondList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::AdjacentAtomList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::ElementLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::IsotopeLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::ChargeLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::BondOrderLayer<RDKit::ROMol>);
//static_assert(Kitimar::Molecule::ImplicitHydrogensLayer<RDKit::ROMol>); // FIXME
static_assert(Kitimar::Molecule::AromaticLayer<RDKit::ROMol>);

using namespace Kitimar;


TEST_CASE("SmilesMolSource")
{
    RDKitSmilesMolSource source{chembl_smi_filename("1K")};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules()) {
        REQUIRE(mol.get());
        ++numMolecules;
    }

    CHECK(numMolecules == 1000);
}

TEST_CASE("PickleMolSource")
{
    RDKitPickleMolSource source{chembl_rdkit_filename("1K")};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules()) {
        REQUIRE(mol.get());
        ++numMolecules;
    }

    CHECK(numMolecules == 1000);
}


#define TEST_MATCH(SMILES, SMARTS) \
    { \
        auto mol = Toolkit::readSmiles<Toolkit::rdkit>(SMILES); \
        CHECK(ctse::match<SMARTS>(*mol) == Toolkit::match<Toolkit::rdkit>(SMARTS, *mol)); \
    }

TEST_CASE("Match")
{
    TEST_MATCH("CCOC(=O)c1cc2cc(C(=O)O)ccc2[nH]1", "[nX2r5]");

    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[cD3H0;r6;R2]");
    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[*;r5,r6;R1]");
    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[cD3H0;r6;R2]");
    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[cD3H0;r6;R2][*;r5,r6;R1]");
    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[*;r5,r6;R1][cD3H0;r6;R2]");
    TEST_MATCH("COc1c2ccoc2cc2oc(C)cc(=O)c12", "[cD3H0;r6;R2][*;r5,r6;R1][cD3H0;r6;R2]");
}

TEST_CASE("Debug")
{
    /*
    auto mol = Toolkit::readSmiles<Toolkit::rdkit>("COc1c2ccoc2cc2oc(C)cc(=O)c12");

    // rdkit
    {
        std::cout << "RDKit:" << std::endl;
        auto maps = Toolkit::maps_unique<Toolkit::rdkit>("[r6]", *mol); //
        //auto maps = Toolkit::maps_unique<Toolkit::rdkit>("[cD3H0]", *mol); // 2 3 7 9 11 14 16
        //auto maps = Toolkit::maps_unique<Toolkit::rdkit>("[cD3H0;r6]", *mol); // 2 9 11 14 16
        //auto maps = Toolkit::maps_unique<Toolkit::rdkit>("[cD3H0;r6;R2]", *mol); // 9 16
        //auto maps = Toolkit::maps_unique<Toolkit::rdkit>("[*;r5,r6;R1]", *mol); // 2 4 5 6 8 10 11 13 14
        for (const auto &map : maps) {
            std::cout << "    [ ";
            for (const auto &p : map)
                std::cout << p.second << " ";
            std::cout << "]" << std::endl;
        }
    }

    // ctse
    {
        std::cout << "CTSmarts:" << std::endl;
        auto maps = ctse::maps_unique<"[r6]">(*mol); //
        //auto maps = ctse::maps_unique<"[cD3H0]">(*mol); // 2 3 7 9 11 14 16
        //auto maps = ctse::maps_unique<"[cD3H0;r6]">(*mol); // 2 3 7 9 11 14 16
        //auto maps = ctse::maps_unique<"[cD3H0;r6;R2]">(*mol); // 3 7 9 16
        //auto maps = ctse::maps_unique<"[*;r5,r6;R1]">(*mol); // 2 4 5 6 8 10 11 13 14
        //auto maps = ctse::maps_unique<"[cD3H0;r6;R2][*;r5,r6;R1][cD3H0;r6;R2]">(*mol);
        for (const auto &map : maps)
            std::cout << "    " << map << std::endl;
    }

    auto smarts = ctse::Smarts<"[cD3H0;r6;R2][*;r5,r6;R1][cD3H0;r6;R2]">{};
    */
}
