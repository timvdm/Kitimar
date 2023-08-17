#include <Kitimar/RDKit/RDKit.hpp>

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
