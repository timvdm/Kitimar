#include <Kitimar/RDKit/RDKit.hpp>

#include "TestData.hpp"

#include <gtest/gtest.h>

static_assert(Kitimar::Molecule::IsAtomList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::IsBondList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::MoleculeGraph<RDKit::ROMol>);
static_assert(Kitimar::Molecule::IncidentBondList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::AdjacentAtomList<RDKit::ROMol>);
static_assert(Kitimar::Molecule::ElementLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::IsotopeLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::ChargeLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::BondOrderLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::ImplicitHydrogensLayer<RDKit::ROMol>);
static_assert(Kitimar::Molecule::AromaticLayer<RDKit::ROMol>);

using namespace Kitimar;


TEST(TestRDKit, SmilesMolSource)
{
    RDKitSmilesMolSource source{chembl_smi_filename("1K")};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules()) {
        ASSERT_TRUE(mol.get());
        ++numMolecules;
    }

    EXPECT_EQ(numMolecules, 1000);
}

TEST(TestRDKit, PickleMolSource)
{
    RDKitPickleMolSource source{chembl_rdkit_filename()};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules()) {
        ASSERT_TRUE(mol.get());
        ++numMolecules;
    }

    EXPECT_EQ(numMolecules, 2327928);
}
