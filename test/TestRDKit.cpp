#include <Kitimar/RDKit/RDKit.hpp>

#include "TestData.hpp"

#include <gtest/gtest.h>

static_assert(Kitimar::IsAtomList<RDKit::ROMol*>);
static_assert(Kitimar::IsBondList<RDKit::ROMol*>);
static_assert(Kitimar::MoleculeGraph<RDKit::ROMol*>);
static_assert(Kitimar::IncidentBondList<RDKit::ROMol*>);
static_assert(Kitimar::AdjacentAtomList<RDKit::ROMol*>);
static_assert(Kitimar::ElementLayer<RDKit::ROMol*>);
static_assert(Kitimar::IsotopeLayer<RDKit::ROMol*>);
static_assert(Kitimar::ChargeLayer<RDKit::ROMol*>);
static_assert(Kitimar::BondOrderLayer<RDKit::ROMol*>);
static_assert(Kitimar::ImplicitHydrogensLayer<RDKit::ROMol*>);
static_assert(Kitimar::AromaticLayer<RDKit::ROMol*>);

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
