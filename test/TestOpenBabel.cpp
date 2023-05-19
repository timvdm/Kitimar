#include <Kitimar/OpenBabel/OpenBabel.hpp>

#include "TestData.hpp"

#include <gtest/gtest.h>

static_assert(Kitimar::AtomList<OpenBabel::OBMol>);
static_assert(Kitimar::BondList<OpenBabel::OBMol>);
static_assert(Kitimar::MoleculeGraph<OpenBabel::OBMol>);
static_assert(Kitimar::IncidentBondList<OpenBabel::OBMol>);
static_assert(Kitimar::AdjacentAtomList<OpenBabel::OBMol>);
static_assert(Kitimar::ElementLayer<OpenBabel::OBMol>);
static_assert(Kitimar::IsotopeLayer<OpenBabel::OBMol>);
static_assert(Kitimar::ChargeLayer<OpenBabel::OBMol>);
static_assert(Kitimar::BondOrderLayer<OpenBabel::OBMol>);
static_assert(Kitimar::ImplicitHydrogensLayer<OpenBabel::OBMol>);
static_assert(Kitimar::AromaticLayer<OpenBabel::OBMol>);

using namespace Kitimar;


TEST(TestOpenBabel, SmilesMolSource)
{
    OpenBabelSmilesMolSource source{chembl_smi_filename("1K")};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules())
        ++numMolecules;

    EXPECT_EQ(numMolecules, 1000);
}
