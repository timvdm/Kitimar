#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>

#include "Test.hpp"
#include "TestData.hpp"

#include <catch2/catch_test_macros.hpp>


static_assert(Kitimar::Molecule::AtomList<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::BondList<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::MoleculeGraph<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::IncidentBondList<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::AdjacentAtomList<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::ElementLayer<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::IsotopeLayer<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::ChargeLayer<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::BondOrderLayer<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::ImplicitHydrogensLayer<OpenBabel::OBMol>);
static_assert(Kitimar::Molecule::AromaticLayer<OpenBabel::OBMol>);

using namespace Kitimar;


TEST_CASE("SmilesMolSource")
{
    OpenBabelSmilesMolSource source{chembl_smi_filename("1K")};

    auto numMolecules = 0;
    for (const auto &mol : source.molecules())
        ++numMolecules;

    CHECK(numMolecules == 1000);
}



TEST_CASE("HasImplicitH")
{
    using OpenBabelConfig = ctse::Config<ctse::DefaultImplicitH::ExactlyOne,
                                         ctse::Specialize::All,
                                         ctse::FullOptimizer,
                                         ctse::InverseMap>;



    auto mol = Toolkit::readSmiles<Toolkit::openbabel>("CC");

    std::cout << "CTSmarts:" << std::endl;
    for (const auto &map : ctse::maps_unique<"[h]", OpenBabelConfig>(mol))
        std::cout << "    " << map << std::endl;

    std::cout << "OpenBabel:" << std::endl;
    for (const auto &map : Toolkit::maps_unique<Toolkit::openbabel>("[h]", mol))
        std::cout << "    " << map << std::endl;


    auto smarts = OpenBabelConfig::transformSmarts(ctse::Smarts<"[h]">{});
    //identify_type<decltype(smarts.atoms)>{};
    //identify_type<decltype(OpenBabelConfig::transformSmarts(smarts).atoms)>{};

}
