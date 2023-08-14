#include "TestCTSmarts.hpp"
#include <Kitimar/Util/Util.hpp>

#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <array>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;
using Molecule::MockButane;

//
// maps_unique(mol)
//

using MapsUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Maps<Map<0>, Map<1>>>,
    TypeTestCase<MockAcetateAnion, "O", Maps<Map<2>, Map<3>>>,
    TypeTestCase<MockAcetateAnion, "[O-]", Maps<Map<3>>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "[O+]", EmptyMaps>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Maps<Map<0, 1>>>,
    TypeTestCase<MockAcetateAnion, "C=O", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "C[O-]", Maps<Map<1, 3>>>,
    TypeTestCase<MockAcetateAnion, "C~O", Maps<Map<1, 2>, Map<1, 3>>>,
    TypeTestCase<MockAcetateAnion, "*~*", Maps<Map<0, 1>, Map<1, 2>, Map<1, 3>>>,

    TypeTestCase<MockAcetateAnion, "C=C", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "C#O", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "C[O+]", EmptyMaps>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC(=O)[O-]", Maps<Map<0, 1, 2, 3>>>,
    TypeTestCase<MockAcetateAnion, "*~*~*", Maps<Map<0, 1, 2>, Map<0, 1, 3>, Map<2, 1, 3>>>,
    TypeTestCase<MockAcetateAnion, "*~*(~*)~*", Maps<Map<0, 1, 2, 3>>>,

    TypeTestCase<MockAcetateAnion, "CC(=O)N", EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_unique()
{
    auto mol = Case::template mol<Mol>();    
    CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, ctse::Unique)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_unique<Case::smarts>(mol)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
}

TEMPLATE_LIST_TEST_CASE("maps_unique", "", MapsUniqueCases)
{
    test_maps_unique<Molecule::MockIndexMolecule, TestType>();
    test_maps_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// maps_all(mol)
//

using MapsAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Maps<Map<0>, Map<1>>>,
    TypeTestCase<MockAcetateAnion, "O", Maps<Map<2>, Map<3>>>,
    TypeTestCase<MockAcetateAnion, "[O-]", Maps<Map<3>>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "[O+]", EmptyMaps>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Maps<Map<0, 1>, Map<1, 0>>>,
    TypeTestCase<MockAcetateAnion, "C=O", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "C[O-]", Maps<Map<1, 3>>>,
    TypeTestCase<MockAcetateAnion, "C~O", Maps<Map<1, 2>, Map<1, 3>>>,
    TypeTestCase<MockAcetateAnion, "*~*", Maps<Map<0, 1>, Map<1, 0>, Map<1, 2>, Map<2, 1>, Map<1, 3>, Map<3, 1>>>,

    TypeTestCase<MockAcetateAnion, "C=C", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "C#O", EmptyMaps>,
    TypeTestCase<MockAcetateAnion, "C[O+]", EmptyMaps>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC(=O)[O-]", Maps<Map<0, 1, 2, 3>>>,
    TypeTestCase<MockAcetateAnion, "O=CO", Maps<Map<2, 1, 3>>>,
    TypeTestCase<MockAcetateAnion, "OC=O", Maps<Map<3, 1, 2>>>,
    TypeTestCase<MockAcetateAnion, "*~*~*", Maps<Map<0, 1, 2>, Map<2, 1, 0>, Map<0, 1, 3>, Map<3, 1, 0>, Map<2, 1, 3>, Map<3, 1, 2>>>,
    TypeTestCase<MockAcetateAnion, "*~*(~*)~*", Maps<Map<0, 1, 2, 3>, Map<0, 1, 3, 2>, Map<2, 1, 0, 3>, Map<2, 1, 3, 0>, Map<3, 1, 0, 2>, Map<3, 1, 2, 0>>>,

    TypeTestCase<MockAcetateAnion, "CC(=O)N", EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_all()
{
    auto mol = Case::template mol<Mol>();
    CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, ctse::All)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_all<Case::smarts>(mol)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
}

TEMPLATE_LIST_TEST_CASE("maps_all", "", MapsAllCases)
{
    test_maps_all<Molecule::MockIndexMolecule, TestType>();
    test_maps_all<Molecule::MockProxyMolecule, TestType>();
}

//
// maps_atom_unique(mol, atom)
//

using MapsAtomUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Maps<Map<1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*", 1, Maps<Map<1, 0>, Map<1, 2>, Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 1, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=C", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C#O", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C[O+]", 1, EmptyMaps>,
    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC([O-])=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "O=C(C)[O-]", 2, Maps<Map<2, 1, 0, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 1, Maps<Map<1, 0, 2>, Map<1, 0, 3>, Map<1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)N", 0, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_atom_unique()
{
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    CHECK_THAT(toVectorMaps(ctse::maps_atom<Case::smarts>(mol, atom)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_atom<Case::smarts>(mol, atom, ctse::Unique)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_atom_unique<Case::smarts>(mol, atom)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, atom)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, atom, ctse::Unique)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps_unique<Case::smarts>(mol, atom)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    }
}

TEMPLATE_LIST_TEST_CASE("maps_atom_unique", "", MapsAtomUniqueCases)
{
    test_maps_atom_unique<Molecule::MockIndexMolecule, TestType>();
    test_maps_atom_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// maps_atom_all(mol, atom)
//

using MapsAtomAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Maps<Map<1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*", 1, Maps<Map<1, 0>, Map<1, 2>, Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 1, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=C", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C#O", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C[O+]", 1, EmptyMaps>,
    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC([O-])=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "O=C(C)[O-]", 2, Maps<Map<2, 1, 0, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 1, Maps<Map<1, 0, 2>, Map<1, 0, 3>, Map<1, 2, 0>, Map<1, 2, 3>, Map<1, 3, 0>, Map<1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 0, Maps<Map<0, 1, 2, 3>, Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)N", 0, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_atom_all()
{
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    CHECK_THAT(toVectorMaps(ctse::maps_atom<Case::smarts>(mol, atom, ctse::All)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_atom_all<Case::smarts>(mol, atom)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, atom, ctse::All)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps_all<Case::smarts>(mol, atom)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    }
}

TEMPLATE_LIST_TEST_CASE("maps_atom_all", "", MapsAtomAllCases)
{
    test_maps_atom_all<Molecule::MockIndexMolecule, TestType>();
    test_maps_atom_all<Molecule::MockProxyMolecule, TestType>();
}

//
// maps_bond_unique(mol, bond)
//

using MapsBondUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC~O", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 1, Maps<Map<2, 1, 0>, Map<2, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 2, Maps<Map<3, 1, 0>, Map<3, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 0, Maps<Map<0, 1, 2, 3>>>,
    //
    // CCCC
    //
    IndexTypeTestCase<MockButane, "****", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockButane, "****", 1, EmptyMaps>,
    IndexTypeTestCase<MockButane, "****", 2, Maps<Map<3, 2, 1, 0>>>,
    IndexTypeTestCase<MockButane, "*(**)*", 0, EmptyMaps>,
    IndexTypeTestCase<MockButane, "*(**)*", 1, Maps<Map<1, 2, 3, 0>>>,
    IndexTypeTestCase<MockButane, "*(**)*", 2, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_bond_unique()
{
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    CHECK_THAT(toVectorMaps(ctse::maps_bond<Case::smarts>(mol, bond)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_bond<Case::smarts>(mol, bond, ctse::Unique)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_bond_unique<Case::smarts>(mol, bond)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, bond)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, bond, ctse::Unique)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps_unique<Case::smarts>(mol, bond)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    }
}

TEMPLATE_LIST_TEST_CASE("maps_bond_unique", "", MapsBondUniqueCases)
{
    test_maps_bond_unique<Molecule::MockIndexMolecule, TestType>();
    test_maps_bond_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// maps_bond_all(mol, bond)
//

using MapsBondAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC~O", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 0, Maps<Map<0, 1, 2>, Map<0, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 1, Maps<Map<2, 1, 0>, Map<2, 1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*~*", 2, Maps<Map<3, 1, 0>, Map<3, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 0, Maps<Map<0, 1, 2, 3>, Map<0, 1, 3, 2>>>,
    //
    // CCCC
    //
    IndexTypeTestCase<MockButane, "****", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockButane, "****", 1, EmptyMaps>,
    IndexTypeTestCase<MockButane, "****", 2, Maps<Map<3, 2, 1, 0>>>,
    IndexTypeTestCase<MockButane, "*(**)*", 0, EmptyMaps>,
    IndexTypeTestCase<MockButane, "*(**)*", 1, Maps<Map<1, 2, 3, 0>, Map<2, 1, 0, 3>>>,
    IndexTypeTestCase<MockButane, "*(**)*", 2, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_maps_bond_all()
{
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    CHECK_THAT(toVectorMaps(ctse::maps_bond<Case::smarts>(mol, bond, ctse::All)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    CHECK_THAT(toVectorMaps(ctse::maps_bond_all<Case::smarts>(mol, bond)),
               Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK_THAT(toVectorMaps(ctse::maps<Case::smarts>(mol, bond, ctse::All)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
        CHECK_THAT(toVectorMaps(ctse::maps_all<Case::smarts>(mol, bond)),
                   Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
    }
}

TEMPLATE_LIST_TEST_CASE("maps_bond_all", "", MapsBondAllCases)
{
    test_maps_bond_all<Molecule::MockIndexMolecule, TestType>();
    test_maps_bond_all<Molecule::MockProxyMolecule, TestType>();
}
