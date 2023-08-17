#include "TestCTSmarts.hpp"
#include <Kitimar/Util/Util.hpp>

#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <array>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;

template<typename Case>
void test_captures_impl(auto &mol, const auto &captures)
{
    CHECK_THAT(toVectorMaps(captures, mol), Catch::Matchers::UnorderedRangeEquals(Case::expected.maps()));
}

//
// captures_unique(mol)
//

using CapturesUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Maps<Map<0>, Map<1>>>,
    TypeTestCase<MockAcetateAnion, "O", Maps<Map<2>, Map<3>>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMaps>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Maps<Map<0, 1>>>,
    TypeTestCase<MockAcetateAnion, "C=O", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "CO", Maps<Map<1, 3>>>,

    TypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", Maps<Map<2, 1>>>,
    TypeTestCase<MockAcetateAnion, "[C:1]=O", Maps<Map<1>>>,
    TypeTestCase<MockAcetateAnion, "C=[O:1]", Maps<Map<2>>>,

    TypeTestCase<MockAcetateAnion, "*~*", Maps<Map<0, 1>, Map<1, 2>, Map<1, 3>>>,
    TypeTestCase<MockAcetateAnion, "[*:1]~*", Maps<Map<0>, Map<1>, Map<2>, Map<3>>>,

    TypeTestCase<MockAcetateAnion, "*#*", EmptyMaps>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC=O", Maps<Map<0, 1, 2>>>,
    TypeTestCase<MockAcetateAnion, "O=CC", Maps<Map<2, 1, 0>>>,
    TypeTestCase<MockAcetateAnion, "C(=O)C", Maps<Map<1, 2, 0>>>,
    TypeTestCase<MockAcetateAnion, "CC(=O)O", Maps<Map<0, 1, 2, 3>>>,
    TypeTestCase<MockAcetateAnion, "CC(O)=O", Maps<Map<0, 1, 3, 2>>>,
    TypeTestCase<MockAcetateAnion, "*~*~*", Maps<Map<0, 1, 2>, Map<0, 1, 3>, Map<2, 1, 3>>>,
    TypeTestCase<MockAcetateAnion, "*~*(~*)~*", Maps<Map<0, 1, 2, 3>>>,
    TypeTestCase<MockAcetateAnion, "[*:1]~[*:2]~*", Maps<Map<0, 1>, Map<2, 1>, Map<3, 1>>>,
    TypeTestCase<MockAcetateAnion, "*~*(~*)~*", Maps<Map<0, 1, 2, 3>>>,

    TypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", Maps<Map<0>>>,
    TypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", Maps<Map<0, 2>>>,
    TypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", Maps<Map<3, 0>>>,



    TypeTestCase<MockAcetateAnion, "CC(=O)N", EmptyMaps>
>;

template<typename Mol, typename Case>
void test_captures_unique()
{
    CASE_INFO("captures_unique");
    auto mol = Case::template mol<Mol>();
    test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol));
    test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, ctse::Unique));
    test_captures_impl<Case>(mol, ctse::captures_unique<Case::smarts>(mol));
}

TEMPLATE_LIST_TEST_CASE("captures_unique", "", CapturesUniqueCases)
{
    test_captures_unique<Molecule::MockIndexMolecule, TestType>();
    test_captures_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// captures_all(mol)
//

using CapturesAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Maps<Map<0>, Map<1>>>,
    TypeTestCase<MockAcetateAnion, "O", Maps<Map<2>, Map<3>>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMaps>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Maps<Map<0, 1>, Map<1, 0>>>,
    TypeTestCase<MockAcetateAnion, "C=O", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "CO", Maps<Map<1, 3>>>,

    TypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", Maps<Map<1, 2>>>,
    TypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", Maps<Map<2, 1>>>,
    TypeTestCase<MockAcetateAnion, "[C:1]=O", Maps<Map<1>>>,
    TypeTestCase<MockAcetateAnion, "C=[O:1]", Maps<Map<2>>>,

    TypeTestCase<MockAcetateAnion, "*#*", EmptyMaps>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC=O", Maps<Map<0, 1, 2>>>,
    TypeTestCase<MockAcetateAnion, "O=CC", Maps<Map<2, 1, 0>>>,
    TypeTestCase<MockAcetateAnion, "C(=O)C", Maps<Map<1, 2, 0>>>,
    TypeTestCase<MockAcetateAnion, "CC(=O)O", Maps<Map<0, 1, 2, 3>>>,
    TypeTestCase<MockAcetateAnion, "CC(O)=O", Maps<Map<0, 1, 3, 2>>>,
    TypeTestCase<MockAcetateAnion, "*~*~*", Maps<Map<0, 1, 2>, Map<2, 1, 0>, Map<0, 1, 3>, Map<3, 1, 0>, Map<2, 1, 3>, Map<3, 1, 2>>>,
    TypeTestCase<MockAcetateAnion, "*~*(~*)~*", Maps<Map<0, 1, 2, 3>, Map<0, 1, 3, 2>, Map<2, 1, 0, 3>, Map<2, 1, 3, 0>, Map<3, 1, 0, 2>, Map<3, 1, 2, 0>>>,

    TypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", Maps<Map<0>>>,
    TypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", Maps<Map<0, 2>>>,
    TypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", Maps<Map<3, 0>>>,

    TypeTestCase<MockAcetateAnion, "CC(=O)N", EmptyMaps>
>;


template<typename Mol, typename Case>
void test_captures_all()
{
    CASE_INFO("captures_all");
    auto mol = Case::template mol<Mol>();
    test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, ctse::All));
    test_captures_impl<Case>(mol, ctse::captures_all<Case::smarts>(mol));
}

TEMPLATE_LIST_TEST_CASE("captures_all", "", CapturesAllCases)
{
    test_captures_all<Molecule::MockIndexMolecule, TestType>();
    test_captures_all<Molecule::MockProxyMolecule, TestType>();
}

//
// captures_atom_unique(mol, atom)
//

using CapturesAtomUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Maps<Map<1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 1, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 3, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Maps<Map<1>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Maps<Map<2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMaps>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 2, Maps<Map<2, 1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Maps<Map<1, 2, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 1, Maps<Map<1, 0, 2>, Map<1, 0, 3>, Map<1, 2, 3>>>,


    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Maps<Map<0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Maps<Map<0, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Maps<Map<3, 0>>>
>;



template<typename Mol, typename Case>
void test_captures_atom_unique()
{
    INDEX_CASE_INFO("captures_atom_unique");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    test_captures_impl<Case>(mol, ctse::captures_atom<Case::smarts>(mol, atom));
    test_captures_impl<Case>(mol, ctse::captures_atom<Case::smarts>(mol, atom, ctse::Unique));
    test_captures_impl<Case>(mol, ctse::captures_atom_unique<Case::smarts>(mol, atom));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, atom));
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, atom, ctse::Unique));
        test_captures_impl<Case>(mol, ctse::captures_unique<Case::smarts>(mol, atom));
    }
}

TEMPLATE_LIST_TEST_CASE("captures_atom_unique", "", CapturesAtomUniqueCases)
{
    test_captures_atom_unique<Molecule::MockIndexMolecule, TestType>();
    test_captures_atom_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// captures_atom_all(mol, atom)
//

using CapturesAtomAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Maps<Map<1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 1, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 3, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Maps<Map<1>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Maps<Map<2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMaps>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 2, Maps<Map<2, 1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Maps<Map<1, 2, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 1, Maps<Map<1, 0, 2>, Map<1, 2, 0>, Map<1, 0, 3>, Map<1, 3, 0>, Map<1, 2, 3>, Map<1, 3, 2>>>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Maps<Map<0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Maps<Map<0, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Maps<Map<3, 0>>>
>;



template<typename Mol, typename Case>
void test_captures_atom_all()
{
    INDEX_CASE_INFO("captures_atom_all");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    test_captures_impl<Case>(mol, ctse::captures_atom<Case::smarts>(mol, atom, ctse::All));
    test_captures_impl<Case>(mol, ctse::captures_atom_all<Case::smarts>(mol, atom));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, atom, ctse::All));
        test_captures_impl<Case>(mol, ctse::captures_all<Case::smarts>(mol, atom));
    }
}

TEMPLATE_LIST_TEST_CASE("captures_atom_all", "", CapturesAtomAllCases)
{
    test_captures_atom_all<Molecule::MockIndexMolecule, TestType>();
    test_captures_atom_all<Molecule::MockProxyMolecule, TestType>();
}

//
// captures_bond_unique(mol, bond)
//

using CapturesBondUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 2, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 0, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Maps<Map<1>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Maps<Map<2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMaps>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 1, Maps<Map<2, 1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Maps<Map<1, 2, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 0, Maps<Map<1, 0, 2>, Map<1, 0, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)(~*)~*", 0, Maps<Map<1, 0, 2, 3>>>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Maps<Map<0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Maps<Map<0, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Maps<Map<3, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 1, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "CCO", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "OCC", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O(C)C", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C(O)=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)O", 2, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_captures_bond_unique()
{
    INDEX_CASE_INFO("captures_bond_unique");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    test_captures_impl<Case>(mol, ctse::captures_bond<Case::smarts>(mol, bond));
    test_captures_impl<Case>(mol, ctse::captures_bond<Case::smarts>(mol, bond, ctse::Unique));
    test_captures_impl<Case>(mol, ctse::captures_bond_unique<Case::smarts>(mol, bond));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, bond));
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, bond, ctse::Unique));
        test_captures_impl<Case>(mol, ctse::captures_unique<Case::smarts>(mol, bond));
    }
}

TEMPLATE_LIST_TEST_CASE("captures_bond_unique", "", CapturesBondUniqueCases)
{
    test_captures_bond_unique<Molecule::MockIndexMolecule, TestType>();
    test_captures_bond_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// captures_bond_all(mol, bond)
//

using CapturesBondAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Maps<Map<0, 1>, Map<1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 2, Maps<Map<1, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 0, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Maps<Map<1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Maps<Map<1>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Maps<Map<2>>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMaps>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMaps>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Maps<Map<0, 1, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 1, Maps<Map<2, 1, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Maps<Map<1, 2, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Maps<Map<0, 1, 2, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Maps<Map<0, 1, 3, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)~*", 0, Maps<Map<1, 0, 2>, Map<1, 0, 3>>>,
    IndexTypeTestCase<MockAcetateAnion, "*(~*)(~*)~*", 0, Maps<Map<1, 0, 2, 3>, Map<1, 0, 3, 2>>>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Maps<Map<0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Maps<Map<0, 2>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Maps<Map<3, 0>>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 1, EmptyMaps>,


    IndexTypeTestCase<MockAcetateAnion, "CCO", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "OCC", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O(C)C", 0, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C(O)=O", 1, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 2, EmptyMaps>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)O", 2, EmptyMaps>
>;

template<typename Mol, typename Case>
void test_captures_bond_all()
{
    INDEX_CASE_INFO("captures_bond_all");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    test_captures_impl<Case>(mol, ctse::captures_bond<Case::smarts>(mol, bond, ctse::All));
    test_captures_impl<Case>(mol, ctse::captures_bond_all<Case::smarts>(mol, bond));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        test_captures_impl<Case>(mol, ctse::captures<Case::smarts>(mol, bond, ctse::All));
        test_captures_impl<Case>(mol, ctse::captures_all<Case::smarts>(mol, bond));
    }
}

TEMPLATE_LIST_TEST_CASE("captures_bond_all", "", CapturesBondAllCases)
{
    test_captures_bond_all<Molecule::MockIndexMolecule, TestType>();
    test_captures_bond_all<Molecule::MockProxyMolecule, TestType>();
}
