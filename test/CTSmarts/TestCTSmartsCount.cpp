#include "TestCTSmarts.hpp"

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;
using Molecule::MockButane;

//
// count_unique(mol)
//

using CountUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    ValueTestCase<MockAcetateAnion, "C", 2>,
    ValueTestCase<MockAcetateAnion, "O", 2>,
    ValueTestCase<MockAcetateAnion, "[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "N", 0>,
    ValueTestCase<MockAcetateAnion, "[O+]", 0>,
    // single bond
    ValueTestCase<MockAcetateAnion, "CC", 1>,
    ValueTestCase<MockAcetateAnion, "C=O", 1>,
    ValueTestCase<MockAcetateAnion, "C[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "C=C", 0>,
    ValueTestCase<MockAcetateAnion, "C#O", 0>,
    ValueTestCase<MockAcetateAnion, "C[O+]", 0>,
    // general case
    ValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "CC([O-])=O", 1>,
    ValueTestCase<MockAcetateAnion, "O=C(C)[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "*~*(~*)~*", 1>,
    ValueTestCase<MockAcetateAnion, "CC(=O)N", 0>
>;

template<typename Mol, typename Case>
void test_count_unique()
{
    CASE_INFO("count_unique");
    auto mol = Case::template mol<Mol>();
    CHECK(ctse::count<Case::smarts>(mol) == Case::expected);
    CHECK(ctse::count<Case::smarts>(mol, ctse::Unique) == Case::expected);
    CHECK(ctse::count_unique<Case::smarts>(mol) == Case::expected);
}

TEMPLATE_LIST_TEST_CASE("count_unique", "", CountUniqueCases)
{
    test_count_unique<Molecule::MockIndexMolecule, TestType>();
    test_count_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// count_all(mol)
//

using CountAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    ValueTestCase<MockAcetateAnion, "C", 2>,
    ValueTestCase<MockAcetateAnion, "O", 2>,
    ValueTestCase<MockAcetateAnion, "[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "N", 0>,
    ValueTestCase<MockAcetateAnion, "[O+]", 0>,
    // single bond
    ValueTestCase<MockAcetateAnion, "CC", 2>,
    ValueTestCase<MockAcetateAnion, "C=O", 1>,
    ValueTestCase<MockAcetateAnion, "C[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "C=C", 0>,
    ValueTestCase<MockAcetateAnion, "C#O", 0>,
    ValueTestCase<MockAcetateAnion, "C[O+]", 0>,
    // general case
    ValueTestCase<MockAcetateAnion, "CC(~O)~O", 2>,
    ValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "CC([O-])=O", 1>,
    ValueTestCase<MockAcetateAnion, "O=C(C)[O-]", 1>,
    ValueTestCase<MockAcetateAnion, "*~*(~*)~*", 6>,
    ValueTestCase<MockAcetateAnion, "CC(=O)N", 0>
>;

template<typename Mol, typename Case>
void test_count_all()
{
    CASE_INFO("count_all");
    auto mol = Case::template mol<Mol>();
    CHECK(ctse::count<Case::smarts>(mol, ctse::All) == Case::expected);
    CHECK(ctse::count_all<Case::smarts>(mol) == Case::expected);
}

TEMPLATE_LIST_TEST_CASE("count_all", "", CountAllCases)
{
    test_count_all<Molecule::MockIndexMolecule, TestType>();
    test_count_all<Molecule::MockProxyMolecule, TestType>();
}

//
// count_atom_unique(mol, atom)
//

using CountAtomUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexValueTestCase<MockAcetateAnion, "CC", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*", 1, 3>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "C[O-]", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C[O-]", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "C=C", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C#O", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C[O+]", 1, 0>,
    // general case
    IndexValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC([O-])=O", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "O=C(C)[O-]", 2, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*(~*)~*", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC(=O)N", 0, 0>
>;

template<typename Mol, typename Case>
void test_count_atom_unique()
{
    INDEX_CASE_INFO("count_atom_unique");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);
    CHECK(ctse::count_atom<Case::smarts>(mol, atom) == Case::expected);
    CHECK(ctse::count_atom<Case::smarts>(mol, atom, ctse::Unique) == Case::expected);
    CHECK(ctse::count_atom_unique<Case::smarts>(mol, atom) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK(ctse::count<Case::smarts>(mol, atom) == Case::expected);
        CHECK(ctse::count<Case::smarts>(mol, atom, ctse::Unique) == Case::expected);
        CHECK(ctse::count_unique<Case::smarts>(mol, atom) == Case::expected);
    }
}

TEMPLATE_LIST_TEST_CASE("count_atom_unique", "", CountAtomUniqueCases)
{
    test_count_atom_unique<Molecule::MockIndexMolecule, TestType>();
    test_count_atom_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// count_atom_all(mol, atom)
//

using CountAtomAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexValueTestCase<MockAcetateAnion, "CC", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*", 1, 3>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "C[O-]", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C[O-]", 1, 1>,
    IndexValueTestCase<MockAcetateAnion, "C=C", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C#O", 0, 0>,
    IndexValueTestCase<MockAcetateAnion, "C[O+]", 1, 0>,
    // general case
    IndexValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC([O-])=O", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "O=C(C)[O-]", 2, 1>,
    IndexValueTestCase<MockAcetateAnion, "*~*(~*)~*", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "CC(=O)N", 0, 0>
>;

template<typename Mol, typename Case>
void test_count_atom_all()
{
    INDEX_CASE_INFO("count_atom_all");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);
    CHECK(ctse::count_atom<Case::smarts>(mol, atom, ctse::All) == Case::expected);
    CHECK(ctse::count_atom_all<Case::smarts>(mol, atom) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK(ctse::count<Case::smarts>(mol, atom, ctse::All) == Case::expected);
        CHECK(ctse::count_all<Case::smarts>(mol, atom) == Case::expected);
    }
}

TEMPLATE_LIST_TEST_CASE("count_atom_all", "", CountAtomAllCases)
{
    test_count_atom_all<Molecule::MockIndexMolecule, TestType>();
    test_count_atom_all<Molecule::MockProxyMolecule, TestType>();
}

//
// count_bond_unique(mol, bond)
//

using CountBondUniqueCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    IndexValueTestCase<MockAcetateAnion, "CC=O", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC~O", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "CC=O", 1, 0>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 1, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 2, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*(~*)~*", 0, 1>,
    //
    // CCCC
    //
    IndexValueTestCase<MockButane, "****", 0, 1>,
    IndexValueTestCase<MockButane, "****", 1, 0>,
    IndexValueTestCase<MockButane, "****", 2, 1>,
    IndexValueTestCase<MockButane, "*(**)*", 0, 0>,
    IndexValueTestCase<MockButane, "*(**)*", 1, 1>,
    IndexValueTestCase<MockButane, "*(**)*", 2, 0>
>;

template<typename Mol, typename Case>
void test_count_bond_unique()
{
    INDEX_CASE_INFO("count_bond_unique");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);
    CHECK(ctse::count_bond<Case::smarts>(mol, bond) == Case::expected);
    CHECK(ctse::count_bond<Case::smarts>(mol, bond, ctse::Unique) == Case::expected);
    CHECK(ctse::count_bond_unique<Case::smarts>(mol, bond) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK(ctse::count<Case::smarts>(mol, bond) == Case::expected);
        CHECK(ctse::count<Case::smarts>(mol, bond, ctse::Unique) == Case::expected);
        CHECK(ctse::count_unique<Case::smarts>(mol, bond) == Case::expected);
    }
}

TEMPLATE_LIST_TEST_CASE("count_bond_unique", "", CountBondUniqueCases)
{
    test_count_bond_unique<Molecule::MockIndexMolecule, TestType>();
    test_count_bond_unique<Molecule::MockProxyMolecule, TestType>();
}

//
// count_bond_all(mol, bond)
//

using CountBondAllCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    IndexValueTestCase<MockAcetateAnion, "CC=O", 0, 1>,
    IndexValueTestCase<MockAcetateAnion, "CC~O", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "CC=O", 1, 0>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 0, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 1, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*~*", 2, 2>,
    IndexValueTestCase<MockAcetateAnion, "*~*(~*)~*", 0, 2>,
    //
    // CCCC
    //
    IndexValueTestCase<MockButane, "****", 0, 1>,
    IndexValueTestCase<MockButane, "****", 1, 0>,
    IndexValueTestCase<MockButane, "****", 2, 1>,
    IndexValueTestCase<MockButane, "*(**)*", 0, 0>,
    IndexValueTestCase<MockButane, "*(**)*", 1, 2>,
    IndexValueTestCase<MockButane, "*(**)*", 2, 0>
>;

template<typename Mol, typename Case>
void test_count_bond_all()
{
    INDEX_CASE_INFO("count_bond_all");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);
    CHECK(ctse::count_bond<Case::smarts>(mol, bond, ctse::All) == Case::expected);
    CHECK(ctse::count_bond_all<Case::smarts>(mol, bond) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        CHECK(ctse::count<Case::smarts>(mol, bond, ctse::All) == Case::expected);
        CHECK(ctse::count_all<Case::smarts>(mol, bond) == Case::expected);
    }
}

TEMPLATE_LIST_TEST_CASE("count_bond_all", "", CountBondAllCases)
{
    test_count_bond_all<Molecule::MockIndexMolecule, TestType>();
    test_count_bond_all<Molecule::MockProxyMolecule, TestType>();
}
