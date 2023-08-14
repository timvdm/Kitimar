#include "TestCTSmarts.hpp"

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;

//
// match(mol)
//

using MatchCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    ValueTestCase<MockAcetateAnion, "C", true>,
    ValueTestCase<MockAcetateAnion, "O", true>,
    ValueTestCase<MockAcetateAnion, "[O-]", true>,

    ValueTestCase<MockAcetateAnion, "N", false>,
    ValueTestCase<MockAcetateAnion, "[O+]", false>,

    // single bond
    ValueTestCase<MockAcetateAnion, "CC", true>,
    ValueTestCase<MockAcetateAnion, "C=O", true>,
    ValueTestCase<MockAcetateAnion, "C[O-]", true>,

    ValueTestCase<MockAcetateAnion, "C=C", false>,
    ValueTestCase<MockAcetateAnion, "C#O", false>,
    ValueTestCase<MockAcetateAnion, "C[O+]", false>,

    // general case
    ValueTestCase<MockAcetateAnion, "CC(=O)[O-]", true>,
    ValueTestCase<MockAcetateAnion, "CC([O-])=O", true>,
    ValueTestCase<MockAcetateAnion, "O=C(C)[O-]", true>,

    ValueTestCase<MockAcetateAnion, "CC(=O)N", false>
>;

template<typename Mol, typename Case>
void test_match()
{
    auto mol = Case::template mol<Mol>();
    CHECK(ctse::match<Case::smarts>(mol) == Case::expected);
}

TEMPLATE_LIST_TEST_CASE("match", "", MatchCases)
{
    test_match<Molecule::MockIndexMolecule, TestType>();
    test_match<Molecule::MockProxyMolecule, TestType>();
}

//
// match_atom(mol, atom)
//

using MatchAtomCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    IndexValueTestCase<MockAcetateAnion, "C", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "C", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "O", 2, true>,
    IndexValueTestCase<MockAcetateAnion, "[O-]", 3, true>,

    IndexValueTestCase<MockAcetateAnion, "C", 2, false>,
    IndexValueTestCase<MockAcetateAnion, "[O-]", 2, false>,

    // single bond
    IndexValueTestCase<MockAcetateAnion, "CC", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "CC", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "CO", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "O=C", 2, true>,
    IndexValueTestCase<MockAcetateAnion, "OC", 3, true>,

    IndexValueTestCase<MockAcetateAnion, "CO", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "CO", 2, false>,
    IndexValueTestCase<MockAcetateAnion, "OC", 1, false>,

    // general case
    IndexValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "C(C)(=O)[O-]", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "O=CO", 2, true>,
    IndexValueTestCase<MockAcetateAnion, "OC=O", 3, true>,

    IndexValueTestCase<MockAcetateAnion, "CC(=O)[O-]", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "C(C)(=O)[O-]", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "O=CO", 3, false>,
    IndexValueTestCase<MockAcetateAnion, "OC=O", 2, false>
>;

template<typename Mol, typename Case>
void test_match_atom()
{
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);
    CHECK(ctse::match_atom<Case::smarts>(mol, atom) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
        CHECK(ctse::match<Case::smarts>(mol, atom) == Case::expected);
}

TEMPLATE_LIST_TEST_CASE("match_atom", "", MatchAtomCases)
{
    test_match_atom<Molecule::MockIndexMolecule, TestType>();
    test_match_atom<Molecule::MockProxyMolecule, TestType>();
}

//
// match_bond(mol, bond)
//

using MatchBondCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexValueTestCase<MockAcetateAnion, "CC", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "O=C", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "CO", 2, true>,
    IndexValueTestCase<MockAcetateAnion, "OC", 2, true>,

    IndexValueTestCase<MockAcetateAnion, "C=C", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "NC", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "CN", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "CO", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "N=O", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "C=N", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "C=O", 2, false>,

    // general case
    IndexValueTestCase<MockAcetateAnion, "CCO", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "C(C)O", 0, true>,
    IndexValueTestCase<MockAcetateAnion, "O=CO", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "C(=O)O", 1, true>,
    IndexValueTestCase<MockAcetateAnion, "OC=O", 2, true>,
    IndexValueTestCase<MockAcetateAnion, "C(O)=O", 2, true>,

    IndexValueTestCase<MockAcetateAnion, "CCO", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "OCC", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "O(C)C", 0, false>,
    IndexValueTestCase<MockAcetateAnion, "OC=O", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "C(O)=O", 1, false>,
    IndexValueTestCase<MockAcetateAnion, "O=CO", 2, false>,
    IndexValueTestCase<MockAcetateAnion, "C(=O)O", 2, false>
>;

template<typename Mol, typename Case>
void test_match_bond()
{
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);
    CHECK(ctse::match_bond<Case::smarts>(mol, bond) == Case::expected);
    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
        CHECK(ctse::match<Case::smarts>(mol, bond) == Case::expected);
}

TEMPLATE_LIST_TEST_CASE("match_bond", "", MatchBondCases)
{
    test_match_bond<Molecule::MockIndexMolecule, TestType>();
    test_match_bond<Molecule::MockProxyMolecule, TestType>();
}
