#include "TestCTSmarts.hpp"
#include <Kitimar/Util/Util.hpp>

#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <array>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;

template<int I = 0, typename ...Ts>
auto get_capture_map(auto &mol, const std::tuple<bool, Ts...> &tuple)
{
    if constexpr (I == sizeof...(Ts))
        return std::array<uint32_t, sizeof...(Ts)>{};
    else {
        auto map = get_capture_map<I + 1>(mol, tuple);
        map[I] = get_index(mol, std::get<I + 1>(tuple));
        return map;
    }
}

template<typename Case>
void test_capture_impl(auto &mol, const auto &tuple)
{
    auto found = std::get<0>(tuple);
    auto map = get_capture_map<0>(mol, tuple);
    CHECK(found == Case::expected.found);
    if (found)
        CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));
    else
        CHECK(std::ranges::count(map, -1) == map.size());
}

//
// capture(mol)
//

using CaptureCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Map<0>>,
    TypeTestCase<MockAcetateAnion, "O", Map<2>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMap>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Map<0, 1>>,
    TypeTestCase<MockAcetateAnion, "C=O", Map<1, 2>>,
    TypeTestCase<MockAcetateAnion, "CO", Map<1, 3>>,

    TypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", Map<1, 2>>,
    TypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", Map<2, 1>>,
    TypeTestCase<MockAcetateAnion, "[C:1]=O", Map<1>>,
    TypeTestCase<MockAcetateAnion, "C=[O:1]", Map<2>>,

    TypeTestCase<MockAcetateAnion, "*#*", EmptyMap>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC=O", Map<0, 1, 2>>,
    TypeTestCase<MockAcetateAnion, "O=CC", Map<2, 1, 0>>,
    TypeTestCase<MockAcetateAnion, "C(=O)C", Map<1, 2, 0>>,
    TypeTestCase<MockAcetateAnion, "CC(=O)O", Map<0, 1, 2, 3>>,
    TypeTestCase<MockAcetateAnion, "CC(O)=O", Map<0, 1, 3, 2>>,

    TypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", Map<0>>,
    TypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", Map<0, 2>>,
    TypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", Map<3, 0>>
>;


template<typename Mol, typename Case>
void test_capture()
{
    CASE_INFO("capture");
    auto mol = Case::template mol<Mol>();
    test_capture_impl<Case>(mol, ctse::capture<Case::smarts>(mol));
}

TEMPLATE_LIST_TEST_CASE("capture", "", CaptureCases)
{
    test_capture<Molecule::MockIndexMolecule, TestType>();
    test_capture<Molecule::MockProxyMolecule, TestType>();
}

//
// capture_atom(mol, atom)
//

using CaptureAtomCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Map<0, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Map<1, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 1, Map<1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 3, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Map<1>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Map<2>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMap>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Map<0, 1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 2, Map<2, 1, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Map<1, 2, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Map<0, 1, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Map<0, 1, 3, 2>>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Map<0>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Map<0, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Map<3, 0>>
>;



template<typename Mol, typename Case>
void test_capture_atom()
{
    INDEX_CASE_INFO("capture_atom");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    test_capture_impl<Case>(mol, ctse::capture_atom<Case::smarts>(mol, atom));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
            test_capture_impl<Case>(mol, ctse::capture<Case::smarts>(mol, atom));
}

TEMPLATE_LIST_TEST_CASE("capture_atom", "", CaptureAtomCases)
{
    test_capture_atom<Molecule::MockIndexMolecule, TestType>();
    test_capture_atom<Molecule::MockProxyMolecule, TestType>();
}

//
// capture_bond(mol, bond)
//

using CaptureBondCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Map<0, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 2, Map<1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 0, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]=[O:2]", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]=[O:1]", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 1, Map<1>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]=O", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 1, Map<2>>,
    IndexTypeTestCase<MockAcetateAnion, "C=[O:1]", 2, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "*#*", 0, EmptyMap>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 0, Map<0, 1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "CC=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O=CC", 1, Map<2, 1, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)C", 1, Map<1, 2, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)O", 0, Map<0, 1, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "CC(O)=O", 0, Map<0, 1, 3, 2>>,

    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0, Map<0>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:1]C(=[O:2])O", 0, Map<0, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 0, Map<3, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "[C:2]C(=O)[O:1]", 1, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "CCO", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "OCC", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O(C)C", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C(O)=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)O", 2, EmptyMap>
>;

template<typename Mol, typename Case>
void test_capture_bond()
{
    INDEX_CASE_INFO("capture_bond");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    test_capture_impl<Case>(mol, ctse::capture_bond<Case::smarts>(mol, bond));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
        test_capture_impl<Case>(mol, ctse::capture<Case::smarts>(mol, bond));
}

TEMPLATE_LIST_TEST_CASE("capture_bond", "", CaptureBondCases)
{
    test_capture_bond<Molecule::MockIndexMolecule, TestType>();
    test_capture_bond<Molecule::MockProxyMolecule, TestType>();
}
