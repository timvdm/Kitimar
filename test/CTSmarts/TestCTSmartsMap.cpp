#include "TestCTSmarts.hpp"
#include <Kitimar/Util/Util.hpp>

#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <array>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;

template<int ...Indexes>
struct TestMap
{
    static constexpr inline bool found = sizeof...(Indexes);
    static constexpr inline std::array<int, sizeof...(Indexes)> map = {Indexes...};
};

template<int ...Indexes>
using Map = TestMap<Indexes...>;

using EmptyMap = TestMap<>;

//
// map(mol)
//

using MapCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    TypeTestCase<MockAcetateAnion, "C", Map<0>>,
    TypeTestCase<MockAcetateAnion, "O", Map<2>>,
    TypeTestCase<MockAcetateAnion, "[O-]", Map<3>>,

    TypeTestCase<MockAcetateAnion, "N", EmptyMap>,
    TypeTestCase<MockAcetateAnion, "[O+]", EmptyMap>,

    // single bond
    TypeTestCase<MockAcetateAnion, "CC", Map<0, 1>>,
    TypeTestCase<MockAcetateAnion, "C=O", Map<1, 2>>,
    TypeTestCase<MockAcetateAnion, "C[O-]", Map<1, 3>>,

    TypeTestCase<MockAcetateAnion, "C=C", EmptyMap>,
    TypeTestCase<MockAcetateAnion, "C#O", EmptyMap>,
    TypeTestCase<MockAcetateAnion, "C[O+]", EmptyMap>,

    // general case
    TypeTestCase<MockAcetateAnion, "CC(=O)[O-]", Map<0, 1, 2, 3>>,
    TypeTestCase<MockAcetateAnion, "CC([O-])=O", Map<0, 1, 3, 2>>,
    TypeTestCase<MockAcetateAnion, "C(C)(=O)[O-]", Map<1, 0, 2, 3>>,
    TypeTestCase<MockAcetateAnion, "C(C)([O-])=O", Map<1, 0, 3, 2>>,
    TypeTestCase<MockAcetateAnion, "C(=O)(C)[O-]", Map<1, 2, 0, 3>>,
    TypeTestCase<MockAcetateAnion, "C(=O)([O-])C", Map<1, 2, 3, 0>>,
    TypeTestCase<MockAcetateAnion, "C([O-])(C)=O", Map<1, 3, 0, 2>>,
    TypeTestCase<MockAcetateAnion, "C([O-])(=O)C", Map<1, 3, 2, 0>>,
    TypeTestCase<MockAcetateAnion, "O=C(C)[O-]", Map<2, 1, 0, 3>>,
    TypeTestCase<MockAcetateAnion, "O=C([O-])C", Map<2, 1, 3, 0>>,
    TypeTestCase<MockAcetateAnion, "[O-]C(C)=O", Map<3, 1, 0, 2>>,
    TypeTestCase<MockAcetateAnion, "[O-]C(=O)C", Map<3, 1, 2, 0>>,

    TypeTestCase<MockAcetateAnion, "CC(=O)N", EmptyMap>

>;

template<typename Mol, typename Case>
void test_map()
{
    auto mol = Case::template mol<Mol>();
    auto [found, map] = ctse::map<Case::smarts>(mol);
    CHECK(found == Case::expected.found);
    if (found)
        CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));
}

TEMPLATE_LIST_TEST_CASE("map", "", MapCases)
{
    test_map<Molecule::MockIndexMolecule, TestType>();
    test_map<Molecule::MockProxyMolecule, TestType>();
}

//
// map_atom(mol, atom)
//

using MapAtomCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Map<0, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "CC", 1, Map<1, 0>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "O=C", 2, Map<2, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 1, Map<1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "[O-]C", 3, Map<3, 1>>,

    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O=C", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C[O-]", 3, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "[O-]C", 1, EmptyMap>,

    IndexTypeTestCase<MockAcetateAnion, "C=C", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C#O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C[O+]", 1, EmptyMap>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CC(=O)[O-]", 0, Map<0, 1, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "C(C)(=O)[O-]", 1, Map<1, 0, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 2, Map<2, 1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 3, Map<3, 1, 2>>,

    IndexTypeTestCase<MockAcetateAnion, "CC(=O)[O-]", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C(C)(=O)[O-]", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 3, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 2, EmptyMap>
>;

template<typename Mol, typename Case>
void test_map_atom()
{
    INFO("map_atom< \"" << Util::toString(Case::smarts) << "\" >( \"" << Case::smiles << "\" , " << Case::index << " )");
    auto mol = Case::template mol<Mol>();
    auto atom = get_atom(mol, Case::index);

    auto [found, map] = ctse::map_atom<Case::smarts>(mol, atom);
    CHECK(found == Case::expected.found);
    if (found)
        CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        auto [found, map] = ctse::map<Case::smarts>(mol, atom);
        CHECK(found == Case::expected.found);
        if (found)
            CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));
    }
}

TEMPLATE_LIST_TEST_CASE("map_atom", "", MapAtomCases)
{
    test_map_atom<Molecule::MockIndexMolecule, TestType>();
    test_map_atom<Molecule::MockProxyMolecule, TestType>();
}

//
// map_bond(mol, bond)
//

using MapBondCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single bond
    IndexTypeTestCase<MockAcetateAnion, "CC", 0, Map<0, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 1, Map<1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "O=C", 1, Map<2, 1>>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 2, Map<1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "OC", 2, Map<3, 1>>,

    IndexTypeTestCase<MockAcetateAnion, "C=C", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "NC", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "CN", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "CO", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "N=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C=N", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C=O", 2, EmptyMap>,

    // general case
    IndexTypeTestCase<MockAcetateAnion, "CCO", 0, Map<0, 1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "C(C)O", 0, Map<1, 0, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 1, Map<2, 1, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)O", 1, Map<1, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 2, Map<3, 1, 2>>,
    IndexTypeTestCase<MockAcetateAnion, "C(O)=O", 2, Map<1, 3, 2>>,

    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 0, Map<0, 1, 2, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 1, Map<2, 1, 0, 3>>,
    IndexTypeTestCase<MockAcetateAnion, "*~*(~*)~*", 2, Map<3, 1, 0, 2>>,

    IndexTypeTestCase<MockAcetateAnion, "CCO", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "OCC", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O(C)C", 0, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "OC=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C(O)=O", 1, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "O=CO", 2, EmptyMap>,
    IndexTypeTestCase<MockAcetateAnion, "C(=O)O", 2, EmptyMap>
>;

template<typename Mol, typename Case>
void test_map_bond()
{
    INFO("map_bond< \"" << Util::toString(Case::smarts) << "\" >( \"" << Case::smiles << "\" , " << Case::index << " )");
    auto mol = Case::template mol<Mol>();
    auto bond = get_bond(mol, Case::index);

    auto [found, map] = ctse::map_bond<Case::smarts>(mol, bond);
    CHECK(found == Case::expected.found);
    if (found)
        CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));

    if constexpr (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) {
        auto [found, map] = ctse::map<Case::smarts>(mol, bond);
        CHECK(found == Case::expected.found);
        if (found)
            CHECK_THAT(map, Catch::Matchers::RangeEquals(Case::expected.map));
    }
}

TEMPLATE_LIST_TEST_CASE("map_bond", "", MapBondCases)
{
    test_map_bond<Molecule::MockIndexMolecule, TestType>();
    test_map_bond<Molecule::MockProxyMolecule, TestType>();
}
