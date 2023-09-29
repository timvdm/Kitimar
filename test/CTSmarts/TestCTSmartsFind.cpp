#include "TestCTSmarts.hpp"
#include <Kitimar/Util/Util.hpp>

#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <array>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

using Molecule::MockAcetateAnion;

/*
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
*/

//
// capture(mol)
//

constexpr auto nullIndex = static_cast<uint32_t>(-1);

using FindCases = std::tuple<
    //
    // CC(=O)[O-]
    //
    // single atom
    ValueTestCase<MockAcetateAnion, "C", 0U>,
    ValueTestCase<MockAcetateAnion, "O", 2U>,

    ValueTestCase<MockAcetateAnion, "N", nullIndex>,

    // single bond
    ValueTestCase<MockAcetateAnion, "CC", 0>,
    ValueTestCase<MockAcetateAnion, "C=O", 1>,
    ValueTestCase<MockAcetateAnion, "CO", 1>,
    ValueTestCase<MockAcetateAnion, "O=C", 2>,
    ValueTestCase<MockAcetateAnion, "OC", 3>,

    ValueTestCase<MockAcetateAnion, "[C:1]C", 0>,
    ValueTestCase<MockAcetateAnion, "[C:1]=O", 1>,
    ValueTestCase<MockAcetateAnion, "[C:1]O", 1>,
    ValueTestCase<MockAcetateAnion, "[O:1]=C", 2>,
    ValueTestCase<MockAcetateAnion, "[O:1]C", 3>,

    ValueTestCase<MockAcetateAnion, "C[C:1]", 1>,
    ValueTestCase<MockAcetateAnion, "C=[O:1]", 2>,
    ValueTestCase<MockAcetateAnion, "C[O:1]", 3>,
    ValueTestCase<MockAcetateAnion, "O=[C:1]", 1>,
    ValueTestCase<MockAcetateAnion, "O[C:1]", 1>,

    ValueTestCase<MockAcetateAnion, "C#C", nullIndex>,

    // general case
    ValueTestCase<MockAcetateAnion, "CCO", 0>,
    ValueTestCase<MockAcetateAnion, "CO", 1>,
    ValueTestCase<MockAcetateAnion, "OCC", 3>,
    ValueTestCase<MockAcetateAnion, "O=CC", 2>,
    ValueTestCase<MockAcetateAnion, "CC(=O)O", 0>,
    ValueTestCase<MockAcetateAnion, "C(C)(=O)O", 1>,

    ValueTestCase<MockAcetateAnion, "[C:1]C(=O)O", 0>,
    ValueTestCase<MockAcetateAnion, "C[C:1](=O)O", 1>,
    ValueTestCase<MockAcetateAnion, "CC(=[O:1])O", 2>,
    ValueTestCase<MockAcetateAnion, "CC(=O)[O:1]", 3>,

    ValueTestCase<MockAcetateAnion, "CC(=O)N", nullIndex>
>;


template<typename Mol, typename Case>
void test_find()
{
    CASE_INFO("find");
    auto mol = Case::template mol<Mol>();
    CHECK(ctse::find<Case::smarts>(mol) == (Case::expected == nullIndex ? null_atom(mol) : get_atom(mol, Case::expected)));
}

TEMPLATE_LIST_TEST_CASE("find", "", FindCases)
{
    test_find<Molecule::MockIndexMolecule, TestType>();
    test_find<Molecule::MockProxyMolecule, TestType>();
}

