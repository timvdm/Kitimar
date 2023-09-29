#pragma once

#include <Kitimar/Util/Util.hpp>

#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecules.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#define CASE_INFO(name) \
    INFO(name << "< \"" << Util::toString(Case::smarts) << "\" >( \"" << Case::smiles << "\" )")

#define INDEX_CASE_INFO(name) \
    INFO(name << "< \"" << Util::toString(Case::smarts) << "\" >( \"" << Case::smiles << "\" , " << Case::index << " )")

template<template<typename> class MolFactory, ctll::fixed_string SMARTS>
struct TestCase
{
    constexpr static inline auto smarts = SMARTS;
    constexpr static inline auto smiles = MolFactory<Kitimar::Molecule::MockIndexMolecule>::smiles;

    template<typename Mol>
    static auto mol()
    {
        return MolFactory<Mol>::create();
    }
};

template<template<typename> class MolFactory, ctll::fixed_string SMARTS, auto Expected>
struct ValueTestCase : TestCase<MolFactory, SMARTS>
{
    constexpr static inline auto expected = Expected;
};

template<template<typename> class MolFactory, ctll::fixed_string SMARTS, typename Expected>
struct TypeTestCase : TestCase<MolFactory, SMARTS>
{
    constexpr static inline auto expected = Expected{};
};

template<template<typename> class MolFactory, ctll::fixed_string SMARTS, int Index, auto Expected>
struct IndexValueTestCase : ValueTestCase<MolFactory, SMARTS, Expected>
{
    constexpr static inline auto index = Index;
};

template<template<typename> class MolFactory, ctll::fixed_string SMARTS, int Index, typename Expected>
struct IndexTypeTestCase : TypeTestCase<MolFactory, SMARTS, Expected>
{
    constexpr static inline auto index = Index;
};

template<uint32_t ...Index>
struct TestMap
{
    static constexpr inline bool found = sizeof...(Index) > 0;
    static constexpr inline std::array<uint32_t, sizeof...(Index)> map = {Index...};
};

template<uint32_t ...Index>
using Map = TestMap<Index...>;

using EmptyMap = TestMap<>;

static auto toVectorMaps(const auto &arrays)
{
    std::vector<std::vector<uint32_t>> vectors;
    for (const auto &array : arrays)
        vectors.push_back({std::begin(array), std::end(array)});
    return vectors;
}

static auto toVectorMaps(const auto &arrays, auto &mol)
{
    std::vector<std::vector<uint32_t>> vectors;
    for (const auto &array : arrays) {
        auto map = array | std::views::transform([&mol] (const auto &atom) {
            return get_index(mol, atom);
        });
        vectors.push_back({std::begin(map), std::end(map)});
    }
    return vectors;
}

template<typename ...MapTs>
struct TestMaps
{

    static std::vector<std::vector<uint32_t>> maps()
    {
        if constexpr (sizeof...(MapTs)) {
            using T = std::common_type_t<decltype(MapTs::map)...>;
            std::array<T, sizeof...(MapTs)> arrays = {MapTs::map...};
            return toVectorMaps(arrays);
        } else
            return {};
    }
};

template<typename ...MapTs>
using Maps = TestMaps<MapTs...>;

using EmptyMaps = TestMaps<>;
