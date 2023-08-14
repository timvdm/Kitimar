#pragma once

#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Molecule/MockMolecules.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

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


