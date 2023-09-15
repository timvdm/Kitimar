#pragma once

#include "UnconditionalFilter.hpp"

#include "../../Molecule/Molecule.hpp"

namespace Kitimar::CTSmarts {

    struct NumAtomBondFilter : UnconditionalFilter
    {
        static constexpr bool reject(auto smarts, Molecule::Molecule auto &mol) noexcept
        {
            return num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds;
        }
    };

} // namespace Kitimar::CTSmarts
