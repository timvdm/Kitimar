#pragma once

#include <Kitimar/CTSmarts/CTSmarts.hpp>

#include <iostream>

using namespace Kitimar;
using Kitimar::Molecule::get_nbr;

//using Atom = Molecule::MoleculeTraits<Mol>::Atom;
//using Bond = Molecule::MoleculeTraits<Mol>::Atom;

//
// match_atom(mol, atom)
//

// single atom

template<typename Mol>
bool match_atom_single_atom_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return get_element(mol, atom) == 6;
}

template<typename Mol>
bool match_atom_single_atom_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match_atom<"[#6]">(mol, atom);
}

// single bond

template<typename Mol>
bool match_atom_single_bond_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    if (get_element(mol, atom) != 6)
        return false;
    for (auto bond : get_bonds(mol, atom))
        if (get_element(mol, get_nbr(mol, bond, atom)) == 7)
            return true;
    return false;
}

template<typename Mol>
bool match_atom_single_bond_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match_atom<"[#6]~[#7]">(mol, atom);
}

// 3 atoms chain

template<typename Mol>
bool match_atom_chain_3_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    if (get_element(mol, atom) != 6)
        return false;
    for (auto bond1 : get_bonds(mol, atom)) {
        auto nbr1 = get_nbr(mol, bond1, atom);
        if (get_element(mol, nbr1) != 7)
            continue;
        for (auto bond2 : get_bonds(mol, nbr1)) {
            auto nbr2 = get_nbr(mol, bond2, nbr1);
            if (nbr2 != atom && get_element(mol, nbr2) == 8)
                return true;
        }
    }
    return false;
}

template<typename Mol>
bool match_atom_chain_3_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match_atom<"[#6]~[#7]~[#8]">(mol, atom);
}

template<typename Mol>
bool match_atom_chain_4_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    if (get_element(mol, atom) != 6)
        return false;
    for (auto bond1 : get_bonds(mol, atom)) {
        auto nbr1 = get_nbr(mol, bond1, atom);
        if (get_element(mol, nbr1) != 7)
            continue;
        for (auto bond2 : get_bonds(mol, nbr1)) {
            auto nbr2 = get_nbr(mol, bond2, nbr1);
            if (nbr2 == atom || get_element(mol, nbr2) != 8)
                continue;
            for (auto bond3 : get_bonds(mol, nbr2)) {
                auto nbr3 = get_nbr(mol, bond3, nbr2);
                if (nbr3 != atom && nbr3 != nbr1 && get_element(mol, nbr3) == 9)
                    return true;
            }
        }
    }
    return false;
}

template<typename Mol>
bool match_atom_chain_4_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match_atom<"[#6]~[#7]~[#8]~[#9]">(mol, atom);
}


//
// match(mol)
//

// single atom

template<typename Mol>
bool match_single_atom_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    for (auto atom : get_atoms(mol))
        if (get_element(mol, atom) == 6)
            return true;
    return false;
}

template<typename Mol>
bool match_single_atom_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match<"[#6]">(mol);
}

// single bond

template<typename Mol>
bool match_single_bond_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    for (auto bond : get_bonds(mol)) {
        auto source = get_source(mol, bond);
        auto target = get_target(mol, bond);
        if ((get_element(mol, source) == 6 && get_element(mol, target) == 7) ||
                (get_element(mol, target) == 6 && get_element(mol, source) == 7))
            return true;
    }
    return false;
}

template<typename Mol>
bool match_single_bond_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    return ctse::match<"[#6]~[#7]">(mol);
}

// 3 atoms chain

template<typename Mol>
bool match_chain_3_cpp(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    //int elements[] = {8, 7, 6}; // optimized
    //int elements[] = {6, 7, 8}; // not optimized
    int elements[] = {6, 6, 6}; // not optimized

    for (auto atom : get_atoms(mol)) {
        if (get_element(mol, atom) != elements[0])
            continue;
        for (auto bond1 : get_bonds(mol, atom)) {
            auto nbr = get_nbr(mol, bond1, atom);
            if (get_element(mol, nbr) == elements[1]) {
                for (auto bond2 : get_bonds(mol, nbr)) {
                    auto nbr2 = get_nbr(mol, bond2, nbr);
                    if (nbr2 != atom && get_element(mol, nbr2) == elements[2])
                        return true;
                }
            }
        }
    }
    return false;
}

template<typename Mol>
bool match_chain_3_ctse(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom = {}, typename Molecule::MoleculeTraits<Mol>::Bond = {})
{
    //return ctse::match<"[#6]~[#7]~[#8]">(mol);
    return ctse::match<"[#6]~[#6]~[#6]">(mol);
}

// 4 atoms chain

bool extend_chain(auto &mol, const auto &atom, auto &visited, int i = 1)
{
    if (i == visited.size())
        return true;
    for (auto bond : get_bonds(mol, atom)) {
        auto nbr = get_nbr(mol, bond, atom);
        auto nbrIndex = get_index(mol, nbr);
        if (std::find(visited.begin(), visited.begin() + i, nbrIndex) != visited.begin() + i)
            continue;
        if (get_element(mol, nbr) != 6)
            continue;
        visited[i] = nbrIndex;
        if (extend_chain(mol, nbr, visited, i + 1))
            return true;
        visited[i] = -1;
    }
    return false;
}

template<int N>
bool match_chain_n_cpp(auto &mol)
{
    std::array<decltype(get_index(mol, get_atom(mol, 0))), N/* - 2*/> visited;
    std::ranges::fill(visited, -1);

    for (auto atom : get_atoms(mol)) {
        if (get_element(mol, atom) != 6)
            continue;
        visited[0] = get_index(mol, atom);
        if (extend_chain(mol, atom, visited))
            return true;
        visited[0] = -1;
    }
    return false;
}


