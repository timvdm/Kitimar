#pragma once

#include <Kitimar/Molecule/Molecule.hpp>
#include <Kitimar/CTLayout/CTLayout.hpp>

#include <mio/mmap.hpp>

#include <vector>
#include <fstream>

namespace Kitimar {

    struct NumAtoms : ArraySize {};
    struct NumBonds : ArraySize {};

    // Object-Oriented design using Struct
    struct ElementValue : Value<uint8_t> {};
    struct IsotopeValue : Value<uint8_t> {};
    struct ChargeValue : Value<int8_t> {};
    struct DegreeValue : Value<uint8_t> {};
    struct ImplicitHydrogensValue : Value<uint8_t> {};
    struct AromaticAtomValue : Value<bool> {};
    struct AtomStruct : Struct<ElementValue, IsotopeValue, ChargeValue, DegreeValue, ImplicitHydrogensValue, AromaticAtomValue> {};
    struct AtomArray : Array<AtomStruct, NumAtoms> {};

    struct SourceValue : Value<SizeT> {};
    struct TargetValue : Value<SizeT> {};
    struct OrderValue : Value<uint8_t> {};
    struct CyclicBondValue : Value<bool> {};
    struct AromaticBondValue : Value<bool> {};
    struct BondStruct : Struct<SourceValue, TargetValue, OrderValue, AromaticBondValue, CyclicBondValue> {};
    struct BondArray : Array<BondStruct, NumBonds> {};

    // Data-Oriented design using no Struct
    struct ElementArray : Array<ElementValue, NumAtoms> {};
    struct IsotopeArray : Array<IsotopeValue, NumAtoms> {};
    struct ChargeArray : Array<ChargeValue, NumAtoms> {};
    struct DegreeArray : Array<DegreeValue, NumAtoms> {};    
    struct ImplicitHydrogensArray : Array<ImplicitHydrogensValue, NumAtoms> {};
    struct AromaticAtomArray : BitArray<NumAtoms> {};

    struct SourceArray : Array<SourceValue, NumBonds> {};
    struct TargetArray : Array<TargetValue, NumBonds> {};
    struct OrderArray : Array<OrderValue, NumBonds> {};
    struct AromaticBondArray : BitArray<NumBonds> {};

    // incident list (i.e. bonds for atom)
    struct IncidentBegin : Value<SizeT> {};
    struct IncidentBeginArray : Array<IncidentBegin, NumAtoms> {};
    struct IncidentIndex : Value<SizeT> {};
    struct IncidentArray : Array<IncidentIndex, NumBonds, 2> {};


    // adjacent list (i.e. nbrs of atom)
    struct AdjacentBegin : Value<SizeT> {};
    struct AdjacentBeginArray : Array<AdjacentBegin, NumAtoms> {};
    struct AdjacentIndex : Value<SizeT> {};
    struct AdjacentArray : Array<AdjacentIndex, NumBonds, 2> {};



    struct StructMolecule : Layout<AtomArray, BondArray> {};
    struct StructMoleculeIncident : Layout<AtomArray, BondArray, IncidentBeginArray, IncidentArray> {};
    struct StructMoleculeAdjacent : Layout<AtomArray, BondArray, AdjacentBeginArray, AdjacentArray> {};
    struct StructMoleculeIncidentAdjacent : Layout<AtomArray, BondArray,
                                                   IncidentBeginArray, IncidentArray,
                                                   AdjacentBeginArray, AdjacentArray> {};

    struct ArrayMolecule : Layout<ElementArray, IsotopeArray, ChargeArray, DegreeArray, ImplicitHydrogensArray,
                                  SourceArray, TargetArray, OrderArray> {};
    struct ArrayMoleculeIncident : Layout<ElementArray, IsotopeArray, ChargeArray, DegreeArray, ImplicitHydrogensArray,
                                          SourceArray, TargetArray, OrderArray,
                                          IncidentBeginArray, IncidentArray> {};
    struct ArrayMoleculeAdjacent : Layout<ElementArray, IsotopeArray, ChargeArray, DegreeArray, ImplicitHydrogensArray,
                                          SourceArray, TargetArray, OrderArray,
                                          AdjacentBeginArray, AdjacentArray> {};
    struct ArrayMoleculeIncidentAdjacent : Layout<ElementArray, IsotopeArray, ChargeArray, DegreeArray, ImplicitHydrogensArray,
                                                  SourceArray, TargetArray, OrderArray,
                                                  IncidentBeginArray, IncidentArray,
                                                  AdjacentBeginArray, AdjacentArray> {};

    template<typename Layout>
    constexpr auto moleculeSize(SizeT numAtoms, SizeT numBonds) noexcept
    {
        return Layout::size({numAtoms, numBonds});
    }

    //
    // Molecule
    //

    template<typename MolObj>
    constexpr auto num_atoms(const MolObj &mol) noexcept
    {
        if constexpr (ctll::exists_in(AtomArray{}, MolObj::Layout::type))
            return mol.template size<AtomArray>();
        else
            return mol.template size<ElementArray>();
    }

    template<typename MolObj>
    constexpr auto num_bonds(const MolObj &mol) noexcept
    {
        if constexpr (ctll::exists_in(BondArray{}, MolObj::Layout::type))
            return mol.template size<BondArray>();
        else
            return mol.template size<SourceArray>();
    }

    constexpr auto get_atoms(const auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_atoms(mol))>(0), num_atoms(mol));
    }

    constexpr auto get_bonds(const auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_bonds(mol))>(0), num_bonds(mol));
    }

    constexpr auto get_atom(const auto &mol, SizeT index) noexcept
    {
        assert(index < num_atoms(mol));
        return index;
    }

    constexpr auto get_bond(const auto &mol, SizeT index) noexcept
    {
        assert(index < num_bonds(mol));
        return index;
    }

    constexpr auto get_index(const auto &mol, SizeT index) noexcept
    {
        return index;
    }


    // Atom

    constexpr auto get_element(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<ElementValue>(atom);
    }

    constexpr auto get_isotope(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<IsotopeValue>(atom);
    }

    constexpr auto get_charge(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<ChargeValue>(atom);
    }

    constexpr auto get_degree(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<DegreeValue>(atom);
    }

    constexpr auto get_implicit_hydrogens(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<ImplicitHydrogensValue>(atom);
    }

    constexpr auto is_aromatic_atom(const auto &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.template value<AromaticAtomValue>(atom);
    }

    template<typename MolObj>
    constexpr auto get_bonds(const MolObj &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (ctll::exists_in(IncidentArray{}, MolObj::Layout::type)) {
            static_assert(ctll::exists_in(IncidentBeginArray{}, MolObj::Layout::type));
            auto begin = mol.template value<IncidentBegin>(atom);
            auto offset = mol.template offset<IncidentIndex>(begin);
            auto ptr = reinterpret_cast<const SizeT*>(mol.data() + offset);
            return std::ranges::subrange(ptr, ptr + get_degree(mol, atom));
        } else {
            return get_bonds(mol) | std::views::filter([&mol, atom] (auto i) {
                auto bond = get_bond(mol, i);
                return atom == get_source(mol, bond) || atom == get_target(mol, bond);
            });
        }
    }

    template<typename MolObj>
    constexpr auto get_nbrs(const MolObj &mol, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (ctll::exists_in(AdjacentArray{}, MolObj::Layout::type)) {
            static_assert(ctll::exists_in(AdjacentBeginArray{}, MolObj::Layout::type));
            auto begin = mol.template value<AdjacentBegin>(atom);
            auto offset = mol.template offset<AdjacentIndex>(begin);
            auto ptr = reinterpret_cast<const SizeT*>(mol.data() + offset);
            return std::ranges::subrange(ptr, ptr + get_degree(mol, atom));
        } else {
            return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
                return get_nbr(mol, bond, atom);
            });
        }
    }



    //
    // Bond
    //

    constexpr auto get_source(const auto &mol, SizeT bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.template value<SourceValue>(bond);
    }

    constexpr auto get_target(const auto &mol, SizeT bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.template value<TargetValue>(bond);
    }

    constexpr auto get_order(const auto &mol, SizeT bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.template value<OrderValue>(bond);
    }

    constexpr auto is_cyclic_bond(const auto &mol, SizeT bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.template value<CyclicBondValue>(bond);
    }

    constexpr auto is_aromatic_bond(const auto &mol, SizeT bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.template value<AromaticBondValue>(bond);
    }

    constexpr auto get_nbr(const auto &mol, SizeT bond, SizeT atom) noexcept
    {
        assert(atom < num_atoms(mol));
        assert(bond < num_bonds(mol));
        auto source = get_source(mol, bond);
        return source != atom ? source : get_target(mol, bond);
    }


    template<typename Layout>
    void serialize(Molecule auto &mol, std::byte *data)
    {
        auto obj = EditableObject<Layout>{data};
        obj.setSize(NumAtoms{}, num_atoms(mol));
        obj.setSize(NumBonds{}, num_bonds(mol));
        obj.setSize(obj.size());

        for (auto atom : get_atoms(mol)) {
            auto index = get_index(mol, atom);
            obj.setValue(ElementValue{}, index, get_element(mol, atom));
            obj.setValue(IsotopeValue{}, index, get_isotope(mol, atom));
            obj.setValue(ChargeValue{}, index, get_charge(mol, atom));
            obj.setValue(DegreeValue{}, index, get_degree(mol, atom));
            obj.setValue(ImplicitHydrogensValue{}, index, get_implicit_hydrogens(mol, atom));
            obj.setValue(AromaticAtomValue{}, index, is_aromatic_atom(mol, atom));
        }

        for (auto bond : get_bonds(mol)) {
            auto index = get_index(mol, bond);
            obj.setValue(SourceValue{}, index, get_index(mol, get_source(mol, bond)));
            obj.setValue(TargetValue{}, index, get_index(mol, get_target(mol, bond)));
            obj.setValue(OrderValue{}, index, get_order(mol, bond));
            obj.setValue(CyclicBondValue{}, index, is_cyclic_bond(mol, bond));
            obj.setValue(AromaticBondValue{}, index, is_aromatic_bond(mol, bond));
        }

        if constexpr (ctll::exists_in(IncidentArray{}, Layout::type)) {
            static_assert(ctll::exists_in(IncidentBeginArray{}, Layout::type));
            SizeT index = 0;
            for (auto atom : get_atoms(mol)) {
                obj.setValue(IncidentBegin{}, get_index(mol, atom), index);
                for (auto bond : get_bonds(mol, atom))
                    obj.setValue(IncidentIndex{}, index++, get_index(mol, bond));
            }
        }

        if constexpr (ctll::exists_in(AdjacentArray{}, Layout::type)) {
            static_assert(ctll::exists_in(AdjacentBeginArray{}, Layout::type));
            SizeT index = 0;
            for (auto atom : get_atoms(mol)) {
                obj.setValue(AdjacentBegin{}, get_index(mol, atom), index);
                for (auto nbr : get_nbrs(mol, atom))
                    obj.setValue(AdjacentIndex{}, index++, get_index(mol, nbr));
            }
        }
    }

    /**
     * @brief Serialize molecule to buffer.
     */
    template<typename Layout>
    void serialize(Molecule auto &mol, std::vector<std::byte> &data)
    {
        data.resize(moleculeSize<Layout>(num_atoms(mol), num_bonds(mol)));
        serialize<Layout>(mol, data.data());
    }

    /**
     * @brief Serialize molecule to output file stream using existing buffer.
     */
    template<typename Layout>
    void serialize(Molecule auto &mol, std::ofstream &ofs, std::vector<std::byte> &data)
    {
        serialize<Layout>(mol, data);
        ofs.write(reinterpret_cast<char*>(data.data()), data.size());
    }

    /**
     * @brief Serialize molecule to output file stream.
     */
    template<typename Layout>
    void serialize(Molecule auto &mol, std::ofstream &ofs)
    {
        std::vector<std::byte> data;
        serialize<Layout>(mol, data);
        ofs.write(reinterpret_cast<char*>(data.data()), data.size());
    }

    template<typename Layout>
    auto deserialize(std::ifstream &ifs, std::vector<std::byte> &data)
    {
        if (!ifs)
            return std::make_tuple(false, Object<Layout>{});
        data.resize(LayoutSize::size());
        ifs.read(reinterpret_cast<char*>(data.data()), data.size());
        if (!ifs)
            return std::make_tuple(false, Object<Layout>{});
        auto size = *reinterpret_cast<LayoutSize::Type*>(data.data());        
        data.resize(size);
        ifs.read(reinterpret_cast<char*>(data.data() + LayoutSize::size()), size - LayoutSize::size());        
        if (!ifs)
            return std::make_tuple(false, Object<Layout>{});
        return std::make_tuple(true, Object<Layout>{data.data()});
    }

    using MemMapSource = mio::basic_mmap_source<std::byte>;

    template<typename Layout>
    auto deserialize(const MemMapSource &source, std::size_t &offset)
    {
        using Ptr = decltype(source.begin());
        auto data = source.begin() + offset;
        auto size = *reinterpret_cast<const LayoutSize::Type*>(data);
        offset += size;
        return Object<Layout, Ptr>{data};        
    }


} // namespace Kitimar
