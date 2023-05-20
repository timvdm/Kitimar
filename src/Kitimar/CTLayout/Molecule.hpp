#pragma once

#include <Kitimar/Molecule/Molecule.hpp>
#include <Kitimar/CTLayout/Serialize.hpp>

namespace Kitimar::CTLayout {

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

    // Type-Oriented design
    struct NumAtomTypes : ArraySize {};
    struct NumBondTypes : ArraySize {};

    struct AtomTypeMap : Array<AtomStruct, NumAtomTypes> {};
    struct AtomTypeIndex : Value<uint16_t> {};
    struct AtomTypeArray : Array<AtomTypeIndex, NumAtoms> {};

    struct BondTypeMap : Array<BondStruct, NumBondTypes> {};
    struct BondTypeIndex : Value<uint8_t> {};
    struct BondTypeArray : Array<BondTypeIndex, NumBonds> {};

    struct AtomType
    {
        uint8_t element;
        uint8_t isotope;
        int8_t charge;
        uint8_t degree;
        uint8_t implicitHydrogens;
        bool aromatic;

        friend auto operator<=>(const AtomType&, const AtomType&) = default;
    };

    struct BondType
    {
        uint8_t order;
        bool cyclic;
        bool aromatic;

        friend auto operator<=>(const BondType&, const BondType&) = default;
    };



    // incident list (i.e. bonds for atom)
    struct IncidentBegin : Value<SizeT> {};
    struct IncidentBeginArray : Array<IncidentBegin, NumAtoms> {};
    struct IncidentIndex : Value<SizeT> {};
    struct IncidentArray : Array<IncidentIndex, NumBonds, 2> {};
    // FIXME IncidentStruct

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


    struct AtomTypeLayout : Layout<AtomTypeMap> {};
    struct AtomTypeMolecule : Layout<AtomTypeArray, BondArray, IncidentBeginArray, IncidentArray> {};


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
    void serializeIncident(Molecule auto &mol, std::byte *data)
    {
        if constexpr (ctll::exists_in(IncidentArray{}, Layout::type)) {
            static_assert(ctll::exists_in(IncidentBeginArray{}, Layout::type));
            auto obj = EditableObject<Layout>{data};
            SizeT index = 0;
            for (auto atom : get_atoms(mol)) {
                obj.setValue(IncidentBegin{}, get_index(mol, atom), index);
                for (auto bond : get_bonds(mol, atom))
                    obj.setValue(IncidentIndex{}, index++, get_index(mol, bond));
            }
        }
    }

    template<typename Layout>
    void serializeAdjacent(Molecule auto &mol, std::byte *data)
    {
        if constexpr (ctll::exists_in(AdjacentArray{}, Layout::type)) {
            static_assert(ctll::exists_in(AdjacentBeginArray{}, Layout::type));
            auto obj = EditableObject<Layout>{data};
            SizeT index = 0;
            for (auto atom : get_atoms(mol)) {
                obj.setValue(AdjacentBegin{}, get_index(mol, atom), index);
                for (auto nbr : get_nbrs(mol, atom))
                    obj.setValue(AdjacentIndex{}, index++, get_index(mol, nbr));
            }
        }
    }

    template<typename Layout>
    void serializeAtomArray(Molecule auto &mol, std::byte *data)
    {
        if constexpr (ctll::exists_in(AtomArray{}, Layout::type)) {
            auto obj = EditableObject<Layout>{data};
            for (auto atom : get_atoms(mol)) {
                auto index = get_index(mol, atom);
                obj.setValue(ElementValue{}, index, get_element(mol, atom));
                obj.setValue(IsotopeValue{}, index, get_isotope(mol, atom));
                obj.setValue(ChargeValue{}, index, get_charge(mol, atom));
                obj.setValue(DegreeValue{}, index, get_degree(mol, atom));
                obj.setValue(ImplicitHydrogensValue{}, index, get_implicit_hydrogens(mol, atom));
                obj.setValue(AromaticAtomValue{}, index, is_aromatic_atom(mol, atom));
            }
        }
    }

    template<typename Layout>
    void serializeBondArray(Molecule auto &mol, std::byte *data)
    {
        if constexpr (ctll::exists_in(BondArray{}, Layout::type)) {
            auto obj = EditableObject<Layout>{data};
            for (auto bond : get_bonds(mol)) {
                auto index = get_index(mol, bond);
                obj.setValue(SourceValue{}, index, get_index(mol, get_source(mol, bond)));
                obj.setValue(TargetValue{}, index, get_index(mol, get_target(mol, bond)));
                obj.setValue(OrderValue{}, index, get_order(mol, bond));
                obj.setValue(CyclicBondValue{}, index, is_cyclic_bond(mol, bond));
                obj.setValue(AromaticBondValue{}, index, is_aromatic_bond(mol, bond));
            }
        }
    }

    template<typename Layout>
    auto findTypes(auto &source)
    {
        std::vector<AtomType> atomTypes;
        std::vector<BondType> bondTypes;

        if constexpr (ctll::exists_in(AtomTypeArray{}, Layout::type)) { // FIXME or BondTypeArray
            for (auto mol : source.molecules()) {
                for (auto atom : get_atoms(mol)) {
                    auto type = AtomType{
                        static_cast<uint8_t>(get_element(mol, atom)),
                        static_cast<uint8_t>(get_isotope(mol, atom)),
                        static_cast<int8_t>(get_charge(mol, atom)),
                        static_cast<uint8_t>(get_degree(mol, atom)),
                        static_cast<uint8_t>(get_implicit_hydrogens(mol, atom)),
                        is_aromatic_atom(mol, atom)
                    };

                    if (std::ranges::find(atomTypes, type) == atomTypes.end())
                        atomTypes.push_back(type);
                }

                for (auto bond : get_bonds(mol)) {
                    auto type = BondType{
                        static_cast<uint8_t>(get_order(mol, bond)),
                        is_cyclic_bond(mol, bond),
                        is_aromatic_bond(mol, bond)
                    };

                    if (std::ranges::find(bondTypes, type) == bondTypes.end())
                        bondTypes.push_back(type);
                }

            }

            source.reset();
        }

        return std::make_tuple(atomTypes, bondTypes);
    }

    void serializeAtomTypeMap(FileStreamSink &sink, std::vector<std::byte> &buffer,
                            const std::vector<AtomType> &types)
    {
        if (types.empty())
            return;

        buffer.resize(AtomTypeLayout::size({static_cast<SizeT>(types.size())}));

        auto obj = EditableObject<AtomTypeLayout>{buffer.data()};
        obj.setSize(NumAtomTypes{}, types.size());
        obj.setSize(obj.size());

        for (auto i = 0; i < types.size(); ++i) {
            const auto &type = types[i];
            obj.setValue(ElementValue{}, i, type.element);
            obj.setValue(IsotopeValue{}, i, type.isotope);
            obj.setValue(ChargeValue{}, i, type.charge);
            obj.setValue(DegreeValue{}, i, type.degree);
            obj.setValue(ImplicitHydrogensValue{}, i, type.implicitHydrogens);
            obj.setValue(AromaticAtomValue{}, i, type.aromatic);
        }

        sink.write(buffer);
    }

    template<typename Layout>
    void serializeAtomTypeArray(Molecule auto &mol, std::byte *data,
                                const std::vector<AtomType> &atomTypes)
    {
        if constexpr (ctll::exists_in(AtomTypeArray{}, Layout::type)) {
            auto obj = EditableObject<Layout>{data};

            for (auto atom : get_atoms(mol)) {
                auto type = AtomType{
                    static_cast<uint8_t>(get_element(mol, atom)),
                    static_cast<uint8_t>(get_isotope(mol, atom)),
                    static_cast<int8_t>(get_charge(mol, atom)),
                    static_cast<uint8_t>(get_degree(mol, atom)),
                    static_cast<uint8_t>(get_implicit_hydrogens(mol, atom)),
                    is_aromatic_atom(mol, atom)
                };

                auto typeIndex = std::ranges::find(atomTypes, type) - atomTypes.begin();
                obj.setValue(AtomTypeIndex{}, get_index(mol, atom), typeIndex);
            }
        }
    }




    template<typename Layout>
    void serializeMolecule(Molecule auto &mol, std::byte *data,
                           const std::vector<AtomType> &atomTypes)
    {            
        auto obj = EditableObject<Layout>{data};
        obj.setSize(NumAtoms{}, num_atoms(mol));
        obj.setSize(NumBonds{}, num_bonds(mol));
        obj.setSize(obj.size());

        serializeAtomTypeArray<Layout>(mol, data, atomTypes);

        serializeAtomArray<Layout>(mol, data);
        serializeBondArray<Layout>(mol, data);

        serializeIncident<Layout>(mol, data);
        serializeAdjacent<Layout>(mol, data);

    }






    template<typename Layout>
    constexpr auto moleculeSize(SizeT numAtoms, SizeT numBonds) noexcept
    {
        return Layout::size({numAtoms, numBonds});
    }


    template<typename Layout>
    void serializeMolSource(auto &source, std::string_view filename)
    {
        FileStreamSink sink{filename};
        std::vector<std::byte> buffer;


        auto [atomTypes, bondTypes] = findTypes<Layout>(source);

        std::cout << "# atom types: " << atomTypes.size() << std::endl;
        std::cout << "# bond types: " << bondTypes.size() << std::endl;




        for (auto mol : source.molecules()) {
            buffer.resize(moleculeSize<Layout>(num_atoms(mol), num_bonds(mol)));
            serializeMolecule<Layout>(mol, buffer.data(), atomTypes);
            sink.writeObject(buffer);
            if (sink.size() % 10000 == 0) std::cout << sink.size() << std::endl;
        }
    }



} // namespace Kitimar::CTLayout
