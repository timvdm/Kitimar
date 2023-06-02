#pragma once

#include <Kitimar/Util/Util.hpp>
#include <Kitimar/Molecule/Molecule.hpp>
#include <Kitimar/CTLayout/Value.hpp>
#include <Kitimar/CTLayout/Struct.hpp>
#include <Kitimar/CTLayout/Vector.hpp>
#include <Kitimar/CTLayout/Sink.hpp>

namespace Kitimar::CTLayout {

    /*
    using NumAtoms = uint32_t;
    using NumBonds = uint32_t;
    using NumNbrs = uint32_t;
    */

    using NumAtoms = uint16_t;
    using NumBonds = uint16_t;
    using NumNbrs = uint8_t;

    // Atom values

    struct ElementValue : Value<uint8_t> {};
    struct IsotopeValue : Value<uint8_t> {};
    struct ChargeValue : Value<int8_t> {};
    struct DegreeValue : Value<uint8_t> {};
    struct ImplicitHydrogensValue : Value<uint8_t> {};
    struct AromaticAtomValue : Value<bool> {};

    // Bond values

    struct SourceValue : Value<NumAtoms> {};
    struct TargetValue : Value<NumAtoms> {};
    struct OrderValue : Value<uint8_t> {};
    struct CyclicBondValue : Value<bool> {};
    struct AromaticBondValue : Value<bool> {};

    // Incident list (i.e. bonds for atom)

    struct IncidentIndex : Value<NumAtoms> {};
    struct IncidentList : Vector<IncidentIndex, NumNbrs> {};
    struct IncidentListList : Vector<IncidentList, NumAtoms> {};

    // Adjacent list (i.e. nbrs of atom)

    struct AdjacentIndex : Value<NumAtoms> {};
    struct AdjacentList : Vector<AdjacentIndex, NumNbrs> {};
    struct AdjacentListList : Vector<AdjacentList, NumAtoms> {};

    // Object-Oriented design using Struct

    struct AtomStruct : Struct<ElementValue, IsotopeValue, ChargeValue, DegreeValue, ImplicitHydrogensValue, AromaticAtomValue> {};
    struct AtomList : Vector<AtomStruct, NumAtoms> {};

    struct BondStruct : Struct<SourceValue, TargetValue, OrderValue, AromaticBondValue, CyclicBondValue> {};
    struct BondList : Vector<BondStruct, NumBonds> {};

    struct StructMolecule : Struct<AtomList, BondList> {};
    struct StructMoleculeIncident : Struct<AtomList, BondList, IncidentListList> {};
    struct StructMoleculeAdjacent : Struct<AtomList, BondList, AdjacentListList> {};
    struct StructMoleculeIncidentAdjacent : Struct<AtomList, BondList,
                                                   IncidentListList,
                                                   AdjacentListList> {};

    struct StructMolecules : Vector<StructMoleculeIncident> {};

    // Data-Oriented design using no Struct

    struct ElementList : Vector<ElementValue, NumAtoms> {};
    struct IsotopeList : Vector<IsotopeValue, NumAtoms> {};
    struct ChargeList : Vector<ChargeValue, NumAtoms> {};
    struct DegreeList : Vector<DegreeValue, NumAtoms> {};
    struct ImplicitHydrogensList : Vector<ImplicitHydrogensValue, NumAtoms> {};
    struct AromaticAtomList : Vector<AromaticAtomValue, NumAtoms> {};
    //struct AromaticAtomList : BitVector<NumAtoms> {};


    struct SourceList : Vector<SourceValue, NumBonds> {};
    struct TargetList : Vector<TargetValue, NumBonds> {};
    struct OrderList : Vector<OrderValue, NumBonds> {};
    struct CyclicBondList : Vector<CyclicBondValue, NumBonds> {};
    struct AromaticBondList : Vector<AromaticBondValue, NumBonds> {};
    //struct AromaticBondList : BitVector<NumBonds> {};


    struct ListMolecule : Struct<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList, AromaticAtomList,
                                 SourceList, TargetList, OrderList, CyclicBondList, AromaticBondList> {};
    struct ListMoleculeIncident : Struct<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList, AromaticAtomList,
                                         SourceList, TargetList, OrderList, CyclicBondList, AromaticBondList, IncidentListList> {};
    struct ListMoleculeAdjacent : Struct<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList, AromaticAtomList,
                                         SourceList, TargetList, OrderList, CyclicBondList, AromaticBondList, AdjacentListList> {};
    struct ListMoleculeIncidentAdjacent : Struct<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList, AromaticAtomList,
                                                 SourceList, TargetList, OrderList, CyclicBondList, AromaticBondList,
                                                 IncidentListList, AdjacentListList> {};

    struct ListMolecules : Vector<ListMoleculeIncident> {};

    // Type-Oriented design    

    struct AtomTypeList : Vector<AtomStruct> {};
    struct AtomTypeIndex : Value<uint16_t> {};
    struct AtomTypeIndexList : Vector<AtomTypeIndex, NumAtoms> {};

    struct BondTypeStruct : Struct<OrderValue, AromaticBondValue, CyclicBondValue> {};
    struct BondTypeList : Vector<BondTypeStruct> {};
    struct BondTypeIndex : Value<uint8_t> {};
    struct BondTypeIndexList : Vector<BondTypeIndex, NumBonds> {};

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



    struct TypeMolecule : Struct<AtomTypeIndexList, BondTypeIndexList, SourceList, TargetList, IncidentListList> {};
    struct TypeMolecules : Struct<AtomTypeList, BondTypeList, Vector<TypeMolecule>> {};


    //
    // Molecule
    //

    template<typename MolObj>
    constexpr auto num_atoms(const MolObj &mol) noexcept
    {
        if constexpr (requires { mol.numAtoms(); })
            return mol.numAtoms();
        else if constexpr (contains(typename MolObj::Type{}, AtomList{}))
            return mol.get(AtomList{}).length();
        else
            return mol.get(ElementList{}).length();

    }

    template<typename MolObj>
    constexpr auto num_bonds(const MolObj &mol) noexcept
    {
        if constexpr (requires { mol.numAtoms(); })
            return mol.numBonds();
        else if constexpr (contains(typename MolObj::Type{}, BondList{}))
            return mol.get(BondList{}).length();
        else
            return mol.get(SourceList{}).length();
    }

    constexpr auto get_atoms(const auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_atoms(mol))>(0), num_atoms(mol));
    }

    constexpr auto get_bonds(const auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_bonds(mol))>(0), num_bonds(mol));
    }

    constexpr auto get_atom(const auto &mol, auto index) noexcept
    {
        assert(index < num_atoms(mol));
        return index;
    }

    constexpr auto get_bond(const auto &mol, auto index) noexcept
    {
        assert(index < num_bonds(mol));
        return index;
    }

    constexpr auto get_index(const auto &mol, auto index) noexcept
    {
        return index;
    }

    //
    // Atom
    //

    template<typename MolObj, typename Property, typename PropertyList>
    constexpr auto get_atom_property(const MolObj &mol, auto atom, Property, PropertyList) noexcept
    {
        if constexpr (contains(typename MolObj::Type{}, AtomList{}))
            return mol.get(AtomList{}).at(atom).get(Property{}).get();
        else
            return mol.get(PropertyList{}).at(atom).get();
    }

    template<typename MolObj>
    constexpr auto get_element(const MolObj &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (requires { mol.degree(0); })
            return mol.element(atom);
        else
            return get_atom_property(mol, atom, ElementValue{}, ElementList{});
    }

    constexpr auto get_isotope(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (requires { mol.isotope(0); })
            return mol.isotope(atom);
        else
            return get_atom_property(mol, atom, IsotopeValue{}, IsotopeList{});
    }

    constexpr auto get_charge(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (requires { mol.charge(0); })
            return mol.charge(atom);
        else
            return get_atom_property(mol, atom, ChargeValue{}, ChargeList{});
    }

    constexpr auto get_degree(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        if constexpr (requires { mol.degree(0); })
            return mol.degree(atom);
        else
            return get_atom_property(mol, atom, DegreeValue{}, DegreeList{});
    }

    constexpr auto get_implicit_hydrogens(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        if constexpr (requires { mol.implicitHydrogens(0); })
            return mol.implicitHydrogens(atom);
        else
            return get_atom_property(mol, atom, ImplicitHydrogensValue{}, ImplicitHydrogensList{});
    }

    constexpr auto is_aromatic_atom(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (requires { mol.isAromaticAtom(0); })
            return mol.isAromaticAtom(atom);
        else
            return get_atom_property(mol, atom, AromaticAtomValue{}, AromaticAtomList{});
    }

    template<typename MolObj>
    constexpr auto get_bonds(const MolObj &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (requires { mol.incident(0); })
            return mol.incident(atom);
        else if constexpr (contains(typename MolObj::Type{}, IncidentListList{})) {
            auto bonds = mol.get(IncidentListList{}).at(atom);
            return bonds.range();
        } else {
            return get_bonds(mol) | std::views::filter([&mol, atom] (auto i) {
                auto bond = get_bond(mol, i);
                return atom == get_source(mol, bond) || atom == get_target(mol, bond);
            });
        }
    }

    template<typename MolObj>
    constexpr auto get_nbrs(const MolObj &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
            return Kitimar::Molecule::get_nbr(mol, bond, atom);
        });
        /*
        if constexpr (contains(typename MolObj::Type{}, AdjacentListList{})) {
            auto nbrs = mol.get(AdjacentListList{}).at(atom);
            return nbrs.range();
        } else {
            return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
                return get_nbr(mol, bond, atom);
            });
        }
        */
    }

    template<typename MolObj>
    constexpr auto null_atom(const MolObj &mol) noexcept
    {
        return -1;
    }

    //
    // Bond
    //

    template<typename MolObj, typename Property, typename PropertyList>
    constexpr auto get_bond_property(const MolObj &mol, auto bond, Property, PropertyList) noexcept
    {        
        if constexpr (contains(typename MolObj::Type{}, BondList{}))
            return mol.get(BondList{}).at(bond).get(Property{}).get();
        else
            return mol.get(PropertyList{}).at(bond).get();
    }

    constexpr auto get_source(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        if constexpr (requires { mol.source(0); })
            return mol.source(bond);
        else
            return get_bond_property(mol, bond, SourceValue{}, SourceList{});
    }

    constexpr auto get_target(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        if constexpr (requires { mol.target(0); })
            return mol.target(bond);
        else
            return get_bond_property(mol, bond, TargetValue{}, TargetList{});
    }

    constexpr auto get_order(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        if constexpr (requires { mol.order(0); })
            return mol.order(bond);
        else
            return get_bond_property(mol, bond, OrderValue{}, OrderList{});
    }

    constexpr auto is_cyclic_bond(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        if constexpr (requires { mol.isCyclicBond(0); })
            return mol.isCyclicBond(bond);
        else
            return get_bond_property(mol, bond, CyclicBondValue{}, CyclicBondList{});
    }

    constexpr auto is_aromatic_bond(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        if constexpr (requires { mol.isAromaticBond(0); })
            return mol.isAromaticBond(bond);
        else
            return get_bond_property(mol, bond, AromaticBondValue{}, AromaticBondList{});
    }

    /*
    constexpr auto get_nbr(const auto &mol, auto bond, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        assert(bond < num_bonds(mol));
        auto source = get_source(mol, bond);
        return source != atom ? source : get_target(mol, bond);
    }
    */

    template<typename MolObj>
    constexpr auto null_bond(const MolObj &mol) noexcept
    {
        return -1;
    }

    namespace detail {

        template<typename Layout>
        void writeIncident(Molecule::Molecule auto &mol, auto &writer)
        {
            if constexpr (contains(Layout{}, IncidentListList{})) {                
                auto incidentList = writer.get(IncidentListList{});
                incidentList.setLength(num_atoms(mol));
                for (auto atom : get_atoms(mol)) {                    
                    auto incident = incidentList.at(get_index(mol, atom));                    
                    incident.setLength(get_degree(mol, atom));
                    int i = 0;
                    for (auto bond : get_bonds(mol, atom)) {
                        assert(i < get_degree(mol, atom));
                        incident.at(i++).write(static_cast<IncidentIndex::Type>(get_index(mol, bond)));
                    }
                }
            }
        }

        template<typename Layout>
        void writeAdjacent(Molecule::Molecule auto &mol, auto &writer)
        {
            if constexpr (contains(Layout{}, AdjacentListList{})) {
                auto adjacentList = writer.get(AdjacentListList{});
                adjacentList.setLength(num_atoms(mol));
                for (auto atom : get_atoms(mol)) {
                    auto adjacent = adjacentList.at(get_index(mol, atom));
                    adjacent.setLength(get_degree(mol, atom));
                    int i = 0;
                    for (auto nbr : get_nbrs(mol, atom))
                        adjacent.at(i++).write(static_cast<IncidentIndex::Type>(get_index(mol, nbr)));
                }
            }
        }

        template<typename Layout>
        void writeAtomList(Molecule::Molecule auto &mol, auto &writer)
        {
            if constexpr (contains(Layout{}, AtomList{})) {
                auto atomList = writer.get(AtomList{});
                atomList.setLength(num_atoms(mol));
                for (auto atom : get_atoms(mol)) {
                    auto index = get_index(mol, atom);
                    auto atomWriter = atomList.at(index);
                    atomWriter.get(ElementValue{}).write(static_cast<ElementValue::Type>(get_element(mol, atom)));
                    atomWriter.get(IsotopeValue{}).write(static_cast<IsotopeValue::Type>(get_isotope(mol, atom)));
                    atomWriter.get(ChargeValue{}).write(static_cast<ChargeValue::Type>(get_charge(mol, atom)));
                    atomWriter.get(DegreeValue{}).write(static_cast<DegreeValue::Type>(get_degree(mol, atom)));
                    atomWriter.get(ImplicitHydrogensValue{}).write(static_cast<ImplicitHydrogensValue::Type>(get_implicit_hydrogens(mol, atom)));
                    atomWriter.get(AromaticAtomValue{}).write(static_cast<AromaticAtomValue::Type>(is_aromatic_atom(mol, atom)));
                }
            } else if constexpr (contains(Layout{}, ElementList{})) {
                auto elementList = writer.get(ElementList{});
                elementList.setLength(num_atoms(mol));
                auto isotopeList = writer.get(IsotopeList{});
                isotopeList.setLength(num_atoms(mol));
                auto chargeList = writer.get(ChargeList{});
                chargeList.setLength(num_atoms(mol));
                auto degreeList = writer.get(DegreeList{});
                degreeList.setLength(num_atoms(mol));
                auto implicitHydrogensList = writer.get(ImplicitHydrogensList{});
                implicitHydrogensList.setLength(num_atoms(mol));
                auto aromaticList = writer.get(AromaticAtomList{});
                aromaticList.setLength(num_atoms(mol));
                for (auto atom : get_atoms(mol)) {
                    auto index = get_index(mol, atom);
                    auto elementWriter = elementList.at(index);
                    elementWriter.write(static_cast<ElementValue::Type>(get_element(mol, atom)));
                    auto isotopeWriter = isotopeList.at(index);
                    isotopeWriter.write(static_cast<IsotopeValue::Type>(get_isotope(mol, atom)));
                    auto chargeWriter = chargeList.at(index);
                    chargeWriter.write(static_cast<ChargeValue::Type>(get_charge(mol, atom)));
                    auto degreeWriter = degreeList.at(index);
                    degreeWriter.write(static_cast<DegreeValue::Type>(get_degree(mol, atom)));
                    auto implicitHydrogensWriter = implicitHydrogensList.at(index);
                    implicitHydrogensWriter.write(static_cast<ImplicitHydrogensValue::Type>(get_implicit_hydrogens(mol, atom)));
                    auto aromaticAtomWriter = aromaticList.at(index);
                    aromaticAtomWriter.write(static_cast<AromaticAtomValue::Type>(is_aromatic_atom(mol, atom)));
                }
            }
        }

        template<typename Layout>
        void writeBondList(Molecule::Molecule auto &mol, auto &writer)
        {
            if constexpr (std::is_same_v<Layout, TypeMolecule>) {
                auto sourceList = writer.get(SourceList{});
                sourceList.setLength(num_bonds(mol));
                auto targetList = writer.get(TargetList{});
                targetList.setLength(num_bonds(mol));
                for (auto bond : get_bonds(mol)) {
                    auto index = get_index(mol, bond);
                    auto sourceWriter = sourceList.at(index);
                    sourceWriter.write(static_cast<SourceValue::Type>(get_index(mol, get_source(mol, bond))));
                    auto targetWriter = targetList.at(index);
                    targetWriter.write(static_cast<TargetValue::Type>(get_index(mol, get_target(mol, bond))));
                }
            } else if constexpr (contains(Layout{}, BondList{})) {
                auto bondList = writer.get(BondList{});
                bondList.setLength(num_bonds(mol));
                for (auto bond : get_bonds(mol)) {
                    auto index = get_index(mol, bond);
                    auto bondWriter = bondList.at(index);
                    bondWriter.get(SourceValue{}).write(static_cast<SourceValue::Type>(get_index(mol, get_source(mol, bond))));
                    bondWriter.get(TargetValue{}).write(static_cast<TargetValue::Type>(get_index(mol, get_target(mol, bond))));
                    bondWriter.get(OrderValue{}).write(static_cast<OrderValue::Type>(get_order(mol, bond)));
                    bondWriter.get(CyclicBondValue{}).write(static_cast<CyclicBondValue::Type>(is_cyclic_bond(mol, bond)));
                    bondWriter.get(AromaticBondValue{}).write(static_cast<AromaticBondValue::Type>(is_aromatic_bond(mol, bond)));
                }
            } else if constexpr (contains(Layout{}, SourceList{})) {
                auto sourceList = writer.get(SourceList{});
                sourceList.setLength(num_bonds(mol));
                auto targetList = writer.get(TargetList{});
                targetList.setLength(num_bonds(mol));
                auto orderList = writer.get(OrderList{});
                orderList.setLength(num_bonds(mol));
                auto cyclicList = writer.get(CyclicBondList{});
                cyclicList.setLength(num_bonds(mol));
                auto aromaticList = writer.get(AromaticBondList{});
                aromaticList.setLength(num_bonds(mol));
                for (auto bond : get_bonds(mol)) {
                    auto index = get_index(mol, bond);
                    auto sourceWriter = sourceList.at(index);
                    sourceWriter.write(static_cast<SourceValue::Type>(get_index(mol, get_source(mol, bond))));
                    auto targetWriter = targetList.at(index);
                    targetWriter.write(static_cast<TargetValue::Type>(get_index(mol, get_target(mol, bond))));
                    auto orderWriter = orderList.at(index);
                    orderWriter.write(static_cast<OrderValue::Type>(get_order(mol, bond)));
                    auto cyclicWriter = cyclicList.at(index);
                    cyclicWriter.write(static_cast<CyclicBondValue::Type>(is_cyclic_bond(mol, bond)));
                    auto aromaticWriter = aromaticList.at(index);
                    aromaticWriter.write(static_cast<AromaticBondValue::Type>(is_aromatic_bond(mol, bond)));
                }
            }            
        }

        template<typename Layout>
        auto findTypes(auto &source)
        {
            std::size_t numMolecules = 0;
            std::vector<AtomType> atomTypes;
            std::vector<BondType> bondTypes;

            if constexpr (contains(Layout{}, AtomTypeList{})) { // FIXME or BondTypeList
                for (auto mol : source.molecules()) {
                    ++numMolecules;

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

            return std::make_tuple(numMolecules, atomTypes, bondTypes);
        }

        void writeAtomTypeList(auto &writer, const std::vector<AtomType> &atomTypes)
        {
            auto atomTypeList = writer.get(AtomTypeList{});
            atomTypeList.setLength(atomTypes.size());
            for (auto i = 0; i < atomTypes.size(); ++i) {
                auto &atomType = atomTypes[i];
                auto atomTypeWriter = atomTypeList.at(i);
                atomTypeWriter.get(ElementValue{}).write(static_cast<ElementValue::Type>(atomType.element)); // FIXME : cast not needed anymore...
                atomTypeWriter.get(IsotopeValue{}).write(static_cast<IsotopeValue::Type>(atomType.isotope));
                atomTypeWriter.get(ChargeValue{}).write(static_cast<ChargeValue::Type>(atomType.charge));
                atomTypeWriter.get(DegreeValue{}).write(static_cast<DegreeValue::Type>(atomType.degree));
                atomTypeWriter.get(ImplicitHydrogensValue{}).write(static_cast<ImplicitHydrogensValue::Type>(atomType.implicitHydrogens));
                atomTypeWriter.get(AromaticAtomValue{}).write(static_cast<AromaticAtomValue::Type>(atomType.aromatic));
            }
        }

        void writeBondTypeList(auto &writer, const std::vector<BondType> &bondTypes)
        {
            auto bondTypeList = writer.get(BondTypeList{});
            bondTypeList.setLength(bondTypes.size());
            for (auto i = 0; i < bondTypes.size(); ++i) {
                auto &bondType = bondTypes[i];
                auto bondTypeWriter = bondTypeList.at(i);
                bondTypeWriter.get(OrderValue{}).write(static_cast<OrderValue::Type>(bondType.order));
                bondTypeWriter.get(CyclicBondValue{}).write(static_cast<CyclicBondValue::Type>(bondType.cyclic));
                bondTypeWriter.get(AromaticBondValue{}).write(static_cast<AromaticBondValue::Type>(bondType.aromatic));
            }
        }

        void writeAtomTypeIndexList(Molecule::Molecule auto &mol, auto &writer, const std::vector<AtomType> &atomTypes)
        {
            auto atomTypeIndexList = writer.get(AtomTypeIndexList{});
            atomTypeIndexList.setLength(num_atoms(mol));

            for (auto atom : get_atoms(mol)) {
                auto index = get_index(mol, atom);
                auto atomType = AtomType{
                    static_cast<uint8_t>(get_element(mol, atom)),
                    static_cast<uint8_t>(get_isotope(mol, atom)),
                    static_cast<int8_t>(get_charge(mol, atom)),
                    static_cast<uint8_t>(get_degree(mol, atom)),
                    static_cast<uint8_t>(get_implicit_hydrogens(mol, atom)),
                    is_aromatic_atom(mol, atom)
                };
                auto atomTypeIndex = Util::indexOf(atomTypes, atomType);
                atomTypeIndexList.at(index).write(atomTypeIndex);
            }
        }

        void writeBondTypeIndexList(Molecule::Molecule auto &mol, auto &writer, const std::vector<BondType> &bondTypes)
        {
            auto bondTypeIndexList = writer.get(BondTypeIndexList{});
            bondTypeIndexList.setLength(num_bonds(mol));
            for (auto bond : get_bonds(mol)) {
                auto index = get_index(mol, bond);
                auto bondType = BondType{
                    static_cast<uint8_t>(get_order(mol, bond)),
                    static_cast<bool>(is_cyclic_bond(mol, bond)),
                    static_cast<bool>(is_aromatic_bond(mol, bond)),

                };
                auto bondTypeIndex = Util::indexOf(bondTypes, bondType);
                bondTypeIndexList.at(index).write(bondTypeIndex);
            }
        }

        /*


        template<typename Layout>
        void serializeAtomTypeList(Molecule::Molecule auto &mol, std::byte *data,
                                    const std::vector<AtomType> &atomTypes)
        {
            if constexpr (ctll::exists_in(AtomTypeList{}, Layout::type)) {
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

        */

        template<typename Layout>
        void writeMolecule(Molecule::Molecule auto &mol, auto &writer,
                           const std::vector<AtomType> &atomTypes = {},
                           const std::vector<BondType> &bondTypes = {})
        {
            //fmt::println("WriteMolecule()");
            //serializeAtomTypeList<Layout>(mol, data, atomTypes);

            if constexpr (std::is_same_v<Layout, TypeMolecule>) {
                writeAtomTypeIndexList(mol, writer, atomTypes);
                writeBondTypeIndexList(mol, writer, bondTypes);
            }

            writeAtomList<Layout>(mol, writer);
            writeBondList<Layout>(mol, writer);

            writeIncident<Layout>(mol, writer);
            writeAdjacent<Layout>(mol, writer);

        }

    } // namespace detail

    template<typename Layout, typename Sink>
    void serializeMolecule(Molecule::Molecule auto &mol, Sink &sink,
                           const std::vector<AtomType> &atomTypes = {})
    {
        auto writer = toWriter(Layout{}, sink);
        detail::writeMolecule<Layout>(mol, writer, atomTypes);
    }

    template<typename Layout>
    void serializeMolSource(auto &source, auto &sink) requires (isVector(Layout{}))
    {
        auto writer = toWriter(Layout{}, sink);
        auto numMolecules =  source.numMolecules();
        writer.setLength(numMolecules);
        for (auto i = 0; i < numMolecules; ++i) {
            if (i && i % 1000 == 0)
                std::cout << i << std::endl;
            auto molWriter = writer.at(i);
            auto mol = source.read();
            detail::writeMolecule<typename Layout::Type>(mol, molWriter);
        }
    }

    template<typename Layout>
    void serializeMolSource(auto &source, auto &sink) requires (std::same_as<Layout, TypeMolecules>)
    {
        auto writer = toWriter(Layout{}, sink);



        auto [numMolecules, atomTypes, bondTypes] = detail::findTypes<Layout>(source);

        /*
        std::cout << "# molecules: " << numMolecules << std::endl;
        std::cout << "# atom types: " << atomTypes.size() << std::endl;
        std::cout << "# bond types: " << bondTypes.size() << std::endl;
        */

        detail::writeAtomTypeList(writer, atomTypes);
        detail::writeBondTypeList(writer, bondTypes);

        auto moleculeList = writer.get(Vector<TypeMolecule>{});
        moleculeList.setLength(numMolecules);

        for (auto i = 0; i < numMolecules; ++i) {
            if (i && i % 1000 == 0)
                std::cout << i << std::endl;
            auto molWriter = moleculeList.at(i);
            auto mol = source.read();
            detail::writeMolecule<TypeMolecule>(mol, molWriter, atomTypes, bondTypes);
        }




        //for (auto mol : molSource.molecules()) {
            //buffer.resize(moleculeSize<Layout>(num_atoms(mol), num_bonds(mol)));
            //serializeMolecule<Layout>(mol, buffer.data(), atomTypes);
            //sink.writeObject(buffer);
            //if (sink.size() % 10000 == 0) std::cout << sink.size() << std::endl;
        //}
    }

    template<typename Object>
    class StructMoleculeObject
    {
        public:
            StructMoleculeObject(const Object &object)
                : m_atoms(object.get(AtomList{})), m_bonds(object.get(BondList{}))
            {
                static_assert(contains(typename Object::Type{}, IncidentListList{}));
                m_incident = object.get(IncidentListList{});
            }

            auto numAtoms() const noexcept
            {
                return m_atoms.length();
            }

            auto numBonds() const noexcept
            {
                return m_bonds.length();
            }

            auto element(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(ElementValue{}).get();
            }

            auto isotope(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(IsotopeValue{}).get();
            }

            auto charge(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(ChargeValue{}).get();
            }

            auto degree(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(DegreeValue{}).get();
            }

            auto implicitHydrogens(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(ImplicitHydrogensValue{}).get();
            }

            auto isAromaticAtom(auto atomIdx) const noexcept
            {
                return m_atoms.at(atomIdx).get(AromaticAtomValue{}).get();
            }

            auto incident(auto atomIdx) const noexcept
            {
                return m_incident.at(atomIdx).range();
            }

            auto source(auto bondIdx) const noexcept
            {
                return m_bonds.at(bondIdx).get(SourceValue{}).get();
            }

            auto target(auto bondIdx) const noexcept
            {
                return m_bonds.at(bondIdx).get(TargetValue{}).get();
            }

            auto order(auto bondIdx) const noexcept
            {
                return m_bonds.at(bondIdx).get(OrderValue{}).get();
            }


            auto isCyclicBond(auto bondIdx) const noexcept
            {
                return m_bonds.at(bondIdx).get(CyclicBondValue{}).get();
            }

            auto isAromaticBond(auto bondIdx) const noexcept
            {
                return m_bonds.at(bondIdx).get(AromaticBondValue{}).get();
            }

        private:
            VectorObject<AtomList, typename Object::Source> m_atoms;
            VectorObject<BondList, typename Object::Source> m_bonds;
            VectorObject<IncidentListList, typename Object::Source> m_incident;
    };

    template<typename Object>
    class ListMoleculeObject
    {
        public:
            ListMoleculeObject(const Object &object)
                : m_elements(object.get(ElementList{})), m_isotopes(object.get(IsotopeList{})),
                  m_charges(object.get(ChargeList{})), m_degrees(object.get(DegreeList{})),
                  m_implicitHydrogens(object.get(ImplicitHydrogensList{})), m_aromaticAtoms(object.get(AromaticAtomList{})),
                  m_sources(object.get(SourceList{})), m_targets(object.get(TargetList{})),
                  m_orders(object.get(OrderList{})), m_cyclicBonds(object.get(CyclicBondList{})),
                  m_aromaticBonds(object.get(AromaticBondList{}))

            {
                static_assert(contains(typename Object::Type{}, IncidentListList{}));
                m_incident = object.get(IncidentListList{});
            }

            auto numAtoms() const noexcept
            {
                return m_elements.length();
            }

            auto numBonds() const noexcept
            {
                return m_sources.length();
            }

            auto element(auto atomIdx) const noexcept
            {
                return m_elements.at(atomIdx).get();
            }

            auto isotope(auto atomIdx) const noexcept
            {
                return m_isotopes.at(atomIdx).get();
            }

            auto charge(auto atomIdx) const noexcept
            {
                return m_charges.at(atomIdx).get();
            }

            auto degree(auto atomIdx) const noexcept
            {
                return m_degrees.at(atomIdx).get();
            }

            auto implicitHydrogens(auto atomIdx) const noexcept
            {
                return m_implicitHydrogens.at(atomIdx).get();
            }

            auto isAromaticAtom(auto atomIdx) const noexcept
            {
                return m_aromaticAtoms.at(atomIdx).get();
            }

            auto incident(auto atomIdx) const noexcept
            {
                return m_incident.at(atomIdx).range();
            }

            auto source(auto bondIdx) const noexcept
            {
                return m_sources.at(bondIdx).get();
            }

            auto target(auto bondIdx) const noexcept
            {
                return m_targets.at(bondIdx).get();
            }

            auto order(auto bondIdx) const noexcept
            {
                return m_orders.at(bondIdx).get();
            }

            auto isCyclicBond(auto bondIdx) const noexcept
            {
                return m_cyclicBonds.at(bondIdx).get();
            }

            auto isAromaticBond(auto bondIdx) const noexcept
            {
                return m_aromaticBonds.at(bondIdx).get();
            }

        private:
            VectorObject<ElementList, typename Object::Source> m_elements;
            VectorObject<IsotopeList, typename Object::Source> m_isotopes;
            VectorObject<ChargeList, typename Object::Source> m_charges;
            VectorObject<DegreeList, typename Object::Source> m_degrees;
            VectorObject<ImplicitHydrogensList, typename Object::Source> m_implicitHydrogens;
            VectorObject<AromaticAtomList, typename Object::Source> m_aromaticAtoms;
            VectorObject<SourceList, typename Object::Source> m_sources;
            VectorObject<TargetList, typename Object::Source> m_targets;
            VectorObject<OrderList, typename Object::Source> m_orders;
            VectorObject<CyclicBondList, typename Object::Source> m_cyclicBonds;
            VectorObject<AromaticBondList, typename Object::Source> m_aromaticBonds;
            VectorObject<IncidentListList, typename Object::Source> m_incident;
    };



    template<typename T, typename ...Ts, typename Ptr = const std::byte*>
    void print(ctll::list<T, Ts...>, Ptr data = nullptr, int depth = 0)
    {
        print(T{}, data, depth);
        if constexpr (sizeof...(Ts)) {
            if constexpr (isFixedSize(T{}))
                print(ctll::list<Ts...>{}, data + sizeOf(T{}), depth);
            else
                print(ctll::list<Ts...>{}, data + sizeOf(T{}, data), depth);
        }
    }



    template<typename T, typename Ptr = const std::byte*>
    void print(T, Ptr data = nullptr, int depth = 0)
    {
        constexpr auto format = "{:.<50}  {: <8}";
        constexpr auto format2 = "{: <50}  {: <8}";

        if (!depth)
            fmt::println(format2, "TYPE", "SIZE");


        auto indent = std::string(4 * depth, ' ');

        if constexpr (isValue(T{})) {
            auto name = fmt::format("{}{}<{}>  ", indent, Util::typeName(T{}), Util::typeName(typename T::Type{}));
            fmt::println(format, name, sizeOf(T{}));
        }

        if constexpr (isStruct(T{})) {
            auto name = fmt::format("{}{}", indent, "Struct  ");
            if constexpr (isFixedSize(T{}))
                fmt::println(format, name, sizeOf(T{}));
            else
                fmt::println(format, name, sizeOf(T{}, data));
            print(T::members, data, depth + 1);
        }

        if constexpr (isVector(T{})) {
            auto name = fmt::format("{}{}", indent, "Vector  ");
            fmt::println(format, name, sizeOf(T{}, data));
            print(typename T::Type{}, data, depth + 1);
        }

    }

    template<typename Object>
    class TypeMoleculeObject
    {
        public:
            TypeMoleculeObject(const std::vector<AtomType> &atomTypes,
                               const std::vector<BondType> &bondTypes,
                               const Object &object)
                : m_atomTypes{atomTypes}, m_bondTypes{bondTypes},
                  m_atomTypeIndexes{object.get(AtomTypeIndexList{})}, m_bondTypeIndexes{object.get(BondTypeIndexList{})},
                  m_sources{object.get(SourceList{})}, m_targets{object.get(TargetList{})}
            {
                static_assert(contains(typename Object::Type{}, IncidentListList{}));
                m_incident = object.get(IncidentListList{});

                //print(TypeMolecule{}, object.data());
            }

            auto numAtoms() const noexcept
            {
                return m_atomTypeIndexes.length();
            }

            auto numBonds() const noexcept
            {
                return m_bondTypeIndexes.length();
            }

            auto element(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].element;
            }

            auto isotope(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].isotope;
            }

            auto charge(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].charge;
            }

            auto degree(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].degree;
            }

            auto implicitHydrogens(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].implicitHydrogens;
            }

            auto isAromaticAtom(auto atomIdx) const noexcept
            {
                return m_atomTypes[m_atomTypeIndexes.at(atomIdx).get()].aromatic;
            }

            auto incident(auto atomIdx) const noexcept
            {
                return m_incident.at(atomIdx).range();
            }

            auto source(auto bondIdx) const noexcept
            {
                return m_sources.at(bondIdx).get();
            }

            auto target(auto bondIdx) const noexcept
            {
                return m_targets.at(bondIdx).get();
            }

            auto order(auto bondIdx) const noexcept
            {
                return m_bondTypes[m_bondTypeIndexes.at(bondIdx).get()].order;
            }

            auto isCyclicBond(auto bondIdx) const noexcept
            {
                return m_bondTypes[m_bondTypeIndexes.at(bondIdx).get()].cyclic;
            }

            auto isAromaticBond(auto bondIdx) const noexcept
            {
                return m_bondTypes[m_bondTypeIndexes.at(bondIdx).get()].aromatic;
            }

        private:
            const std::vector<AtomType> &m_atomTypes;
            const std::vector<BondType> &m_bondTypes;
            VectorObject<AtomTypeIndexList, typename Object::Source> m_atomTypeIndexes;
            VectorObject<BondTypeIndexList, typename Object::Source> m_bondTypeIndexes;
            VectorObject<SourceList, typename Object::Source> m_sources;
            VectorObject<TargetList, typename Object::Source> m_targets;
            VectorObject<IncidentListList, typename Object::Source> m_incident;
    };

    template<typename Layout, typename Source>
    class LayoutMolSource
    {
        public:
            LayoutMolSource(std::string_view path) : m_source{path}
            {
                initTypes();
            }

            LayoutMolSource(const Source &source) : m_source{source}
            {
                initTypes();
            }


            auto numMolecules()
            {
                if constexpr (std::is_same_v<Layout, TypeMolecules>)
                    return toObject(Layout(), m_source).get(Vector<TypeMolecule>{}).length();
                else
                    return toObject(Layout(), m_source).length();
            }

            auto read()
            {
                assert (m_index < numMolecules());
                //std::cout << m_index << std::endl;

                if constexpr (contains(Layout{}, AtomList{})) {
                    auto object = toObject(Layout{}, m_source).at(m_index++);
                    return StructMoleculeObject{object};
                } else if constexpr (contains(Layout{}, ElementList{})) {
                    auto object = toObject(Layout{}, m_source).at(m_index++);
                    return ListMoleculeObject{object};
                } else if constexpr (std::is_same_v<Layout, TypeMolecules>) {
                    auto object = toObject(Layout{}, m_source).get(Vector<TypeMolecule>{}).at(m_index++);
                    return TypeMoleculeObject{m_atomTypes, m_bondTypes, object};
                }
                //return toObject(Layout{}, m_source).at(m_index++);
            }

            void reset()
            {
                m_index = 0;
            }

            auto molecules()
            {
                return std::views::iota(static_cast<decltype(numMolecules())>(0), numMolecules()) |
                        std::views::transform([this] (auto) {
                            return read();
                        });
            }

        protected:
            auto initTypes()
            {
                if constexpr (std::is_same_v<Layout, TypeMolecules>) {
                    auto atomTypes = toObject(Layout{}, m_source).get(AtomTypeList{});
                    auto numAtomTypes = atomTypes.length();
                    for (auto i = 0; i < numAtomTypes; ++i) {
                        auto atomType = atomTypes.at(i);
                        m_atomTypes.push_back(AtomType{
                            static_cast<uint8_t>(atomType.get(ElementValue{}).get()),
                            static_cast<uint8_t>(atomType.get(IsotopeValue{}).get()),
                            static_cast<int8_t>(atomType.get(ChargeValue{}).get()),
                            static_cast<uint8_t>(atomType.get(DegreeValue{}).get()),
                            static_cast<uint8_t>(atomType.get(ImplicitHydrogensValue{}).get()),
                            atomType.get(AromaticAtomValue{}).get()
                        });
                    }

                    auto bondTypes = toObject(Layout{}, m_source).get(BondTypeList{});
                    auto numBondTypes = bondTypes.length();
                    for (auto i = 0; i < numBondTypes; ++i) {
                        auto bondType = bondTypes.at(i);
                        m_bondTypes.push_back(BondType{
                            static_cast<uint8_t>(bondType.get(OrderValue{}).get()),
                            bondType.get(CyclicBondValue{}).get(),
                            bondType.get(AromaticBondValue{}).get()
                        });
                    }
                }
            }

            Source m_source;
            std::vector<AtomType> m_atomTypes;
            std::vector<BondType> m_bondTypes;
            std::size_t m_index = 0;
    };

    template<typename Layout>
    using FileStreamMolSource = LayoutMolSource<Layout, FileStreamSource>;

    template<typename Layout>
    using MemoryMappedMolSource = LayoutMolSource<Layout, MemoryMappedSource>;

    template<typename Layout>
    using InMemoryMolSource = LayoutMolSource<Layout, InMemorySource>;

    template<typename Layout>
    using BytePtrMolSource = LayoutMolSource<Layout, PtrSource>;

} // namespace Kitimar::CTLayout
