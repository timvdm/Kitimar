#pragma once

#include <Kitimar/Util/Util.hpp>
#include <Kitimar/Molecule/Molecule.hpp>
#include <Kitimar/CTLayout/Value.hpp>
#include <Kitimar/CTLayout/Struct.hpp>
#include <Kitimar/CTLayout/Vector.hpp>
#include <Kitimar/CTLayout/Sink.hpp>

namespace Kitimar::CTLayout {

    using SizeT2 = uint32_t;

    // Object-Oriented design using Struct
    struct ElementValue : Value<uint8_t> {};
    struct IsotopeValue : Value<uint8_t> {};
    struct ChargeValue : Value<int8_t> {};
    struct DegreeValue : Value<uint8_t> {};
    struct ImplicitHydrogensValue : Value<uint8_t> {};
    struct AromaticAtomValue : Value<bool> {};
    struct AtomStruct : Struct<ElementValue, IsotopeValue, ChargeValue, DegreeValue, ImplicitHydrogensValue, AromaticAtomValue> {};
    struct AtomList : Vector<AtomStruct> {};

    struct SourceValue : Value<AtomList::Index> {};
    struct TargetValue : Value<AtomList::Index> {};
    struct OrderValue : Value<uint8_t> {};
    struct CyclicBondValue : Value<bool> {};
    struct AromaticBondValue : Value<bool> {};
    struct BondStruct : Struct<SourceValue, TargetValue, OrderValue, AromaticBondValue, CyclicBondValue> {};
    struct BondList : Vector<BondStruct> {};

    // Data-Oriented design using no Struct    
    struct ElementList : Vector<ElementValue> {};
    struct IsotopeList : Vector<IsotopeValue> {};
    struct ChargeList : Vector<ChargeValue> {};
    struct DegreeList : Vector<DegreeValue> {};
    struct ImplicitHydrogensList : Vector<ImplicitHydrogensValue> {};
    //struct AromaticAtomList : BitVector<NumAtoms> {};
    struct AromaticAtomList : Vector<AromaticAtomValue> {};

    struct SourceList : Vector<SourceValue> {};
    struct TargetList : Vector<TargetValue> {};
    struct OrderList : Vector<OrderValue> {};
    //struct AromaticBondList : BitVector<NumBonds> {};
    struct CyclicBondList : Vector<CyclicBondValue> {};
    struct AromaticBondList : Vector<AromaticBondValue> {};


    // Type-Oriented design    
    struct AtomTypeList : Vector<AtomStruct> {};
    struct AtomTypeIndex : Value<uint16_t> {};
    struct AtomTypeIndexList : Vector<AtomTypeIndex> {};

    struct BondTypeList : Vector<BondStruct> {};
    struct BondTypeIndex : Value<uint8_t> {};
    struct BondTypeIndexList : Vector<BondTypeIndex> {};

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
    struct IncidentIndex : Value<SizeT2> {};
    struct IncidentList : Vector<IncidentIndex> {};
    struct IncidentListList : Vector<IncidentList> {};


    // adjacent list (i.e. nbrs of atom)
    struct AdjacentIndex : Value<SizeT2> {};
    struct AdjacentList : Vector<AdjacentIndex> {};
    struct AdjacentListList : Vector<AdjacentList> {};



    struct StructMolecule : Struct<AtomList, BondList> {};
    struct StructMoleculeIncident : Struct<AtomList, BondList, IncidentListList> {};
    struct StructMoleculeAdjacent : Struct<AtomList, BondList, AdjacentListList> {};
    struct StructMoleculeIncidentAdjacent : Struct<AtomList, BondList,
                                                   IncidentListList,
                                                   AdjacentListList> {};

    /*
    struct ListMolecule : Layout<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList,
                                  SourceList, TargetList, OrderList> {};
    struct ListMoleculeIncident : Layout<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList,
                                          SourceList, TargetList, OrderList,
                                          IncidentBeginList, IncidentList> {};
    struct ListMoleculeAdjacent : Layout<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList,
                                          SourceList, TargetList, OrderList,
                                          AdjacentBeginList, AdjacentList> {};
    struct ListMoleculeIncidentAdjacent : Layout<ElementList, IsotopeList, ChargeList, DegreeList, ImplicitHydrogensList,
                                                  SourceList, TargetList, OrderList,
                                                  IncidentBeginList, IncidentList,
                                                  AdjacentBeginList, AdjacentList> {};
    */



    struct AtomTypeMolecule : Struct<AtomTypeIndexList, BondList, IncidentListList> {};


    //
    // Molecule
    //

    template<typename MolObj>
    constexpr auto num_atoms(const MolObj &mol) noexcept
    {
        if constexpr (contains(typename MolObj::Type{}, AtomList{}))
            return mol.get(AtomList{}).length();
        else
            return mol.get(ElementList{}).length();
    }

    template<typename MolObj>
    constexpr auto num_bonds(const MolObj &mol) noexcept
    {
        if constexpr (contains(typename MolObj::Type{}, BondList{}))
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


    // Atom

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
        return get_atom_property(mol, atom, ElementValue{}, ElementList{});
    }

    constexpr auto get_isotope(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_atom_property(mol, atom, IsotopeValue{}, IsotopeList{});
    }

    constexpr auto get_charge(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_atom_property(mol, atom, ChargeValue{}, ChargeList{});
    }

    constexpr auto get_degree(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        return get_atom_property(mol, atom, DegreeValue{}, DegreeList{});
    }

    constexpr auto get_implicit_hydrogens(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        return get_atom_property(mol, atom, ImplicitHydrogensValue{}, ImplicitHydrogensList{});
    }

    constexpr auto is_aromatic_atom(const auto &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));        
        return get_atom_property(mol, atom, AromaticAtomValue{}, AromaticAtomList{});
    }

    template<typename MolObj>
    constexpr auto get_bonds(const MolObj &mol, auto atom) noexcept
    {
        assert(atom < num_atoms(mol));
        if constexpr (contains(typename MolObj::Type{}, IncidentListList{})) {
            auto bonds = mol.get(IncidentListList{}).at(atom);
            auto ptr = reinterpret_cast<const IncidentIndex::Type*>(static_cast<const std::byte*>(bonds.data()));
            return std::ranges::subrange(ptr, ptr + bonds.length());
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
        if constexpr (contains(typename MolObj::Type{}, AdjacentListList{})) {
            auto nbrs = mol.get(AdjacentListList{}).at(atom);
            auto ptr = reinterpret_cast<const AdjacentIndex::Type*>(static_cast<const std::byte*>(nbrs.data()));
            return std::ranges::subrange(ptr, ptr + nbrs.length());
        } else {
            return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
                return get_nbr(mol, bond, atom);
            });
        }
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
        return get_bond_property(mol, bond, SourceValue{}, SourceList{});
    }

    constexpr auto get_target(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return get_bond_property(mol, bond, TargetValue{}, TargetList{});
    }

    constexpr auto get_order(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));        
        return get_bond_property(mol, bond, OrderValue{}, OrderList{});
    }

    constexpr auto is_cyclic_bond(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));        
        return get_bond_property(mol, bond, CyclicBondValue{}, CyclicBondList{});
    }

    constexpr auto is_aromatic_bond(const auto &mol, auto bond) noexcept
    {
        assert(bond < num_bonds(mol));        
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

    template<typename Layout>
    void serializeIncident(Molecule auto &mol, auto &writer)
    {
        if constexpr (contains(Layout{}, IncidentListList{})) {
            auto incidentList = writer.get(IncidentListList{});
            incidentList.setLength(num_atoms(mol));
            for (auto atom : get_atoms(mol)) {
                auto incident = incidentList.at(get_index(mol, atom));
                incident.setLength(get_degree(mol, atom));
                int i = 0;
                for (auto bond : get_bonds(mol, atom))
                    incident.at(i++).write(static_cast<IncidentIndex::Type>(get_index(mol, bond)));
            }
        }
    }

    template<typename Layout>
    void serializeAdjacent(Molecule auto &mol, auto &writer)
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
    void serializeAtomList(Molecule auto &mol, auto &writer)
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
        }
    }

    template<typename Layout>
    void serializeBondList(Molecule auto &mol, auto &writer)
    {
        if constexpr (contains(Layout{}, BondList{})) {
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
        }
    }

    template<typename Layout>
    auto findTypes(auto &source)
    {
        std::vector<AtomType> atomTypes;
        std::vector<BondType> bondTypes;

        if constexpr (contains(Layout{}, AtomTypeList{})) { // FIXME or BondTypeList
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

    /*
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
    void serializeAtomTypeList(Molecule auto &mol, std::byte *data,
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
    void writeMolecule(Molecule auto &mol, auto &writer,
                           const std::vector<AtomType> &atomTypes = {})
    {
        //serializeAtomTypeList<Layout>(mol, data, atomTypes);

        serializeAtomList<Layout>(mol, writer);
        serializeBondList<Layout>(mol, writer);

        serializeIncident<Layout>(mol, writer);
        serializeAdjacent<Layout>(mol, writer);

    }


    template<typename Layout, typename Sink>
    void serializeMolecule(Molecule auto &mol, Sink &sink,
                           const std::vector<AtomType> &atomTypes = {})
    {            
        auto writer = toWriter(Layout{}, sink);
        writeMolecule<Layout>(mol, writer, atomTypes);
    }


/*



    template<typename Layout>
    constexpr auto moleculeSize(SizeT numAtoms, SizeT numBonds) noexcept
    {
        return Layout::size({numAtoms, numBonds});
    }

*/

    // FIXME: reaneme FileStreamMolSink
    template<typename Layout>
    void serializeMolSource(auto &molSource, std::string_view filename)
    {
        FileStreamSink sink{filename};
        auto writer = toWriter(Layout{}, sink);

        std::cout << "filename: " << filename << std::endl;

        fmt::println("counting molecules....");
        std::cout << std::endl;
        auto n =  molSource.numMolecules();
        fmt::println("# molecules: {}", n);


        writer.setLength(n);
        for (auto i = 0; i < n; ++i) {
            //if (i % 1000 == 0)
                std::cout << i << std::endl;
            auto molWriter = writer.at(i);
            auto mol = molSource.read();
            writeMolecule<typename Layout::Type>(mol, molWriter);

        }

        //auto [atomTypes, bondTypes] = findTypes<Layout>(molSource);

        //std::cout << "# atom types: " << atomTypes.size() << std::endl;
        //std::cout << "# bond types: " << bondTypes.size() << std::endl;


        //for (auto mol : molSource.molecules()) {
            //buffer.resize(moleculeSize<Layout>(num_atoms(mol), num_bonds(mol)));
            //serializeMolecule<Layout>(mol, buffer.data(), atomTypes);
            //sink.writeObject(buffer);
            //if (sink.size() % 10000 == 0) std::cout << sink.size() << std::endl;
        //}
    }

    template<typename Layout, typename Source>
    class LayoutMolSource
    {
        public:
            LayoutMolSource(std::string_view path) : m_source{path}
            {
            }

            auto numMolecules()
            {
                return toObject(Layout(), m_source).length();
            }

            auto read()
            {
                assert (m_index < numMolecules());
                std::cout << m_index << std::endl;
                return toObject(Layout(), m_source).at(m_index++);
            }

            void reset()
            {
                m_index = 0;
            }

            auto molecules()
            {
                return std::views::iota(static_cast<Layout::Length>(0), numMolecules()) |
                        std::views::transform([this] (auto) {
                            return read();
                        });
            }

        protected:
            Source m_source;
            std::size_t m_index = 0;
    };

    template<typename Layout>
    using FileStreamMolSource = LayoutMolSource<Layout, FileStreamSource>;

    template<typename Layout>
    using MemoryMappedMolSource = LayoutMolSource<Layout, MemoryMappedSource>;

    /*
    template<typename Layout>
    class FileStreamMolSource
    {
        public:
            FileStreamMolSource(std::string_view path) : m_source{path}
            {
            }

            auto numMolecules()
            {
                return toObject(Layout(), m_source).length();
            }

            auto read()
            {
                assert (m_index < numMolecules());
                std::cout << m_index << std::endl;
                return toObject(Layout(), m_source).at(m_index++);
            }

            void reset()
            {
                m_index = 0;
            }

            auto molecules()
            {
                return std::views::iota(static_cast<Layout::Length>(0), numMolecules()) |
                        std::views::transform([this] (auto) {
                            return read();
                        });
            }




        protected:
            //Source m_source;
            FileStreamSource m_source;
            std::size_t m_index = 0;
    };

    template<typename Layout>
    class MemoryMappedMolSource
    {
        public:
            MemoryMappedMolSource(std::string_view path) : m_source{path}
            {
            }

            auto numMolecules()
            {
                return toObject(Layout(), m_source).length();
            }

            auto read()
            {
                assert (m_index < numMolecules());
                if (m_index % 1000 == 0)
                    std::cout << m_index << std::endl;
                return toObject(Layout(), m_source).at(m_index++);
            }

            void reset()
            {
                m_index = 0;
            }

            auto molecules()
            {
                return std::views::iota(static_cast<Layout::Length>(0), numMolecules()) |
                        std::views::transform([this] (auto) {
                            return read();
                        });
            }




        protected:
            //Source m_source;
            MemoryMappedSource m_source;
            std::size_t m_index = 0;
    };
    */

    /*
    template<typename Layout>
    class FileStreamMolSource : public LayoutMolSource<Layout, FileStreamSource>
    {
        public:
            FileStreamMolSource(std::string_view path) : m_source()
            {
            }
        private:
    }
    */



} // namespace Kitimar::CTLayout
