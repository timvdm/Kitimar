#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>
//#include <Kitimar/Util/Util.hpp>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;




template<template<int> class Primitive, int Begin, int End>
constexpr auto addPrimitives(auto primitives) noexcept
{
    if constexpr (Begin == End)
        return primitives;
    else
        return ctll::push_front(Primitive<Begin>{}, addPrimitives<Primitive, Begin + 1, End>(primitives));
}


constexpr auto makeAtomPrimitives() noexcept
{
    auto primitives1 = ctll::list<
        AromaticAtom<5>,
        AromaticAtom<6>,
        AromaticAtom<7>,
        AromaticAtom<8>,
        AromaticAtom<15>,
        AromaticAtom<16>,
        AromaticAtom<33>,
        AromaticAtom<34>,
        AnyAtom,
        AnyAliphatic,
        AnyAromatic,
        Cyclic,
        Acyclic
    >{};

    return addPrimitives<AliphaticAtom, 1, 105>(
                addPrimitives<Element, 1, 105>(
                    addPrimitives<Isotope, 0, 260>(
                        addPrimitives<Degree, 0, 10>(
                            addPrimitives<Valence, 0, 10>(
                                addPrimitives<Connectivity, 0, 10>(
                                    addPrimitives<TotalH, 0, 10>(
                                        addPrimitives<ImplicitH, 0, 10>(
                                            addPrimitives<RingCount, 0, 10>(
                                                addPrimitives<RingSize, 0, 10>(
                                                    addPrimitives<RingConnectivity, 0, 10>(
                                                        addPrimitives<Charge, -15, 16>(primitives1)
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
    );
}

constexpr auto makeBondPrimitives() noexcept
{
    return ctll::list<
        ImplicitBond,
        AnyBond,
        AromaticBond,
        RingBond,
        BondOrder<1>,
        BondOrder<2>,
        BondOrder<3>,
        BondOrder<4>
    >{};
}




class PrimitiveFequency
{
    public:
        static constexpr inline auto atomPrimitives = makeAtomPrimitives();
        static constexpr inline auto bondPrimitives = makeBondPrimitives();


        void addAtom(auto &mol, const auto &atom)
        {
            ++numAtoms;
            updateAtomPrimitive(mol, atom, atomPrimitives);
        }

        void addBond(auto &mol, const auto &bond)
        {
            ++numBonds;
            updateBondPrimitive(mol, bond, bondPrimitives);
        }

        void addMolecule(auto &mol)
        {
            ++numMolecules;
            for (auto atom : get_atoms(mol))
                addAtom(mol, atom);
            for (auto bond : get_bonds(mol))
                addBond(mol, bond);
        }

        std::array<double, ctll::size(atomPrimitives)> atomFrequency() const
        {
            std::array<double, ctll::size(atomPrimitives)> freq;
            for (auto i = 0; i < freq.size(); ++i)
                freq[i] = static_cast<double>(atomMatches[i]) / numAtoms;
            return freq;
        }

        std::array<double, ctll::size(bondPrimitives)> bondFrequency() const
        {
            std::array<double, ctll::size(bondPrimitives)> freq;
            for (auto i = 0; i < freq.size(); ++i)
                freq[i] = static_cast<double>(bondMatches[i]) / numBonds;
            return freq;
        }

        std::array<std::size_t, ctll::size(atomPrimitives)> atomMatches = {};
        std::array<std::size_t, ctll::size(bondPrimitives)> bondMatches = {};
        std::size_t numMolecules = 0;
        std::size_t numAtoms = 0;
        std::size_t numBonds = 0;

    private:
        template<int I = 0>
        void updateAtomPrimitive(auto &mol, const auto &atom, auto primitives)
        {
            if constexpr (ctll::empty(primitives))
                return;
            else {
                auto [primitive, tail] = ctll::pop_and_get_front(primitives);
                if (matchAtomExpr(mol, atom, primitive))
                    ++atomMatches[I];
                updateAtomPrimitive<I + 1>(mol, atom, tail);
            }
        }

        template<int I = 0>
        void updateBondPrimitive(auto &mol, const auto &bond, auto primitives)
        {
            if constexpr (ctll::empty(primitives))
                return;
            else {
                auto [primitive, tail] = ctll::pop_and_get_front(primitives);
                if (matchBondExpr(mol, bond, primitive))
                    ++bondMatches[I];
                updateBondPrimitive<I + 1>(mol, bond, tail);
            }
        }




};

std::string toStringHelper(std::string_view name, int value)
{
    std::stringstream ss;
    ss << name << "<" << value << ">";
    return ss.str();
}

std::string toString(AnyAtom) { return "AnyAtom"; }
std::string toString(AnyAliphatic) { return "AnyAliphatic"; }
std::string toString(AnyAromatic) { return "AnyAromatic"; }
std::string toString(Cyclic) { return "Cyclic"; }
std::string toString(Acyclic) { return "Acyclic"; }
template<int N> std::string toString(AliphaticAtom<N>) { return toStringHelper("AliphaticAtom", N); }
template<int N> std::string toString(AromaticAtom<N>) { return toStringHelper("AromaticAtom", N); }
template<int N> std::string toString(Isotope<N>) { return toStringHelper("Isotope", N); }
template<int N> std::string toString(Element<N>) { return toStringHelper("Element", N); }
template<int N> std::string toString(Degree<N>) { return toStringHelper("Degree", N); }
template<int N> std::string toString(Valence<N>) { return toStringHelper("Valence", N); }
template<int N> std::string toString(Connectivity<N>) { return toStringHelper("Connectivity", N); }
template<int N> std::string toString(TotalH<N>) { return toStringHelper("TotalH", N); }
template<int N> std::string toString(ImplicitH<N>) { return toStringHelper("ImplicitH", N); }
template<int N> std::string toString(RingCount<N>) { return toStringHelper("RingCount", N); }
template<int N> std::string toString(RingSize<N>) { return toStringHelper("RingSize", N); }
template<int N> std::string toString(RingConnectivity<N>) { return toStringHelper("RingConnectivity", N); }
template<int N> std::string toString(Charge<N>) { return toStringHelper("Charge", N); }

std::string toString(AnyBond) { return "AnyBond"; }
std::string toString(ImplicitBond) { return "ImplicitBond"; }
std::string toString(AromaticBond) { return "AromaticBond"; }
std::string toString(RingBond) { return "RingBond"; }
template<int N> std::string toString(BondOrder<N>) { return toStringHelper("BondOrder", N); }







template<typename T, auto N, int I = 0>
void printPrimitiveProperty(auto primitives, const std::array<T, N> &properties)
{
    if constexpr (!ctll::empty(primitives)) {
        auto [primitive, tail] = ctll::pop_and_get_front(primitives);
        std::cout << toString(primitive) << ": " << properties[I] << std::endl;
        printPrimitiveProperty<T, N, I + 1>(tail, properties);
    }
}

template<auto N, int I = 0>
void writePrimitiveFrequencyCodeFunction(std::ostream &os, auto primitives, const std::array<double, N> &frequency)
{
    if constexpr (!ctll::empty(primitives)) {
        auto [primitive, tail] = ctll::pop_and_get_front(primitives);
        os << "    consteval double expressionFrequency(" << toString(primitive) << ") noexcept { return " << frequency[I] << "; }\n";
        writePrimitiveFrequencyCodeFunction<N, I + 1>(os, tail, frequency);
    }
}

void writePrimitiveFrequencyCodeFile(std::ostream &os, std::string_view filename, const PrimitiveFequency &freq)
{
    os << "#pragma once\n";
    os << "\n";
    os << "#include <Kitimar/CTSmarts/AST.hpp>\n";
    os << "\n";
    os << "namespace Kitimar::CTSmarts {\n";
    os << "\n";
    os << "    // Input file: " << filename << "\n";
    os << "    //\n";
    os << "    //     # molecules: " << freq.numMolecules << "\n";
    os << "    //     # atoms:     " << freq.numAtoms << "\n";
    os << "    //     # bonds:     " << freq.numBonds << "\n";
    os << "\n";
    os << "    // Atom Primitives\n";
    os << "\n";
    writePrimitiveFrequencyCodeFunction(os, freq.atomPrimitives, freq.atomFrequency());
    os << "\n";
    os << "    // Bond Primitives\n";
    os << "\n";
    writePrimitiveFrequencyCodeFunction(os, freq.bondPrimitives, freq.bondFrequency());
    os << "\n";
    os << "} // namespace Kitimar::CTSmarts\n";
}


template<auto N, int I = 0>
void writePrimitiveFrequencyData(std::ostream &os, auto primitives, const std::array<double, N> &frequency)
{
    if constexpr (!ctll::empty(primitives)) {
        auto [primitive, tail] = ctll::pop_and_get_front(primitives);
        os << toString(primitive) << '\t' << frequency[I] << '\n';
        writePrimitiveFrequencyData<N, I + 1>(os, tail, frequency);
    }
}

void writePrimitiveFrequencyDataFile(std::ostream &os, std::string_view filename, const PrimitiveFequency &freq)
{
    os << freq.numMolecules << '\t' << freq.numAtoms << '\t' << freq.numBonds << '\n';
    writePrimitiveFrequencyData(os, freq.atomPrimitives, freq.atomFrequency());
    writePrimitiveFrequencyData(os, freq.bondPrimitives, freq.bondFrequency());
}





int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <code or data> <molecule file>" << std::endl;
        return -1;
    }

    PrimitiveFequency freq;

    auto p = freq.atomPrimitives;

    auto filename = argv[2];
    OpenBabelSmilesMolSource source{filename};


    for (const auto &mol : source.molecules())
        freq.addMolecule(const_cast<OpenBabel::OBMol&>(mol));

    std::string output{argv[1]};
    if (output == "code")
        writePrimitiveFrequencyCodeFile(std::cout, filename, freq);
    else if (output == "data")
        writePrimitiveFrequencyDataFile(std::cout, filename, freq);
    else {
        std::cerr << "Invalid output type: " << output << std::endl;
        return -1;
    }



}
