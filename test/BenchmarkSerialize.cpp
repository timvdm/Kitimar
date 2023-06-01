#include <Kitimar/Util/Util.hpp>
#include <Kitimar/CTLayout/CTLayout.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>

#ifdef KITIMAR_WITH_RDKIT
#include <Kitimar/RDKit/RDKit.hpp>
#include <GraphMol/Substruct/SubstructMatch.h>
#endif

#include "TestData.hpp"
#include "Benchmark.hpp"

#include <openbabel/parsmart.h>


#include <algorithm>


using namespace Kitimar;
using namespace Kitimar::CTLayout;
using namespace Kitimar::CTSmarts;

//
// Generic single match using callback
//

// OpenBabel

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesOpenBabel(auto &source, Callback callback)
{
    OpenBabel::OBSmartsPattern smarts;
    smarts.Init(Util::toString(SMARTS));
    for (auto mol : source.molecules())
        callback(mol, smarts.Match(mol));
}

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesOpenBabelSmiles(Callback callback)
{
    OpenBabelSmilesMolSource source{chembl_smi_filename()};
    matchesOpenBabel<SMARTS>(source, std::forward<Callback>(callback));
}

// RDKit

#ifdef KITIMAR_WITH_RDKIT

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesRDKit(auto &source, Callback callback)
{
    std::unique_ptr<RDKit::RWMol> smarts{RDKit::SmartsToMol(Util::toString(SMARTS))};
    RDKit::MatchVectType res;
    for (const auto &mol : source.molecules())
        callback(mol, RDKit::SubstructMatch(*mol, *smarts, res));
}

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesRDKitSmiles(Callback callback)
{
    RDKitSmilesMolSource source{chembl_smi_filename()};
    matchesRDKit<SMARTS>(source, std::forward<Callback>(callback));
}

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesRDKitPickle(Callback callback)
{
    RDKitPickleMolSource source{chembl_rdkit_filename()};
    matchesRDKit<SMARTS>(source, std::forward<Callback>(callback));
}

#endif // KITIMAR_WITH_RDKIT

// Kitimar

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesKitimar(auto &source, Callback callback)
{
    SingleIsomorphism<SMARTS> smarts{};
    for (auto mol : source.molecules())
        callback(mol, smarts.match(mol));
}


template<ctll::fixed_string SMARTS, typename Layout, typename Callback>
auto matchesKitimarFileStream(Callback callback)
{
    FileStreamMolSource<Layout> source{chembl_serialized_filename(Layout{})};
    matchesKitimar<SMARTS>(source, std::forward<Callback>(callback));
}

template<ctll::fixed_string SMARTS, typename Layout, typename Callback>
auto matchesKitimarMemoryMapped(Callback callback)
{    
    MemoryMappedSource source{chembl_serialized_filename(Layout{})};
    BytePtrMolSource<Layout> molSource{source.toPtrSource()};
    matchesKitimar<SMARTS>(molSource, std::forward<Callback>(callback));
}

//
// Count of matching molecules
//

template<ctll::fixed_string SMARTS>
Benchmark countMatchesOpenBabelSmiles()
{
    Util::Timer timer;
    auto n = 0;
    matchesOpenBabelSmiles<SMARTS>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"OpenBabel SMILES", n, timer.elapsed()};
}

#ifdef KITIMAR_WITH_RDKIT

template<ctll::fixed_string SMARTS>
Benchmark countMatchesRDKitSmiles()
{
    Util::Timer timer;
    auto n = 0;
    matchesRDKitSmiles<SMARTS>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"RDKit SMILES", n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS>
Benchmark countMatchesRDKitPickle()
{
    Util::Timer timer;
    auto n = 0;
    matchesRDKitPickle<SMARTS>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"RDKit Pickle", n, timer.elapsed()};
}

#endif // KITIMAR_WITH_RDKIT

template<ctll::fixed_string SMARTS, typename Layout>
Benchmark countMatchesKitimarFileStream()
{
    Util::Timer timer;
    auto n = 0;
    matchesKitimarFileStream<SMARTS, Layout>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"Kitimar FileStream", n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS, typename Layout>
Benchmark countMatchesKitimarInMemory()
{    
    InMemorySource source{chembl_serialized_filename(Layout{})};
    BytePtrMolSource<Layout> molSource{source.toPtrSource()};
    Util::Timer timer;
    auto n = 0;
    matchesKitimar<SMARTS>(molSource, [&n] (auto &mol, bool match) { if (match) ++n; });
    return {"Kitimar InMemory", n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS, typename Layout>
Benchmark countMatchesKitimarMemoryMapped()
{
    Util::Timer timer;
    auto n = 0;
    matchesKitimarMemoryMapped<SMARTS, Layout>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"Kitimar MemoryMapped " + Util::typeName(Layout{}), n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS, ctll::fixed_string ...Tail>
void countMatches(std::vector<BenchmarkGroup> &groups)
{
    std::cout << "Count matches: " << Util::toString(SMARTS) << std::endl;
    groups.push_back({ Util::toString(SMARTS), {
        //countMatchesOpenBabelSmiles<SMARTS>(),
        //Benchmark{"OpenBabel SMILES"},
#ifdef KITIMAR_WITH_RDKIT
        //countMatchesRDKitSmiles<SMARTS>(),
        //countMatchesRDKitPickle<SMARTS>(),
        //Benchmark{"RDKit SMILES"},
        //Benchmark{"RDKit Pickle"},
#endif // KITIMAR_WITH_RDKIT
        //countMatchesKitimarFileStream<SMARTS, Vector<StructMoleculeIncident>>(),
        countMatchesKitimarMemoryMapped<SMARTS, StructMolecules>(),
        countMatchesKitimarMemoryMapped<SMARTS, ListMolecules>(),
        countMatchesKitimarMemoryMapped<SMARTS, TypeMolecules>(),
        //countMatchesKitimarInMemory<SMARTS, Vector<StructMoleculeIncident>>()
    }});
    if constexpr (sizeof...(Tail))
        countMatches<Tail...>(groups);
}


template<ctll::fixed_string ...SMARTS>
BenchmarkGroups countMatches()
{
    BenchmarkGroups result{"Count Single Matches"};
    countMatches<SMARTS...>(result.groups);
    return result;
}



//
// Deserialize
//

/*
template<typename Layout>
auto deserializeKitimar()
{
    Util::Timer timer;
    InMemorySource<Layout> source{chembl_serialized_filename(Layout{})};
    return Benchmark{"Kitimar Deserialize", source.size(), timer.elapsed()};
}
*/















int main()
{
    //auto benchmark1 = deserializeKitimar<StructMoleculeIncident>();
    //std::cout << benchmark1.toString() << std::endl;



    auto groups = countMatches<
        //"*1~*~*~*~*~*~1"
        //"c1ccccc1CCCCc1ccccc1"
        /**/
        "*",
        "*1~*~*~*~*~*~1",
        "c1ccccc1-c2ccccc2",
        "c1ccccc1CCCCc1ccccc1",
        "c1cc2c(cc1)cccc2",
        "Clc1ccccc1",
        "BrCCCCCCCCC",        
        "O1CCOC12CCNCC2"
        /**/

    >();
    std::cout << groups.toString() << std::endl;

    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "[nH]1ccc2c1cccc2">));        //
    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "[nH]">));        //
    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "c1cc(=O)cc[nH]1">));        //

    for (auto &[name, total] : groups.totalElapsed())
        std::cout << Util::pad(name, 50) << total << '\n';
    std::cout << std::endl;

}
