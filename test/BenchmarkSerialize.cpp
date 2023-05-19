#include <Kitimar/Util/Util.hpp>
#include <Kitimar/CTLayout/CTLayout.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/RDKit/RDKit.hpp>

#include "TestData.hpp"
#include "Benchmark.hpp"

#include <openbabel/parsmart.h>

#include <GraphMol/Substruct/SubstructMatch.h>

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

// Kitimar

template<ctll::fixed_string SMARTS, typename Callback>
auto matchesKitimar(auto &source, Callback callback)
{
    SingleIsomorphism<SMARTS> smarts{};
    for (auto mol : source.objects())
        callback(mol, smarts.match(mol));
}


template<ctll::fixed_string SMARTS, typename Layout, typename Callback>
auto matchesKitimarFileStream(Callback callback)
{
    FileStreamSource<Layout> source{chembl_serialized_filename(Layout{})};
    matchesKitimar<SMARTS>(source, std::forward<Callback>(callback));
}

template<ctll::fixed_string SMARTS, typename Layout, typename Callback>
auto matchesKitimarMemoryMapped(Callback callback)
{
    MemoryMappedSource<Layout> source{chembl_serialized_filename(Layout{})};
    matchesKitimar<SMARTS>(source, std::forward<Callback>(callback));
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
    InMemorySource<Layout> source{chembl_serialized_filename(Layout{})};
    Util::Timer timer;
    auto n = 0;
    matchesKitimar<SMARTS>(source, [&n] (auto &mol, bool match) { if (match) ++n; });
    return {"Kitimar InMemory", n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS, typename Layout>
Benchmark countMatchesKitimarMemoryMapped()
{
    Util::Timer timer;
    auto n = 0;
    matchesKitimarMemoryMapped<SMARTS, Layout>([&n] (auto &mol, bool match) { if (match) ++n; });
    return {"Kitimar MemoryMapped", n, timer.elapsed()};
}

template<ctll::fixed_string SMARTS, ctll::fixed_string ...Tail>
void countMatches(std::vector<BenchmarkGroup> &groups)
{
    groups.push_back({ Util::toString(SMARTS), {
        //countMatchesOpenBabelSmiles<SMARTS>(),
        //countMatchesRDKitSmiles<SMARTS>(),
        //countMatchesRDKitPickle<SMARTS>(),
        Benchmark{"OpenBabel SMILES"},
        Benchmark{"RDKit SMILES"},
        Benchmark{"RDKit Pickle"},
        countMatchesKitimarFileStream<SMARTS, StructMoleculeIncident>(),
        countMatchesKitimarInMemory<SMARTS, StructMoleculeIncident>(),
        countMatchesKitimarMemoryMapped<SMARTS, StructMoleculeIncident>()
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

template<typename Layout>
auto deserializeKitimar()
{
    Util::Timer timer;
    InMemorySource<Layout> source{chembl_serialized_filename(Layout{})};
    return Benchmark{"Kitimar Deserialize", source.size(), timer.elapsed()};
}
















int main()
{




    //auto benchmark1 = deserializeKitimar<StructMoleculeIncident>();
    //std::cout << benchmark1.toString() << std::endl;




    //auto groups = countMatches<"*1~*~*~*~*~*~1", "BrCCCCCCCCC">();
    auto groups = countMatches<"*1~*~*~*~*~*~1">();
    //auto groups = countMatches<"*">();
    std::cout << groups.toString() << std::endl;



//    FileStream      : 17.6198s
//    MemoryMapped    : 16.0925s
//    OpenBabel       : 612.724s    10.21m
//    RDKit           : 1394.68s    23.24m

    //std::vector<Benchmark> benchmarks;



    /*
    benchmark<"*1~*~*~*~*~*~1">(benchmarks);
    benchmark<"c1ccccc1-c2ccccc2">(benchmarks);
    benchmark<"c1ccccc1CCCCc1ccccc1">(benchmarks);
    benchmark<"c1cc2c(cc1)cccc2">(benchmarks);
    benchmark<"Clc1ccccc1">(benchmarks);
    benchmark<"BrCCCCCCCCC">(benchmarks);
    benchmark<"O1CCOC12CCNCC2">(benchmarks);
    */



    for (auto &[name, total] : groups.totalElapsed())
        std::cout << Util::pad(name) << total << '\n';


    /*
    std::map<std::string, std::vector<Benchmark>> byName;
    for (auto &benchmark : benchmarks)
        byName[benchmark.name].push_back(benchmark);


    for (auto &[name, nameBenchmarks] : byName) {
        Util::Timer::Duration total = {};
        for (auto &benchmark : nameBenchmarks)
            total += benchmark.elapsed;

        std::cout << Util::pad(name) << total << '\n';
    }
    */

    std::cout << std::endl;

}








/*
template<typename F>
auto benchmark_OpenBabel(const std::string &filename, F f)
{
    auto i = 0;
    readMolecules(filename, [&i, &f] (auto &mol) {
        if (++i % 10000 == 0) std::cout << i << std::endl;
        f(mol);
        return true;
    });
}
*/

/* RECENT
template<typename Layout, typename F>
auto benchmark_Kitimar(F f)
{
    auto filename = chembl_serialized_filename(Layout{});
    std::vector<std::byte> data;
    std::ifstream ifs(filename, std::ios_base::binary | std::ios_base::out); // FIXME : in??
    auto i = 0;
    while (ifs) {
        //if (++i % 1000 == 0) std::cout << i << std::endl;
        auto [ok, mol] = deserialize<Layout>(ifs, data);
        if (!ok)
            break;
        f(mol);
    }
}
*/

/*
auto numAtoms_OpenBabel(const std::string &filename)
{
    std::size_t numAtoms = 0;
    benchmark_OpenBabel(filename, [&numAtoms] (auto &mol) {
        numAtoms += num_atoms(mol);
    });
    std::cout << "# atoms: " << numAtoms << std::endl;
}


auto numAtoms_OpenBabel_Sdf()
{
    numAtoms_OpenBabel(chembl_sdf_filename());
}

auto numAtoms_OpenBabel_Smi()
{
    numAtoms_OpenBabel(chembl_smi_filename());
}
*/

/* RECENT
template<typename Layout>
auto numAtoms_Kitimar()
{
    std::size_t n = 0;
    benchmark_Kitimar<Layout>([&n] (auto &mol) {
        n += num_atoms(mol);
    });
    std::cout << "# atoms: " << n << std::endl;
}

template<typename Layout>
auto numNitrogens_Kitimar()
{
    std::size_t n = 0;
    benchmark_Kitimar<Layout>([&n] (auto &mol) {
        for (auto atom : get_atoms(mol))
            if (get_element(mol, atom) == 7)
                ++n;
    });
    std::cout << "# nitrogens: " << n << std::endl;
}



template<typename Layout, ctll::fixed_string SMARTS>
auto singleMatch_Kitimar()
{
    std::size_t n = 0;
    SingleIsomorphism<SMARTS> iso;
    benchmark_Kitimar<Layout>([&n, &iso] (auto &mol) {
        if (iso.match(mol))
            ++n;
    });
    std::cout << "# macthes: " << n << std::endl;
    return n;
}








#define STRINGIFY(s) #s
#define BENCHMARK(function) benchmark(STRINGIFY(function), &function)


int main()
{
    serialize_chembl<StructMoleculeIncident>();
    //serialize_chembl<StructMolecule, StructMoleculeIncident>();

    //BENCHMARK(numAtoms_OpenBabel_Sdf);                          // 71233387     650 seconds ~ 11 minutes
    //BENCHMARK(numAtoms_OpenBabel_Smi);                          // 70820782     139 seconds ~ 2.3 minutes
    //BENCHMARK(numAtoms_Kitimar<StructMolecule>);                // 71233387     0.344 seconds
    //BENCHMARK(numNitrogens_Kitimar<StructMolecule>);            // 8141624      0.479 seconds
    //BENCHMARK(numRing6_Kitimar<StructMolecule>);                // 18-19 seconds



    auto smarts = Smarts<"**1*2**12">{};
    auto bonds = smarts.bonds;
    auto dfsBonds = getDfsBonds(smarts);

    auto dfsAtoms = getDfsAtoms(smarts);
    //std::cout << bonds << std::endl;
    auto adj = adjacencyList<smarts.numAtoms, smarts.numBonds>(smarts.bonds);
    //std::cout << adj << std::endl;
    //std::cout << dfsBonds << std::endl;

    auto cycleMembership = getCycleMembership(smarts);

    //std::cout << cycleMembership << std::endl;

    //std::cout << dfsBonds << std::endl;



    //std::cout << dfsAtoms << std::endl;






    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "*1*2*1*2">));        //



    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "*1~*~*~*~*~*~1">));        // 3.7 seconds
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "c1ccccc1-c2ccccc2">));        // 2.54 seconds
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "c1ccccc1CCCCc1ccccc1">));        // 2.79 seconds
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "c1cc2c(cc1)cccc2">));        //
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "Clc1ccccc1">));        //
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "BrCCCCCCCCC">));        //
    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "[nH]1ccc2c1cccc2">));        //
    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "[nH]">));        //
    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "c1cc(=O)cc[nH]1">));        //
    BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "O1CCOC12CCNCC2">));        //




    //BENCHMARK((singleMatch_Kitimar<StructMoleculeIncident, "">));        //



}
*/
