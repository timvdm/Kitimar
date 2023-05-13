#include <Kitimar/Serialize/Serialize.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Util/Test.hpp>
#include <Kitimar/Util/Util.hpp>

#include <openbabel/parsmart.h>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;



template<ctll::fixed_string SMARTS>
auto countOBMatches()
{
    OpenBabel::OBSmartsPattern smarts;
    smarts.Init(Util::toString(SMARTS));

    auto n = 0;
    readMolecules(chembl_smi_filename(), [&] (auto &mol) {
        if (smarts.Match(mol))
            ++n;
        return true;
    });

    return n;
}

template<ctll::fixed_string SMARTS, typename Layout>
auto countCTMatchesFileStream()
{
    SingleIsomorphism<SMARTS> smarts{};


    std::ifstream ifs(chembl_serialized_filename(Layout{}), std::ios_base::binary | std::ios_base::in);
    LayoutSize::Type numMolecules;
    ifs.read(reinterpret_cast<char*>(&numMolecules), LayoutSize::size());

    auto n = 0;
    std::vector<std::byte> data;
    for (auto i = 0UL; i < numMolecules; ++i) {
        auto [ok, mol] = deserialize<Layout>(ifs, data);
        if (!ok)
            break;
        if (smarts.match(mol))
            ++n;
    }

    return n;
}

template<ctll::fixed_string SMARTS, typename Layout>
auto countCTMatchesMemoryMapped()
{
    SingleIsomorphism<SMARTS> smarts{};



    auto source = MemMapSource{chembl_serialized_filename(Layout{})};
    auto numMolecules = *reinterpret_cast<const LayoutSize::Type*>(source.begin());

    std::cout << "# molecules: " << numMolecules << std::endl;

    auto n = 0;
    std::size_t offset = LayoutSize::size();
    for (auto i = 0UL; i < numMolecules; ++i) {
        auto mol = deserialize<Layout>(source, offset);        
        if (smarts.match(mol))
            ++n;
    }

    return n;
}


template<ctll::fixed_string SMARTS>
auto benchmark()
{
    //auto [obCount, obElapsed] = std::make_pair(0, 0);
    auto [obCount, obElapsed] = benchmark(&countOBMatches<SMARTS>);
    auto [ifsCount, ifsElapsed] = benchmark(&countCTMatchesFileStream<SMARTS, StructMoleculeIncident>);
    auto [mmapCount, mmapElapsed] = benchmark(&countCTMatchesMemoryMapped<SMARTS, StructMoleculeIncident>);


    std::cout << Util::toString(SMARTS) << '\n';
    std::cout << "    OpenBabel       " << obCount << "    " << obElapsed << '\n';
    std::cout << "    FileStream      " << ifsCount << "    " << ifsElapsed << '\n';
    std::cout << "    MemoryMapped    " << mmapCount << "    " << mmapElapsed << '\n';
}



void test_mio()
{
    auto filename = "test_mio.bin";

    std::ofstream ofs{filename, std::ios_base::binary | std::ios_base::out};
    for (int i = 0; i < 10; ++i)
        ofs.write(reinterpret_cast<char*>(&i), sizeof(int));
    ofs.close();


    MemMapSource source{filename};

    for (int i = 0; i < 10; ++i)
        std::cout << *reinterpret_cast<const int*>(source.begin() + i * sizeof(i)) << std::endl;

    //auto *data = reinterpret_cast<const int*>(source.begin());
    //for (int i = 0; i < 10; ++i)
    //    std::cout << *(data + i * sizeof(i)) << std::endl;

}



int main()
{
    //serialize_chembl<StructMoleculeIncident>();

    //test_mio();
    benchmark<"c1ccccc1-c2ccccc2">();
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
