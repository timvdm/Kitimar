#include <benchmark/benchmark.h>

#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>

#include <cstring>

namespace ctse = Kitimar::CTSmarts;

auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");

#define BM_COUNT(SMARTS, name) \
    void BM_count_##name(benchmark::State& state) { \
        for (auto _ : state) { \
            auto n = ctse::count<SMARTS>(mol, ctse::All); \
            benchmark::DoNotOptimize(n); \
            benchmark::ClobberMemory(); \
        } \
    } \
    BENCHMARK(BM_count_##name);

template<ctll::fixed_string SMARTS>
auto count_multi(auto &mol)
{
    int n = 0;
    for (const auto &map : ctse::multi<SMARTS>(mol, ctse::All))
        ++n;
    return n;
}

#define BM_ALL(SMARTS, name) \
    void BM_all_##name(benchmark::State& state) { \
        for (auto _ : state) { \
            auto maps = ctse::multi<SMARTS>(mol, ctse::All); \
            auto n = std::ranges::distance(maps); \
            benchmark::DoNotOptimize(n); \
            benchmark::ClobberMemory(); \
        } \
    } \
    BENCHMARK(BM_all_##name);

/*

Callback + Recursive (A)
========================

* Optimize multi for single bond case (BM_count_CC vs BM_all_CC)
* BM_all slow when there are many mappings (BM_count_xxx vs BM_all_xxx) -> avoid constructing std::vector<std::vector<int>>

    ------------------------------------------------------------------------
    Benchmark                              Time             CPU   Iterations
    ------------------------------------------------------------------------
    BM_count_CC                          505 ns          505 ns      1362260
    BM_count_xxx                        3151 ns         3151 ns       220864 ---+
    BM_count_xxxx                       5483 ns         5483 ns       127496    |
    BM_count_xxxxx                      8331 ns         8330 ns        84118    |
    BM_count_c1ccccc1                    309 ns          309 ns      2239155    |
    BM_count_c1ccccc1CCN                 299 ns          299 ns      2338269    |
    BM_count_CCCCCCCCCCCCCCCCCCCC       3405 ns         3405 ns       205118    | 31 %
    BM_count_c1ccccc1CCCc1ccccc1         284 ns          284 ns      2467194    |
    BM_count_AAAAAAAAAAAAAAAAAAA       65952 ns        65949 ns        10446    |
    ------------------------------------------------------------------------    |
    BM_all_CC                           4242 ns         4242 ns       166175    |
    BM_all_xxx                         10108 ns        10107 ns        69461 <--+
    BM_all_xxxx                        12752 ns        12751 ns        54528
    BM_all_xxxxx                       17252 ns        17252 ns        40653
    BM_all_c1ccccc1                      286 ns          286 ns      2440548
    BM_all_c1ccccc1CCN                   264 ns          264 ns      2743077
    BM_all_CCCCCCCCCCCCCCCCCCCC         3561 ns         3561 ns       193032
    BM_all_c1ccccc1CCCc1ccccc1           277 ns          277 ns      2555265
    BM_all_AAAAAAAAAAAAAAAAAAA         83651 ns        83647 ns         8402


Coroutine + Recursive (B)
=========================

* Overall slower vs A -> reduce # coroutines (stackless coroutine -> dynamically allocations)
* BM_count_xxx vs BM_all_xxx -> same (avoids std::vector<std::vector<int>>)

    ------------------------------------------------------------------------
    Benchmark                              Time             CPU   Iterations
    ------------------------------------------------------------------------
    BM_count_CC                          575 ns          575 ns      1212800
    BM_count_xxx                       15856 ns        15855 ns        43474 ---+
    BM_count_xxxx                      24499 ns        24498 ns        28308    |
    BM_count_xxxxx                     33018 ns        33017 ns        21206    |
    BM_count_c1ccccc1                    473 ns          473 ns      1483195    |
    BM_count_c1ccccc1CCN                 410 ns          410 ns      1706955    |
    BM_count_CCCCCCCCCCCCCCCCCCCC       9134 ns         9133 ns        76511    | 101 %
    BM_count_c1ccccc1CCCc1ccccc1         482 ns          482 ns      1456031    |
    BM_count_AAAAAAAAAAAAAAAAAAA      310949 ns       310935 ns         2249    |
    ------------------------------------------------------------------------    |
    BM_all_CC                           5444 ns         5444 ns       128007    |
    BM_all_xxx                         15910 ns        15909 ns        43989 <--+
    BM_all_xxxx                        24377 ns        24376 ns        28699
    BM_all_xxxxx                       33161 ns        33159 ns        21107
    BM_all_c1ccccc1                      510 ns          510 ns      1372926
    BM_all_c1ccccc1CCN                   435 ns          435 ns      1597545
    BM_all_CCCCCCCCCCCCCCCCCCCC         9107 ns         9107 ns        76726
    BM_all_c1ccccc1CCCc1ccccc1           501 ns          501 ns      1000000
    BM_all_AAAAAAAAAAAAAAAAAAA        316334 ns       316317 ns         2213

Callback + Iterative (C)
========================

* Overall slower vs A -> optimize access to bond info w/o using with_n

    ------------------------------------------------------------------------
    Benchmark                              Time             CPU   Iterations
    ------------------------------------------------------------------------
    BM_count_CC                          515 ns          515 ns      1349903
    BM_count_xxx                       11230 ns        11230 ns        61668 ---+
    BM_count_xxxx                      20463 ns        20462 ns        33853    |
    BM_count_xxxxx                     29033 ns        29032 ns        24172    |
    BM_count_c1ccccc1                   1637 ns         1637 ns       419762    |
    BM_count_c1ccccc1CCN                1685 ns         1685 ns       410006    |
    BM_count_CCCCCCCCCCCCCCCCCCCC      12677 ns        12676 ns        55408    | 62 %
    BM_count_c1ccccc1CCCc1ccccc1        1678 ns         1678 ns       416996    |
    BM_count_AAAAAAAAAAAAAAAAAAA      236101 ns       236090 ns         2961    |
    ------------------------------------------------------------------------    |
    BM_all_CC                           6315 ns         6315 ns       110905    |
    BM_all_xxx                         18210 ns        18209 ns        38386 <--+
    BM_all_xxxx                        27456 ns        27455 ns        25544
    BM_all_xxxxx                       38552 ns        38550 ns        18188
    BM_all_c1ccccc1                     1670 ns         1670 ns       418305
    BM_all_c1ccccc1CCN                  1668 ns         1668 ns       419563
    BM_all_CCCCCCCCCCCCCCCCCCCC        12937 ns        12936 ns        54038
    BM_all_c1ccccc1CCCc1ccccc1          1697 ns         1697 ns       412294
    BM_all_AAAAAAAAAAAAAAAAAAA        254367 ns       254355 ns         2752

Coroutine + Iterative (D)
=========================

* Same vs C, Slower vs A -> optimize iterative (see C)

    ------------------------------------------------------------------------
    Benchmark                              Time             CPU   Iterations
    ------------------------------------------------------------------------
    BM_count_CC                          539 ns          539 ns      1300290
    BM_count_xxx                       13416 ns        13416 ns        49996 ---+
    BM_count_xxxx                      22336 ns        22335 ns        31294    |
    BM_count_xxxxx                     32375 ns        32373 ns        21652    |
    BM_count_c1ccccc1                   1812 ns         1812 ns       384197    |
    BM_count_c1ccccc1CCN                1787 ns         1787 ns       396095    |
    BM_count_CCCCCCCCCCCCCCCCCCCC      14738 ns        14737 ns        47379    | 104 %
    BM_count_c1ccccc1CCCc1ccccc1        1858 ns         1858 ns       374683    |
    BM_count_AAAAAAAAAAAAAAAAAAA      277141 ns       277123 ns         2527    |
    ------------------------------------------------------------------------    |
    BM_all_CC                           4377 ns         4377 ns       159758    |
    BM_all_xxx                         13363 ns        13362 ns        52308 <--+
    BM_all_xxxx                        22398 ns        22397 ns        31163
    BM_all_xxxxx                       32177 ns        32174 ns        21730
    BM_all_c1ccccc1                     1838 ns         1838 ns       380523
    BM_all_c1ccccc1CCN                  1815 ns         1815 ns       387326
    BM_all_CCCCCCCCCCCCCCCCCCCC        14777 ns        14776 ns        47126
    BM_all_c1ccccc1CCCc1ccccc1          1937 ns         1937 ns       360951
    BM_all_AAAAAAAAAAAAAAAAAAA        279467 ns       279452 ns         2538




*/



//BM_COUNT("C", C);
BM_COUNT("CC", CC);
BM_COUNT("***", xxx);
BM_COUNT("****", xxxx);
BM_COUNT("*****", xxxxx);
BM_COUNT("c1ccccc1", c1ccccc1);
BM_COUNT("c1ccccc1CCN", c1ccccc1CCN);
BM_COUNT("CCCCCCCCCCCCCCCCCCC", CCCCCCCCCCCCCCCCCCCC);
BM_COUNT("c1ccccc1CCCc1ccccc1", c1ccccc1CCCc1ccccc1);
BM_COUNT("AAAAAAAAAAAAAAAAAAA", AAAAAAAAAAAAAAAAAAA);

//BM_ALL("C", C);
BM_ALL("CC", CC);
BM_ALL("***", xxx);
BM_ALL("****", xxxx);
BM_ALL("*****", xxxxx);
BM_ALL("c1ccccc1", c1ccccc1);
BM_ALL("c1ccccc1CCN", c1ccccc1CCN);
BM_ALL("CCCCCCCCCCCCCCCCCCC", CCCCCCCCCCCCCCCCCCCC);
BM_ALL("c1ccccc1CCCc1ccccc1", c1ccccc1CCCc1ccccc1);
BM_ALL("AAAAAAAAAAAAAAAAAAA", AAAAAAAAAAAAAAAAAAA);


// Run the benchmark
BENCHMARK_MAIN();
