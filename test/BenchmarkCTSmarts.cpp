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

First iterative implementation
==============================

    git commit: 57c02ee2058d2ee47b421b3a7f38da38915f2ea2

    Callback + Recursive (A)
    ------------------------

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
    -------------------------

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
    ------------------------

    * Overall slower vs A -> optimize access to bond info w/o using with_n

        ------------------------------------------------------------------------         A
        Benchmark                              Time             CPU   Iterations         .
        ------------------------------------------------------------------------         . 28 %
        BM_count_CC                          515 ns          515 ns      1349903         |
        BM_count_xxx                       11230 ns        11230 ns        61668 ---+ <--+
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
    -------------------------

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


Optimize getQueryBondInfo, matchAtom & matchBond
================================================

    git commit: 7f0e5e96dfb90fed74074e17ef7261da2832f6f1

    Callback + Iterative (E)
    ------------------------

        ------------------------------------------------------------------------    A       C
        Benchmark                              Time             CPU   Iterations    .       .
        ------------------------------------------------------------------------    . 43 %  . 153 %
        BM_count_CC                          520 ns          520 ns      1349761    |       |
        BM_count_xxx                        7432 ns         7431 ns        94431 <--+-------+
        BM_count_xxxx                      13097 ns        13097 ns        53541
        BM_count_xxxxx                     19807 ns        19806 ns        35189
        BM_count_c1ccccc1                   1621 ns         1621 ns       426993
        BM_count_c1ccccc1CCN                1584 ns         1584 ns       439589
        BM_count_CCCCCCCCCCCCCCCCCCCC      11536 ns        11536 ns        60775
        BM_count_c1ccccc1CCCc1ccccc1        1638 ns         1638 ns       426888
        BM_count_AAAAAAAAAAAAAAAAAAA      203106 ns       203094 ns         3429
        ------------------------------------------------------------------------
        BM_all_CC                           6260 ns         6260 ns       111388
        BM_all_xxx                         14385 ns        14384 ns        48662
        BM_all_xxxx                        21285 ns        21284 ns        32724
        BM_all_xxxxx                       28928 ns        28926 ns        24164
        BM_all_c1ccccc1                     1583 ns         1583 ns       441394
        BM_all_c1ccccc1CCN                  1687 ns         1686 ns       417778
        BM_all_CCCCCCCCCCCCCCCCCCCC        11751 ns        11750 ns        59451
        BM_all_c1ccccc1CCCc1ccccc1          1684 ns         1684 ns       414877
        BM_all_AAAAAAAAAAAAAAAAAAA        223613 ns       223601 ns         3103


Store range/iterators in BondIters
==================================

    git commit: 940bb1d81dc28379cbd8c6db71826dbc4878ecda

    Callback + Iterative (F)
    ------------------------

        ------------------------------------------------------------------------
        Benchmark                              Time             CPU   Iterations
        ------------------------------------------------------------------------
        BM_count_CC                          509 ns          509 ns      1337878
        BM_count_xxx                        7455 ns         7455 ns        93307
        BM_count_xxxx                      13138 ns        13137 ns        53270
        BM_count_xxxxx                     19719 ns        19717 ns        35573
        BM_count_c1ccccc1                   1675 ns         1675 ns       418384
        BM_count_c1ccccc1CCN                1647 ns         1647 ns       424455
        BM_count_CCCCCCCCCCCCCCCCCCCC      11533 ns        11532 ns        60669
        BM_count_c1ccccc1CCCc1ccccc1        1672 ns         1672 ns       416413
        BM_count_AAAAAAAAAAAAAAAAAAA      210795 ns       210780 ns         3317
        ------------------------------------------------------------------------
        BM_all_CC                           6072 ns         6072 ns       113878
        BM_all_xxx                         14722 ns        14720 ns        47494
        BM_all_xxxx                        20893 ns        20892 ns        33548
        BM_all_xxxxx                       28993 ns        28991 ns        24021
        BM_all_c1ccccc1                     1666 ns         1665 ns       420468
        BM_all_c1ccccc1CCN                  1649 ns         1649 ns       424377
        BM_all_CCCCCCCCCCCCCCCCCCCC        11657 ns        11656 ns        60012
        BM_all_c1ccccc1CCCc1ccccc1          1713 ns         1713 ns       409208
        BM_all_AAAAAAAAAAAAAAAAAAA        218055 ns       218039 ns         3199


Use std::array to store IsmorphismMap
=====================================

    git commit: 72303176cf8625998488c0ba6fab286ca9d7d5ac

    Callback + Recursive (G)
    ------------------------

        ------------------------------------------------------------------------
        Benchmark                              Time             CPU   Iterations
        ------------------------------------------------------------------------
        BM_count_CC                          520 ns          520 ns      1336128
        BM_count_xxx                        3082 ns         3082 ns       227438 ---+
        BM_count_xxxx                       5695 ns         5694 ns       119936    |
        BM_count_xxxxx                      8240 ns         8239 ns        84883    |
        BM_count_c1ccccc1                    285 ns          285 ns      2455715    |
        BM_count_c1ccccc1CCN                 286 ns          286 ns      2454315    |
        BM_count_CCCCCCCCCCCCCCCCCCCC       3375 ns         3374 ns       207754    | 89 % (was 31% for A)
        BM_count_c1ccccc1CCCc1ccccc1         268 ns          268 ns      2603496    |
        BM_count_AAAAAAAAAAAAAAAAAAA       65181 ns        65176 ns        10611    |
        ------------------------------------------------------------------------    |
        BM_all_CC                           1458 ns         1458 ns       479019    |
        BM_all_xxx                          3474 ns         3473 ns       203396 <--+
        BM_all_xxxx                         6334 ns         6333 ns       110512
        BM_all_xxxxx                        8787 ns         8787 ns        79166
        BM_all_c1ccccc1                      274 ns          274 ns      2561168
        BM_all_c1ccccc1CCN                   277 ns          277 ns      2523717
        BM_all_CCCCCCCCCCCCCCCCCCCC         3437 ns         3437 ns       204048
        BM_all_c1ccccc1CCCc1ccccc1           286 ns          286 ns      2463427
        BM_all_AAAAAAAAAAAAAAAAAAA         69791 ns        69786 ns         9985


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
