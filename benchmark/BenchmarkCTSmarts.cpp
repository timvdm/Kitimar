#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
//#include <Kitimar/RDKit/RDKit.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <cstring>

#include <Kitimar/Util/Util.hpp>

#include "TestData.hpp"

#include "Benchmarks.hpp"

using namespace Kitimar;

template<ctll::fixed_string SMARTS>
struct FixedString
{
    static constexpr inline auto smarts = SMARTS;
};

using ChainBenchmarks = std::tuple<
    //FixedString<"C">,
    FixedString<"CC">,
    FixedString<"***">,
    FixedString<"CCC">,
    FixedString<"CCO">,
    FixedString<"CNO">,
    FixedString<"******Cl">,
    FixedString<"Cl******">
    /*
    FixedString<"*">,
    FixedString<"**">,

    FixedString<"***">,
    FixedString<"****">,
    FixedString<"*****">,
    FixedString<"******">,
    FixedString<"*******">,
    FixedString<"********">,
    FixedString<"********">,
    FixedString<"*********">,
    FixedString<"**********">,
    FixedString<"***********">
    */
>;

enum class BenchmarkArgs {
    Mol,
    Atom,
    Bond
};

#define BENCHMARK_STR(function) #function

#define BENCHMARK_FUNCTION(function) \
    struct benchmark_##function { \
        static constexpr inline auto args = BenchmarkArgs::Mol; \
        static const char* name() noexcept { return BENCHMARK_STR(function); } \
        template<ctll::fixed_string SMARTS, typename Config> \
        static constexpr auto call(auto &mol) \
        { return ctse::function<SMARTS, Config>(mol); } \
    }

#define BENCHMARK_FUNCTION_ARGS(function, Args) \
    struct benchmark_##function { \
        static constexpr inline auto args = BenchmarkArgs::Args; \
        static const char* name() noexcept { return BENCHMARK_STR(function); } \
        template<ctll::fixed_string SMARTS, typename Config> \
        static constexpr auto call(auto &mol, const auto &x) \
        { return ctse::function<SMARTS, Config>(mol, x); } \
    }

#define BENCHMARK_FUNCTION_ATOM(function) BENCHMARK_FUNCTION_ARGS(function, Atom)
#define BENCHMARK_FUNCTION_BOND(function) BENCHMARK_FUNCTION_ARGS(function, Bond)


BENCHMARK_FUNCTION(match);
BENCHMARK_FUNCTION_ATOM(match_atom);
BENCHMARK_FUNCTION_BOND(match_bond);

BENCHMARK_FUNCTION(count_unique);
BENCHMARK_FUNCTION(count_all);
BENCHMARK_FUNCTION_ATOM(count_atom_unique);
BENCHMARK_FUNCTION_ATOM(count_atom_all);
BENCHMARK_FUNCTION_BOND(count_bond_unique);
BENCHMARK_FUNCTION_BOND(count_bond_all);

BENCHMARK_FUNCTION(map);
BENCHMARK_FUNCTION_ATOM(map_atom);
BENCHMARK_FUNCTION_BOND(map_bond);

BENCHMARK_FUNCTION(maps_unique);
BENCHMARK_FUNCTION(maps_all);
BENCHMARK_FUNCTION_ATOM(maps_atom_unique);
BENCHMARK_FUNCTION_ATOM(maps_atom_all);
BENCHMARK_FUNCTION_BOND(maps_bond_unique);
BENCHMARK_FUNCTION_BOND(maps_bond_all);

BENCHMARK_FUNCTION(capture);
BENCHMARK_FUNCTION_ATOM(capture_atom);
BENCHMARK_FUNCTION_BOND(capture_bond);

BENCHMARK_FUNCTION(captures_unique);
BENCHMARK_FUNCTION(captures_all);
BENCHMARK_FUNCTION_ATOM(captures_atom_unique);
BENCHMARK_FUNCTION_ATOM(captures_atom_all);
BENCHMARK_FUNCTION_BOND(captures_bond_unique);
BENCHMARK_FUNCTION_BOND(captures_bond_all);


using SingleAtomFunctions = std::tuple<
    benchmark_match,
    benchmark_match_atom,

    benchmark_count_unique,
    benchmark_count_all,
    benchmark_count_atom_unique,
    benchmark_count_atom_all,

    benchmark_map,

    benchmark_maps_unique,
    benchmark_maps_all,

    benchmark_capture,

    benchmark_captures_unique,
    benchmark_captures_all
>;

using AllFunctions = std::tuple<
    benchmark_match,
    benchmark_match_atom,
    benchmark_match_bond,

    benchmark_count_unique,
    benchmark_count_all,
    benchmark_count_atom_unique,
    benchmark_count_atom_all,
    benchmark_count_bond_unique,
    benchmark_count_bond_all,

    benchmark_map,
    benchmark_map_atom,
    benchmark_map_bond,

    benchmark_maps_unique,
    benchmark_maps_all,
    benchmark_maps_atom_unique,
    benchmark_maps_atom_all,
    benchmark_maps_bond_unique,
    benchmark_maps_bond_all,

    benchmark_capture,
    benchmark_capture_atom,
    benchmark_capture_bond,

    benchmark_captures_unique,
    benchmark_captures_all,
    benchmark_captures_atom_unique,
    benchmark_captures_atom_all,
    benchmark_captures_bond_unique,
    benchmark_captures_bond_all
>;



using Configs = std::tuple<ctse::DefaultConfig, ctse::NoOptimizeConfig>;

std::ostream& operator<<(std::ostream &os, const ctse::DefaultConfig) { os << "DefaultConfig"; return os; }
std::ostream& operator<<(std::ostream &os, const ctse::NoOptimizeConfig) { os << "NoOptimizeConfig"; return os; }
//std::ostream& operator<<(std::ostream &os, const ctse::NoSpecializeConfig) { os << "NoSpecializeConfig"; return os; }


template<typename Functions, typename SMARTS, typename Configs, int I, int J, int K>
struct BenchmarkCase
{
    static constexpr inline auto function = std::get<I>(Functions{});
    static constexpr inline auto smarts = std::get<J>(SMARTS{}).smarts;
    static constexpr inline auto config = std::get<K>(Configs{});

    using Config = std::remove_cvref_t<decltype(config)>;

    static void run(auto &mol, const auto &atom, const auto &bond)
    {
        std::stringstream ss;
        ss << "ctse::" << function.name()
           << "<\"" << Util::toString(smarts) << "\", "
           << config << ">(mol"
           << (function.args == BenchmarkArgs::Atom ? ", atom" : "")
           << (function.args == BenchmarkArgs::Bond ? ", atom" : "")
           << ")";

        BENCHMARK(ss.str()) {
            if constexpr (function.args == BenchmarkArgs::Mol)
                return function.template call<smarts, Config>(mol);
            else if constexpr (function.args == BenchmarkArgs::Atom)
                return function.template call<smarts, Config>(mol, atom);
            else if constexpr (function.args == BenchmarkArgs::Bond)
                return function.template call<smarts, Config>(mol, bond);
        };
    }
};

template<typename Functions, typename SMARTS, typename Configs, int I = 0, int J = 0, int K = 0>
constexpr auto create_benchmarks()
{
    if constexpr (I == std::tuple_size_v<Functions>)
        return std::tuple<>{};
    else if constexpr (J == std::tuple_size_v<SMARTS>)
        return create_benchmarks<Functions, SMARTS, Configs, I + 1, 0, 0>();
    else if constexpr (K == std::tuple_size_v<Configs>)
        return create_benchmarks<Functions, SMARTS, Configs, I, J + 1, 0>();
    else {
        auto tail = create_benchmarks<Functions, SMARTS, Configs, I, J, K + 1>();
        return std::tuple_cat(std::make_tuple(BenchmarkCase<Functions, SMARTS, Configs, I, J, K>{}), tail);
    }
}

template<typename T>
struct identify_type;

TEST_CASE("Debug")
{





}

//using Benchmarks = decltype(create_benchmarks());

/*
using Benchmarks = decltype(create_benchmarks<std::tuple<benchmark_maps_all>, ChainBenchmarks, Configs>());

TEMPLATE_LIST_TEST_CASE("maps_all", "", Benchmarks)
{
    using Benchmark = TestType;
    auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");
    auto atom = get_atom(mol, 0);
    auto bond = get_bond(mol, 0);
    Benchmark::run(mol, atom, bond);
}
*/

#define BENCHMARK_CASE(name, Functions, SMARTS, Configs) \
using Benchmarks_##name = decltype(create_benchmarks<Functions, SMARTS, Configs>()); \
TEMPLATE_LIST_TEST_CASE(BENCHMARK_STR(name), "", Benchmarks_##name) \
{ \
    using Benchmark = TestType; \
    auto mol = Kitimar::Toolkit::readSmiles<Toolkit::openbabel>("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1"); \
    auto atom = get_atom(mol, 0); \
    auto bond = get_bond(mol, 0); \
    Benchmark::run(mol, atom, bond); \
}


#define BENCHMARK_FUNCTION_CASE(function) \
    using Functions_##function = std::tuple<benchmark_##function>; \
    BENCHMARK_CASE(function, Functions_##function, ChainBenchmarks, Configs)

BENCHMARK_FUNCTION_CASE(match)
BENCHMARK_FUNCTION_CASE(match_atom)
BENCHMARK_FUNCTION_CASE(match_bond)

BENCHMARK_FUNCTION_CASE(count_unique)
BENCHMARK_FUNCTION_CASE(count_all)
BENCHMARK_FUNCTION_CASE(count_atom_unique)
BENCHMARK_FUNCTION_CASE(count_atom_all)
BENCHMARK_FUNCTION_CASE(count_bond_unique)
BENCHMARK_FUNCTION_CASE(count_bond_all)

BENCHMARK_FUNCTION_CASE(map)
BENCHMARK_FUNCTION_CASE(map_atom)
BENCHMARK_FUNCTION_CASE(map_bond)

BENCHMARK_FUNCTION_CASE(capture)
BENCHMARK_FUNCTION_CASE(capture_atom)
BENCHMARK_FUNCTION_CASE(capture_bond)

using SMARTS_single_atom = std::tuple<FixedString<"C">>;
BENCHMARK_CASE(single_atom, SingleAtomFunctions, SMARTS_single_atom, Configs)

using SMARTS_single_bond = std::tuple<FixedString<"CC">>;
BENCHMARK_CASE(single_bond, AllFunctions, SMARTS_single_bond, Configs)


/*
 * Variables:
 *
 * - function
 * - SMARTS
 * - Specialize
 * - Optimize
 *
 *
 *
 *
 *
 *
 *
 */





/*

TEST_CASE("optimized_match_atom")
{
    auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");
    auto atom = get_atom(mol, 0);

    BENCHMARK("single atom cpp")  { return match_atom_single_atom_cpp(mol, atom); };
    BENCHMARK("single atom ctse") { return match_atom_single_atom_ctse(mol, atom); };

    BENCHMARK("single bond cpp")  { return match_atom_single_bond_cpp(mol, atom); };
    BENCHMARK("single bond ctse") { return match_atom_single_bond_ctse(mol, atom); };

    BENCHMARK("chain 3 cpp")  { return match_atom_chain_3_cpp(mol, atom); };
    BENCHMARK("chain 3 ctse") { return match_atom_chain_3_ctse(mol, atom); };

    BENCHMARK("chain 4 cpp")  { return match_atom_chain_4_cpp(mol, atom); };
    BENCHMARK("chain 4 ctse") { return match_atom_chain_4_ctse(mol, atom); };
}


TEST_CASE("optimized_match")
{
    auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");

    BENCHMARK("single atom cpp")  { return match_single_atom_cpp(mol); };
    BENCHMARK("single atom ctse") { return match_single_atom_ctse(mol); };

    BENCHMARK("single bond cpp")  { return match_single_bond_cpp(mol); };
    BENCHMARK("single bond ctse") { return match_single_bond_ctse(mol); };

    BENCHMARK("chain 3 cpp")  { return match_chain_3_cpp(mol); };
    BENCHMARK("chain 3 ctse") { return match_chain_3_ctse(mol); };

    BENCHMARK("chain n 3 cpp")  { return match_chain_n_cpp<3>(mol); };
    BENCHMARK("chain n 4 cpp")  { return match_chain_n_cpp<4>(mol); };
    BENCHMARK("chain n 5 cpp")  { return match_chain_n_cpp<5>(mol); };
    BENCHMARK("chain n 6 cpp")  { return match_chain_n_cpp<6>(mol); };
    BENCHMARK("chain n 7 cpp")  { return match_chain_n_cpp<7>(mol); };
}


TEST_CASE("match_any_chain")
{
    auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");

    BENCHMARK("chain 5")  { return ctse::match<"*****">(mol); };
    BENCHMARK("chain 6")  { return ctse::match<"******">(mol); };
    BENCHMARK("chain 7")  { return ctse::match<"*******">(mol); };
    BENCHMARK("chain 8")  { return ctse::match<"********">(mol); };
    BENCHMARK("chain 9")  { return ctse::match<"*********">(mol); };
    BENCHMARK("chain 10")  { return ctse::match<"**********">(mol); };
}

TEST_CASE("maps_all_any_chain")
{
    auto mol = Kitimar::readSmilesOpenBabel("N1(C(=O)[C@H]2N(C(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CO)CCC2)[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N[C@H](C(=O)O)[C@H](CC)C)C)C)Cc2nc[nH]c2)CC(=O)N)CO)C)CCC(=O)N)CC(C)C)C)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCCN=C(N)N)CCC(=O)N)CC(C)C)CCC1");

    BENCHMARK("chain 5")  { return ctse::maps_all<"*****">(mol); };
    BENCHMARK("chain 6")  { return ctse::maps_all<"******">(mol); };
    BENCHMARK("chain 7")  { return ctse::maps_all<"*******">(mol); };
    BENCHMARK("chain 8")  { return ctse::maps_all<"********">(mol); };
    BENCHMARK("chain 9")  { return ctse::maps_all<"*********">(mol); };
    BENCHMARK("chain 10")  { return ctse::maps_all<"**********">(mol); };
}




*/



















