# Kitimar

Kitimar is a container project for CTSmarts providing utilities for development, testing and benchmarking.

## Table of concents

- [CTSmarts: C++ compile-time SMARTS expressions](#ctsmarts-c-compile-time-smarts-expressions)
    - [STATUS: Proof of concept](#status-proof-of-concept)
    - [Features](#features)
    - [Build instructions](#build-instructions)
    - [Examples](#examples)
    - [API](#api)
    - [Performance](#performance)
    - [Supported compilers](#supported-compilers)
    - [Planned features](#planned-features)
    - [History](#history)
- [Molecule: C++20 Molecule concept](#molecule-c20-molecule-concept)
- [Modules](#modules)
- [Thanks](#thanks)

## CTSmarts: C++ compile-time SMARTS expressions

`CTSmarts` provides fast compile-time [SMARTS expressions](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification) with support for matching, counting, mapping and capturing.

The SMARTS expressions are parsed at compile-time resulting in a minimal amount of assembly code and need for dynamically alocated data structures at run-time.

It is implemented as a C++20 header-only library and inspired by the [Compile time regular expressions](https://github.com/hanickadot/compile-time-regular-expressions) and uses it's `ctll` compile-time LL(1) parser.

### STATUS: Proof of concept

This is currently a working proof of concept. The aim was to see if this would possible in C++20 using current compilers.
Many features are still missing and/or untested.

#### Known limitations

- Multiple connected components
- Stereochemistry
- Invalid SMARTS
    - Compile errors need to be improved using `static_assert` and easy to understand messages.
    - Workaround: validate SMARTS expressions using other tools first
- Performance
    - gcc++11 has no `std::ranges::owning_view` &rarr; maps are copied before returning ([range-v3](https://github.com/ericniebler/range-v3)?)

### Features

- No parsing/optimization overhead
- Performance same or better than to handwritten code
- Capture atoms using SMARTS atom classes
- No need for global state &rarr; no threading issues
- No need for external SMARTS to C++ conversion tool and build system setup
- Compile time optimizations:
    - Add cyclicity check to cyclic SMARTS atoms and bonds
    - More planned...

See [planned features](#planned-features) below.

### Build instructions

Build without any external toolkit support:

```shell
~$ git clone https://github.com/timvdm/Kitimar.git
~$ cd Kitimar
~/Kitimar$ mkdir build
~/Kitimar$ cd build
~/Kitimar/build$ cmake -DCMAKE_BUILD_TYPE=Release ..
~/Kitimar/build$ make
~/Kitimar/build$ make test
```

Add `-DKitimar_WITH_OPENBABEL=ON -DOpenBabel3_DIR=/path/to/lib/cmake/openbabel3` to build with OpenBabel3 support
or with `-DKitimar_WITH_RDKIT=ON -DRDKit_DIR=/path/to/lib/cmake/rdkit` for RDKit support.


### Examples

The examples can be found in the `examples/CTSmarts` directory.

```shell
~/Kitimar/build$ ./examples/CTSmarts/Match
Matching SMARTS "C[O-]" in SMILES "CC(=O)[O-]": 1
~/Kitimar/build$ ./examples/CTSmarts/Mapping
...
~/Kitimar/build$ ./examples/CTSmarts/Capture
...
```

### API

The code below is a simplified version of the real API which makes it easier to read form a users perspective. Examples of how to use the API can be found in the `examples/CTSmarts` directory.

```c++
// Context

#include <Kitimar/CTSmarts/CTSmarts.hpp>

auto mol = ...; // Molecule that satifies the Molecule::Molecule concept.

// ctse::match -> bool

ctse::match(mol);
ctse::match_atom(mol, atom);
ctse::match_bond(mol, bond);
ctse::match(mol, atom/bond);

// ctse::count -> int

ctse::count(mol, type = Unique);
ctse::count_unique(mol);
ctse::count_all(mol);

ctse::count_atom(mol, atom, type = Unique);
ctse::count_atom_unique(mol, atom);
ctse::count_atom_all(mol, atom);

ctse::count_bond(mol, bond, type = Unique);
ctse::count_bond_unique_bond(mol, bond);
ctse::count_bond_all(mol, bond);

ctse::count(mol, atom/bond, type = Unique);
ctse::count_unique(mol, atom/bond);
ctse::count_all(mol, atom/bond);

// ctse::map -> std::tuple<bool, std::array<int, N>>

ctse::map(mol);
ctse::map_atom(mol, atom);
ctse::map_bond(mol, bond);
ctse::map(mol, atom/bond);

auto [found, map] = ctse::map<"C=O">(mol);
if (found) {
    // Use map...
}

// ctse::maps -> std::vector<std::array<int, N>>

ctse::maps(mol, type = Unique);
ctse::maps_unique(mol);
ctse::maps_all(mol);

ctse::maps_atom(mol, atom, type = Unqiue);
ctse::maps_atom_unique(mol, atom);
ctse::maps_atom_all(mol, atom);

ctse::maps_bond(mol, bond, type = Unique);
ctse::maps_bond_unique(mol, bond);
ctse::maps_bond_all(mol, bond);

ctse::maps(mol, atom/bond, type = Unique);
ctse::maps_unique(mol, atom/bond);
ctse::maps_all(mol, atom/bond);

// ctse::capture -> std::tuple<bool, Atom...>

cst::capture(mol);
cst::capture_atom(mol, atom);
cst::capture_bond(mol, bond);
cst::capture(mol, atom/bond);

// ctse::captures -> std::vector<std::tuple<Atom...>>

cst::captures(mol, type = Unique);
cst::captures_unique(mol);
cst::captures_all(mol);

cst::captures_atom(mol, atom, type = Unique);
cst::captures_atom_unique(mol, atom);
cst::captures_atom_all(mol, atom);

cst::captures_bond(mol, bond, type = Unique);
cst::captures_bond_unique(mol);
cst::captures_bond_all(mol);

cst::captures(mol, atom/bond, type = Unique);
cst::captures_unique(mol, atom/bond);
cst::captures_all(mol, atom/bond);

for (auto [C, N] : ctse::captures<"C-N">(mol)) {
    // Use atoms C and N
}
```

#### Capture

The `ctse::capture(s)` API is meant to be used with [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding) and [range-based for loops](https://en.cppreference.com/w/cpp/language/range-for).

```c++
// Capture all atoms
auto [match, C, O, N] = ctse::capture<"C(=O)N">(mol);
if (match) {
   // Use C, O, N atoms...
}
```

```c++
// Capture atoms specified using SMARTS atom classes
auto [match, O, N] = ctse::capture<"C(=[O:1])[N:2]">(mol);
if (match) {
   // Use O, N atoms...
}
```

```c++
// Capture all atoms for every unique mapping
for (auto [C, O] : ctse::captures<"C=O">(mol)) {
   // Use C, O atoms...
}
```

### Performance

#### Optimized cases

Simple cases are optimized to match the perormance of handwritten code.
In most of these cases the generated assembly is actually the same.

The optimized cases are currently:
- A single atom (e.g. `C`, `[O-]`, `[CD2,CD3;R]`)
- A single bond (e.g. `*=*`, `C=O`, `C-,=N`)
- ***Not yet implemented***
    - Chains (e.g. `CCCCCO`)
    - Star graphs (e.g. `CC=O`, `CC(=O)O`, `C(C)([O-])=O`)

##### Example of single atom SMARTS

The two functions generate the same assembly code. More complex examples can be tested on the
[compiler exlorer website](https://godbolt.org/z/cYYjnY4Me).

```c++
bool isCarbonDegree3_v1(auto &mol, auto atom)
{
    return get_element(mol, atom) == 6 && get_degree(mol, atom) == 3;
}

bool isCarbonDegree3_v2(auto &mol, auto atom)
{
    return ctse::atom<"[#6D3]">(mol, atom);
}
```

```assembly
isCarbonDegree3_v1
 xor    eax,eax
 cmp    BYTE PTR [rdx+0x7],0x6
 jne    4010f6 <main+0x16>
 cmp    BYTE PTR [rdx+0xa],0x3
 sete   al
 movzx  eax,al
 ret

isCarbonDegree3_v2
 xor    eax,eax
 cmp    BYTE PTR [rdx+0x7],0x6
 jne    4010f6 <main+0x16>
 cmp    BYTE PTR [rdx+0xa],0x3
 sete   al
 movzx  eax,al
 ret
```

#### General case

In the general case, a [subgraph-isomorphism search](https://en.wikipedia.org/wiki/Subgraph_isomorphism_problem) is needed.

A very simple subgraph-isomorphism algorithm is used which requires a **minimal amount of memory**.
At **compile-time**, the **matching order** of the SMARTS expression is **optimized**.
The result is a list of SMARTS bonds in the order of a depth-first-search tree. Molecules are matched against this tree using a backtracking.
To **prune** the search space, a simple check is used to ensure the degree of a candidate atom is greater or equal to the degree of the SMARTS atom.

More advanced algorithms with better pruning rules (e.g. [Ullmann](https://dl.acm.org/doi/abs/10.1145/321921.321925),
[VF2](https://ieeexplore.ieee.org/abstract/document/1323804), [VF2++](https://www.sciencedirect.com/science/article/pii/S0166218X18300829#b8)) require additional memory and overhead to apply these pruning rules. In the case of molecular graphs, these algorithms often perform worse because of this. See [The secrets of fast SMARTS matching](https://www.nextmovesoftware.com/talks/Mayfield_SecretsOfFastSmartsMatching_Sheffield_201906.pdf) for a more detailed explanation.

##### Memory usage

| Name                    | Type                              | Size (N)       | Allocation | Conditional
|-------------------------|-----------------------------------|----------------|------------|------------
| SMARTS atom degrees     | `std::array<uint8_t, N>`          | # SMARTS atoms | Static     | -
| Current mapping         | `std::array<int, N>`              | # SMARTS atoms | Static     | -
| Current mapped atom set | `std::vector<bool>`               | num_atoms(mol) | Dynamic    | -
| Unique atom set hashes  | `std::unordered_set<std::size_t>` | # mappings     | Dynamic    | ctse::Unique


### Supported compilers

CTSmarts requires a `C++20` compiler and associated standard template library.

Known **supported** compilers:

- GCC 12.1+ (12.0 untested)
- Clang 16.0+
- GCC 11.3+ (see [known limitations](#known-limitations) above)

Known **unsupported** compilers:

- Clang 15.0

### Planned features

- More compile-time optimizations
  - Match less common atoms first ( DONE )
  - Optimize atom and bond expressions ( WIP )
  - API for custom optimizers ( DONE )
- API for custom SMARTS extensions
- Improved error reporting
- Optimize special cases (match single atom/bond, central atom with bonds, ...) ( WIP )
    - No need for dynamic memory allocations
- ...

### History

Although I have only been working on this for about a month, the goal of implementing something like this is not new.

When I was working on [OpenBabel](https://github.com/openbabel/openbabel) 10 to 15 years ago, I often had to write code to find atoms and bonds that satisfy some criteria using if statements.
If the criteria were more complex, for loops would be needed to iterate over adjacent atoms and bonds. As soon as there are 2 for loops, you also need to keep track of visited atoms etc.
This repetitive code which is hard to read, not fun to write and distracts from the broader context.

The [C++ Standard Template Library (STL)](https://en.cppreference.com/w/cpp/standard_library) solved a similar problem
long ago using [containers](https://en.cppreference.com/w/cpp/container), [iterators](https://en.cppreference.com/w/cpp/iterator) and [algorithms](https://en.cppreference.com/w/cpp/algorithm).
A *"C++ Molecule Template Library (MTL)"* could have saved me a lot of time. Unfortunately this did not exist.
However, writing this hypothetical *"MTL"* in C++98 is no easy task. A replacement for all 3 would be needed.
[Helium](https://github.com/timvdm/Helium/) was my first attempt to implement this *"MTL"*.

The **container** was the easy part. While molecules are labeled graphs and no iterator pairs over a range of values,
a replacement for this was already available. The [Boost Graph Library (BGL)](boost.org/doc/libs/1_78_0/libs/graph/doc/index.html),
which is also used by [RDKit](https://www.rdkit.org/), solved this by defining a graph concept (in documentation) and a set of free functions.
Developers could implement these functions for their own custom graph data structure and use all of the algorithms provided by the BGL.
**Iterators** are still iterators, but these iterate over the atoms/bonds of a molecule or incident/adjacent to another atom. Dereferencing an iterator no longer returns a value but an atom or a bond.
These atoms/bonds have a vertex/edge label containing their properties (element, charge, bond order, ...) which can also be retrieved using a free function.

The **Molecule concept** was now defined. In the documentation at least, C++98 did not have concepts. Older compilers would also print errors hundreds of lines long when you tried to call std::sort on an std::list.
Template meta-programming which was required to implement this was also a time consuming task.
To make things worse, C++98 for loops using iterators were not exactly elegant code.
The developer using the *"MTL"* should be able to write compact, readable, performant code. All of this without being a template meta-programming wizard.

Nevertheless, I started working on what would become Helium to see how far I would get. Work on the new C++0x standard was in progress and seemed to contain many features that would make all of this easier.
The new [Clang](https://en.wikipedia.org/wiki/Clang) compiler was a major advancement regarding the compiler errors.
It produced colored error messages which were easier to read, less verbose and better than GCC at the time.
Progress was initially slow, the C++0x standard was delayed and became C++11.
Compiler support for the new standard would still have to catch on.
The first versions of Helium (2012-2014) still used C++98 with optional support C++11 `std:thread` to do similarity searches in parallel.
A year later the project was updated to require C++11 and use range based for loops to replace the macros using `BOOST_TYPEOF` and associated boilerplate code.

The **algorithms** found in Helium include a number of basic graph algorithms (e.g. breath-first-search, depth-first-search, cycles, Dijkstra's shortest path, connected components, canonicalization, ...)
and cheminformatics algorithms (aromatize, kekulize, subgraph isomorphisms using SMARTS expressions, stereochemistry, ...).
These are **high-level algorithms**, written using the same repetitive code I mentioned earlier.
These blocks of code are the low-level algorithms to query a molecule and could easily be expressed as SMARTS patterns.

Implementing these **low-level algorithms** efficiently with little to no overhead compared to the handwritten code was harder.
Parsing the **SMARTS expression at run-time** creates **overhead**.
A shared instance between threads which is only initialized once requires synchronisation between threads.
Furthermore, this results in a graph representation and abstract syntax trees (AST) for the atom and bond expressions existing in memory. This layer of indirection will never be as fast as the handwritten code.

For simple cases, Helium tried to provide predefined predicates which could be used and combined to do some of these tasks.
This never resulted in an elegant syntax or provided any real advantage. More recently, I also tried to use the ranges library which is now part of the C++20 standard for this.
Again, the results were never what I was looking for. SMARTS expressions are just a very compact way of representing these queries.
An equivalent way of expressing these **queries using valid C++** syntax will always be **more verbose** and requires a **new, unfamiliar syntax** to be learned.
Existing SMARTS expressions would also have to be converted to this syntax.

Another approach trying to combine the advantages of SMARTS with the performance of handwritten code are tools that **generate code from SMARTS**.
Examples of this include [NextMove Software's Patsy](https://www.nextmovesoftware.com/patsy.html) and the [SmartsCompiler](https://github.com/timvdm/SmartsCompiler) I started writing.
Patsy is definitely a major improvement over simple SMARTS implementations. The preprocessing also allows for optimizations which would otherwise create extra overhead.
However, to use tools like this, you have to specify the SMARTS patterns somewhere, setup a build system to convert the SMARTS to code and call this generated code.
This is great for doing things like atom typing, fingerprints, matching PAINS patterns, etc. Still, this does not solve the original problem described above.
I'm sure there are ways to get this working with macros, parsing C++ files and so on but this feels like a hack.

The goal of replacing this:

```c++
bool OBAtom::IsAromaticNOxide()
{
  if (GetAtomicNum() != OBElements::Nitrogen || !IsAromatic())
    return(false);

  OBAtom *atom;
  OBBondIterator i;

  for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
    if (atom->GetAtomicNum() == OBElements::Oxygen && !(*i)->IsInRing() && (*i)->GetBondOrder() == 2)
      return(true);

  return(false);
}
```

With the code below seemed impossible using valid C++11 syntax alone.

```c++
bool OBAtom::IsAromaticNOxide()
{
    return <magic>("n=[OR0]", this);
}
```

*The next few years life happened and I did not have time to work on this. During this time, the new C++17/20 standards were released.
For my current job I also use C++ and for some time now, I thought revisiting this idea would be a great way to play with the new C++ features I was not already using.*

I knew [**constexpr**](https://en.cppreference.com/w/cpp/language/constexpr) might be capable of this now, but it was not immediately clear how this would work.
There are still limitations that make it difficult to transfer the result of compile-time computations to run time.
The final piece I needed to put this altogether for this was the [Compile time regular expressions (ctre) library](https://github.com/hanickadot/compile-time-regular-expressions).
Using this as a example and basis I quickly had a working prototype.

## Molecule: C++20 Molecule concept

This module specifies the `Molecule` [C++20 concept](https://en.cppreference.com/w/cpp/language/constraints).
The concept requires a simple API consisting of free functions to be implemented. This allows molecule
data structures from different cheminformatics toolkits to be used.

```c++

// mol = Molecule object (e.g. OpenBabel::OBMol, RDKit::ROMol, ...)
// Atom atom, Bond bond =  (proxy) object, pointer, integer, ...
// index, size = int, std::size_t, ...

num_atoms(mol) -> std::integral;
num_bonds(mol) -> std::integral;

get_atom(mol, index) -> Atom;
get_bond(mol, index) -> Bond;

get_atoms(mol) -> std::ranges::input_range<Atom>;
get_bonds(mol) -> std::ranges::input_range<Bond>;

get_bonds(mol, atom) -> std::ranges::input_range<Bond>;
get_nbrs(mol, atom)  -> std::ranges::input_range<Atom>;

get_index(mol, atom)              -> std::integral;
get_degree(mol, atom)             -> std::integral;
get_element(mol, atom)            -> std::integral;
get_isotope(mol, atom)            -> std::integral;
get_charge(mol, atom)             -> std::signed_integral;
get_valence(mol, atom)            -> std::integral;
get_implicit_hydrogens(mol, atom) -> std::integral;
get_total_hydrogens(mol, atom)    -> std::integral;
is_ring_atom(mol, atom)           -> bool;
is_aromatic_atom(mol, atom)       -> bool;
is_in_ring_size(mol, atom, size)  -> bool;
get_ring_count(mol, atom)         -> std::integral;
get_ring_degree(mol, atom)        -> std::integral;
null_atom(mol)                    -> Atom; // e.g. nullptr, -1, ...

get_index(mol, bond)        -> std::integral;
get_source(mol, bond)       -> Atom;
get_target(mol, bond)       -> Atom;
get_order(mol, bond)        -> std::integral;
is_ring_bond(mol, bond)     -> bool;
is_aromatic_bond(mol, bond) -> bool;
null_bond(mol)              -> Bond;  // e.g. nullptr, -1, ...
```

## Modules

| Name          | Description                               | Dependencies        |
|---------------|-------------------------------------------|---------------------|
| Molecule      | Generic molecule interface                | -                   |
| CTSmarts      | Compile time SMARTS expressions           | Molecule, ctre      |
| Util          | Utility functions                         | fmt                 |
| CTLayout      | Compile time serialized data structures   | Util, Molecule, mio |
| OpenBabel     | Molecule interface for OpenBabel          | Molecule, OpenBabel |
| RDKit         | Molecule interface for RDKit              | Molecule, RDKit     |

## Thanks

For their collaboration and mentorship:

- Marcus Hanwell
- Geoffrey Hutchison
- Craig James
- Chris Morley
- Noel O'Boyle

For their inspiring ideas:

- Hana Dusíková
- Greg Landrum
- John Mayfield
- Roger Sayle

This project is dedicated to my wife
