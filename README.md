# Kitimar

Kitimar is a container project for CTSmarts providing utilities for development, testing and benchmarking.

## CTSmarts: C++ compile-time SMARTS expressions

`CTSmarts` provides fast compile-time SMARTS expressions with support for matching, mapping and capturing.

The SMARTS expressions are parsed at compile-time resulting in a minimal amount of assembly code and need for dynamic data structures at run-time.

It is implemented as a C++20 header-only library and inspired by the [Compile time regular expressions](https://github.com/hanickadot/compile-time-regular-expressions) and uses it's `ctll` compile-time LL(1) parser.

### Proof of concept

This is currently a working proof of concept. The aim was to see if this would possible in C++20 using current compilers.

### Features

- No parsing/optimization overhead
- Performance nearly identical to handwritten code
- Capture atoms using SMARTS atom classes
- No need for global state &rarr; no threading issues
- No need for external SMARTS to C++ conversion tool and build system setup

### API

The code below is a simplified version of the real API which makes it easier to read form a users perspective. Examples of how to use the API can be found in the `examples/CTSmarts` directory.

```c++
//
// Context
//
#include <Kitimar/CTSmarts/CTSmarts.hpp>

using namespace Kitimar;

auto mol = ...; // Molecule that satifies the Molecule::Molecule concept.

//
// API 
//
// <MapType> = CTSmarts::Unique or CTSmarts::All
// - Two mappings are considered unique if their atom set if different.
// - CTSmarts::Unique is always the default.

// Check if molecule contains SMARTS.

CTSmarts::contains<"SMARTS">(mol) -> bool;

// Find the first mapping or an empty vector if there is no mapping.

CTSmarts::single<"SMARTS">(mol) -> std::vector<int>;

// Find the first mapping starting from a specified atom or an empty vector if there is no mapping.

CTSmarts::single<"SMARTS">(mol) -> std::vector<int>;

// Count the number of mappings.

CTSmarts::count<"SMARTS">(mol, <MapType>) -> std::integral;

// Find multiple mappings.

CTSmarts::multi<"SMARTS">(mol, <MapType>) -> std::vector<std::vector<int>>;

// Capture the mapped atoms from a single mapping (starting from aspecified atom).

CTSmarts::capture<"SMARTS">(mol) -> std::tuple<bool, Atom...>;
CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>;

// Capture the mapped atoms from multiple mappings (starting from aspecified atom).

CTSmarts::captures<"SMARTS">(mol, <MapType>) -> std::range<std::array<Atom, N>>;
CTSmarts::captures<"SMARTS">(mol, atom, <MapType>) -> std::range<std::array<Atom, N>>;
```

#### Capture

The `CTSmarts::capture(s)` API is meant to be used with [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding) and [range-based for loops](https://en.cppreference.com/w/cpp/language/range-for).

```c++
// Capture all atoms
auto [match, C, O, N] = CTSmarts::capture<"C(=O)N">(mol);
if (match) {
   // Use C, O, N atoms...
}
```

```c++
// Capture atoms specified using SMARTS atom classes
auto [match, O, N] = CTSmarts::capture<"C(=[O:1])[N:2]">(mol);
if (match) {
   // Use O, N atoms...
}
```

```c++
// Capture all atoms for every unique mapping
for (auto [C, O] : CTSmarts::captures<"C=O">(mol)) {
   // Use C, O atoms...
}
```

### Memory usage

CTSmarts uses a ***minimal amount of memory*** at runtime.
Currently a `std::vector<bool>` is used to keep track of mapped atoms in the molecule.
When searching for unique matches, there is also a `std::unordered_set<std::size_t>` to keep track of found matches using a hash of the vector of bools.
An `std::array<int, NumSmartsAtoms>` is used to store the current mapping since the size is known at compile time. This avoids dynamic allocations.

NOTES: 
- See `Isomorphism` class
- `Isomorphism::m_degrees` will be removed
- Optimizations for special cases could be added to reduce this even further (see features).

### Planned features

- More compile-time optimizations
  - Match less common atoms first
- Improved error reporting
- Optimize special cases (match single atom/bond, central atom with bonds, ...)
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

The ***container*** was the easy part. While molecules are labeled graphs and no iterator pairs over a range of values,
a replacement for this was already available. The [Boost Graph Library (BGL)](boost.org/doc/libs/1_78_0/libs/graph/doc/index.html), 
which is also used by [RDKit](https://www.rdkit.org/), solved this by defining a graph concept (in documentation) and a set of free functions. 
Developers could implement these functions for their own custom graph data structure and use all of the algorithms provided by the BGL.
***Iterators*** are still iterators, but these iterate over the atoms/bonds of a molecule or incident/adjacent to another atom. Dereferencing an iterator no longer returns a value but an atom or a bond. 
These atoms/bonds have a vertex/edge label containing their properties (element, charge, bond order, ...) which can also be retrieved using a free function. 

The ***Molecule concept*** was now defined. In the documentation at least, C++98 did not have concepts. Older compilers would also print errors hundreds of lines long when you tried to call std::sort on an std::list.
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

The ***algorithms*** found in Helium include a number of basic graph algorithms (e.g. breath-first-search, depth-first-search, cycles, Dijkstra's shortest path, connected components, canonicalization, ...)
and cheminformatics algorithms (aromatize, kekulize, subgraph isomorphisms using SMARTS expressions, stereochemistry, ...).
These are ***high-level algorithms***, written using the same repetitive code I mentioned earlier.
These blocks of code are the low-level algorithms to query a molecule and could easily be expressed as SMARTS patterns.

Implementing these ***low-level algorithms*** efficiently with little to no overhead compared to the handwritten code was harder.
Parsing the ***SMARTS expression at run-time*** creates ***overhead***.
A shared instance between threads which is only initialized once requires synchronisation between threads.
Furthermore, this results in a graph representation and abstract syntax trees (AST) for the atom and bond expressions existing in memory. This layer of indirection will never be as fast as the handwritten code.

For simple cases, Helium tried to provide predefined predicates which could be used and combined to do some of these tasks.
This never resulted in an elegant syntax or provided any real advantage. More recently, I also tried to use the ranges library which is now part of the C++20 standard for this. 
Again, the results were never what I was looking for. SMARTS expressions are just a very compact way of representing these queries.
An equivalent way of expressing these ***queries using valid C++*** syntax will always be ***more verbose*** and requires a ***new, unfamiliar syntax*** to be learned.
Existing SMARTS expressions would also have to be converted to this syntax.

Another approach trying to combine the advantages of SMARTS with the performance of handwritten code are tools that ***generate code from SMARTS***.
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

I knew ***constexpr*** might be capable of this, but it was not immediately clear how this would work.
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
get_implicit_hydrogens(mol, atom) -> std::integral;
is_cyclic_atom(mol, atom)         -> bool;
is_aromatic_atom(mol, atom)       -> bool;
null_atom(mol)                    -> Atom; // e.g. nullptr, -1, ...

get_index(mol, bond)        -> std::integral;
get_source(mol, bond)       -> Atom;
get_target(mol, bond)       -> Atom;
get_order(mol, bond)        -> std::integral;
is_cyclic_bond(mol, bond)   -> bool;
is_aromatic_bond(mol, bond) -> bool;
null_bond(mol)              -> Bond;  // e.g. nullptr, -1, ...
```

## Modules

| Name          | Description                               | Dependencies        |
|---------------|-------------------------------------------|---------------------|
| Molecule      | Generic molecule interface                | -                   |
| CTSmarts      | Compile time SMARTS expressions           | Molecule, ctre      |
| Util          | Utility functions                         | fmt, ctre           |
| CTLayout      | Compile time serialized data structures   | Util, Molecule, mio |
| OpenBabel     | Molecule interface for OpenBabel          | Molecule, OpenBabel |
| RDKit         | Molecule interface for OpenBabel          | Molecule, RDKit     |

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
