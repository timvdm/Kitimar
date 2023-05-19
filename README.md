Kitimar
=======


Name            Description                                 Dependencies
------------------------------------------------------------------------
Util            Utility functions                           -
Molecule        Generic molecule interface                  -
CTSmarts        Compile time SMARTS matching                Molecule, ctre
CTLayout        Compile time data structures                ctre, mio
Serialize       Molecule serialization using CTLayout       Molecule, CTLayout

OpenBabel       Molecule interface for OpenBabel            Molecule, openbabel





- readMoleculesOpenBabel -> OpenBabelSmilesMolSource
- remove readMoleculesOpenBabel
- modify serializeChemblKitimar(source)
- add RDKit


- add MolSource concept to Molecule





TODO
====

- Move Isomorphism API to Smarts
- Test case names/files/targets: TestFoo
- CTLayout recursive arrays?
