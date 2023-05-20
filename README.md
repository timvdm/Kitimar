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



CTLayout

- Types
    - Value
        - fixed size
    - Struct
        - fixed size members -> fixed size
        - variable size members -> store size
    - Array
        - store length
        - fixed size value type -> store stride
        - variable size value type -> offset table
    

- Special Types
    - ArraySize<Array>

- isX (X = Value, Struct, ...)
- isFixedSize
- sizeOf







TODO
====

- Move Isomorphism API to Smarts
- Test case names/files/targets: TestFoo
- CTLayout recursive arrays?
- Add MolSource concept to Molecule
