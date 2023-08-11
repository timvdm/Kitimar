add_subdirectory(Util)
add_subdirectory(Molecule)
add_subdirectory(CTLayout)
add_subdirectory(CTSmarts)

if(Kitimar_WITH_OPENBABEL)
    add_subdirectory(OpenBabel)
endif()

if(Kitimar_WITH_RDKIT)
    add_subdirectory(RDKit)
endif()
