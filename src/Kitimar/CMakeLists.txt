add_subdirectory(Util)
add_subdirectory(Molecule)
add_subdirectory(CTLayout)
add_subdirectory(CTSmarts)
add_subdirectory(OpenBabel)

if(KITIMAR_WITH_RDKIT)
    add_subdirectory(RDKit)
endif()
