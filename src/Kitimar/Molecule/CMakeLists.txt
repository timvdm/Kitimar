

add_library(Molecule INTERFACE Molecule.hpp)
target_include_directories(Molecule INTERFACE ${CMAKE_SOURCE_DIR}/src)



add_executable(MoleculeTest Test.cpp)
target_link_libraries(MoleculeTest Molecule)
