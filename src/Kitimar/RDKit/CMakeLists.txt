#set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
#set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)
find_package(RDKit)

add_library(KitimarRDKit INTERFACE RDKit.hpp)
target_include_directories(KitimarRDKit INTERFACE
    ${CMAKE_SOURCE_DIR}/src
    ${RDKIT_INCLUDE_DIR}
)
target_compile_definitions(KitimarRDKit INTERFACE -DKITIMAR_WITH_RDKIT)
target_link_libraries(KitimarRDKit INTERFACE Molecule ${RDKIT_LIBRARIES}) # FIXME add RDKit...



