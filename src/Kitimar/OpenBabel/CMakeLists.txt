

find_package(OpenBabel3)



add_library(KitimarOpenBabel INTERFACE OpenBabel.hpp)
target_include_directories(KitimarOpenBabel INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(KitimarOpenBabel INTERFACE Molecule openbabel)




add_executable(OpenBabelTest Test.cpp)
target_link_libraries(OpenBabelTest KitimarOpenBabel)
