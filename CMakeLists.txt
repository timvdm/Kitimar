cmake_minimum_required(VERSION 3.16)

project(Kitimar LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall)
add_compile_options(-Werror)
add_compile_options(-pedantic)
add_compile_options(-Wno-sign-compare)
add_compile_options(-Wno-unused-variable)
add_compile_options(-Wno-unused-but-set-variable)

#add_compile_options(-fconcepts-diagnostics-depth=3)

add_compile_options(-march=native -mtune=native)
#set(CMAKE_UNITY_BUILD ON)

#add_compile_options(-flto)



#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)


include(cmake/CPM.cmake)

CPMAddPackage("gh:fmtlib/fmt#10.0.0")
CPMAddPackage("gh:hanickadot/compile-time-regular-expressions@3.7.2")
CPMAddPackage("gh:mandreyel/mio#8b6b7d878c89e81614d05edca7936de41ccdd2da")
CPMAddPackage("gh:google/googletest@1.13.0")

add_custom_target(AdditionalFiles ALL SOURCES
    README.md
)

include(GoogleTest)
enable_testing()
add_custom_target(RunTests ctest)

option(KITIMAR_WITH_OPENBABEL "Build with OpenBabel support" OFF)
option(KITIMAR_WITH_RDKIT "Build with RDKit support" OFF)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(tools)
add_subdirectory(examples)
