cmake_minimum_required(VERSION 3.16)

project(Kitimar LANGUAGES CXX)

option(KITIMAR_WITH_OPENBABEL "Build with OpenBabel support" ON)
option(KITIMAR_WITH_RDKIT "Build with RDKit support" OUTPUT_NAME)
option(KITIMAR_WITH_COROUTINES "Build with c++20 coroutines support" ON)

option(KITIMAR_WITH_VALIDATION "Build validation executable (WARNING: uses lots of RAM during compilation)" OFF)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall)
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Werror)
else()
    add_compile_options(-fbracket-depth=1024)
endif()
add_compile_options(-pedantic)
add_compile_options(-Wno-sign-compare)
add_compile_options(-Wno-unused-variable)
add_compile_options(-Wno-unused-but-set-variable)

#add_compile_options(-fconcepts-diagnostics-depth=3)

add_compile_options(-march=native -mtune=native)
#set(CMAKE_UNITY_BUILD ON)

#add_compile_options(-flto)



#add_compile_definitions(_GLIBCXX_DEBUG)
#add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
#add_link_options(-fsanitize=address)


include(cmake/CPM.cmake)

CPMAddPackage("gh:fmtlib/fmt#10.0.0")
CPMAddPackage("gh:hanickadot/compile-time-regular-expressions@3.7.2")
CPMAddPackage("gh:mandreyel/mio#8b6b7d878c89e81614d05edca7936de41ccdd2da")
CPMAddPackage("gh:google/googletest@1.13.0")

CPMAddPackage(
    NAME benchmark
    GITHUB_REPOSITORY google/benchmark
    VERSION 1.8.0
    OPTIONS "BENCHMARK_ENABLE_TESTING OFF"
)

if(KITIMAR_WITH_COROUTINES)
    add_compile_definitions(KITIMAR_WITH_COROUTINES)
    CPMAddPackage(
        NAME cppcoro
        GITHUB_REPOSITORY andreasbuhr/cppcoro
        GIT_TAG 10bbcdbf2be3ad3aa56febcf4c7662d771460a99
        OPTIONS "BUILD_TESTING OFF"
    )
endif()

add_custom_target(AdditionalFiles ALL SOURCES
    README.md
)

include(GoogleTest)
enable_testing()
add_custom_target(RunTests ctest)


add_subdirectory(src)
add_subdirectory(test)
if(KITIMAR_WITH_VALIDATION)
    add_subdirectory(validate)
endif()
add_subdirectory(tools)
add_subdirectory(examples)
add_subdirectory(single-header)
