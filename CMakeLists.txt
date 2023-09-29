# CMake setup based on https://github.com/cpp-best-practices/cmake_template

cmake_minimum_required(VERSION 3.16)

project(Kitimar
    VERSION 0.1.0
    DESCRIPTION "Kitimar: C++20 Cheminformatics"
    HOMEPAGE_URL "https://github.com/timvdm/Kitimar"
    LANGUAGES CXX
)

option(Kitimar_WITH_OPENBABEL "Build with OpenBabel support" OFF)
option(Kitimar_WITH_RDKIT "Build with RDKit support" OFF)
option(Kitimar_WITH_COROUTINES "Build with c++20 coroutines support" OFF)
option(Kitimar_WITH_VALIDATION_BASIC "Build basic validation" OFF)
option(Kitimar_WITH_VALIDATION_SQC "Build substructure query collection validation (WARNING: uses lots of RAM during compilation)" OFF)

# Including <iostream> results in code being added to empty executables.
# Without this, generated assembly code is easier to implement.
option(Kitimar_WITH_IOSTREAM "Use STL iostream" ON)

# Only set the cxx_standard if it is not set by someone else
if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# Strongly encouraged to enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20 for example
# when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

# Global configuration
include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)
Kitimar_setup_options()
Kitimar_global_options()

# Dependencies
include(Dependencies.cmake)
Kitimar_setup_dependencies()

# Local configuration
Kitimar_local_options()

target_compile_features(KitimarOptions INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

set_target_properties(KitimarOptions PROPERTIES EXPORT_NAME "Options")
set_target_properties(KitimarWarnings PROPERTIES EXPORT_NAME "Warnings")

add_library(Kitimar::Options ALIAS KitimarOptions)
add_library(Kitimar::Warnings ALIAS KitimarWarnings)

if(Kitimar_WITH_IOSTREAM)
    target_compile_definitions(KitimarOptions INTERFACE -DKITIMAR_WITH_IOSTREAM)
endif()


if(UNIX)

    target_compile_options(KitimarOptions INTERFACE
        -march=native
        -mtune=native
        -pedantic
    )

    target_compile_options(KitimarWarnings INTERFACE
        #-Wno-unused-parameter
        #-Wno-sign-compare
        #-Wno-unused-but-set-variable
        #-Wno-missing-field-initializers
    )

    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(KitimarOptions INTERFACE -ftemplate-depth=1024)
    else()
        target_compile_options(KitimarOptions INTERFACE -fbracket-depth=1024)
        #-fconcepts-diagnostics-depth=5
    endif()

endif() # UNIX


install(TARGETS KitimarOptions
    EXPORT KitimarOptionsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT KitimarOptionsTargets
    FILE KitimarOptionsTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)

install(TARGETS KitimarWarnings EXPORT KitimarWarningsTargets)
install(EXPORT KitimarWarningsTargets
    FILE KitimarWarningsTargets.cmake
    NAMESPACE Kitimar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)






add_custom_target(AdditionalFiles ALL SOURCES
    README.md
    #.github/workflows/build_test.yml
    Config.cmake.in
    ProjectOptions.cmake
    Dependencies.cmake
    cmake/PreventInSourceBuilds.cmake
    cmake/SystemLink.cmake
    cmake/LibFuzzer.cmake
    cmake/InterproceduralOptimization.cmake
    cmake/Hardening.cmake
    cmake/StandardProjectSettings.cmake
    cmake/VCEnvironment.cmake
    cmake/Utilities.cmake
    cmake/CompilerWarnings.cmake
    cmake/Sanitizers.cmake
    cmake/Cache.cmake
    cmake/StaticAnalyzers.cmake
)




add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(examples)
add_subdirectory(single-header)

# Tests
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(benchmark)
    add_subdirectory(validate)
endif()


include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/KitimarConfigVersion.cmake
  VERSION ${Kitimar_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/KitimarConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)


install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/KitimarConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/KitimarConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kitimar
)

# Packaging

set(CPACK_SOURCE_PACKAGE_FILE_NAME "Kitimar-${Kitimar_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES
    "/build"
    ".git"
)

include(CPack)
