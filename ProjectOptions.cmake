include(cmake/SystemLink.cmake)
include(cmake/LibFuzzer.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(Kitimar_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(Kitimar_setup_options)
  option(Kitimar_ENABLE_HARDENING "Enable hardening" ON)
  option(Kitimar_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    Kitimar_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    ON
    Kitimar_ENABLE_HARDENING
    OFF)

  Kitimar_supports_sanitizers()

  if(NOT PROJECT_IS_TOP_LEVEL OR Kitimar_PACKAGING_MAINTAINER_MODE)
    option(Kitimar_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(Kitimar_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(Kitimar_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(Kitimar_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(Kitimar_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(Kitimar_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(Kitimar_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(Kitimar_ENABLE_PCH "Enable precompiled headers" OFF)
    option(Kitimar_ENABLE_CACHE "Enable ccache" OFF)
  else()
    option(Kitimar_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(Kitimar_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(Kitimar_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(Kitimar_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" ${SUPPORTS_ASAN})
    option(Kitimar_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" ${SUPPORTS_UBSAN})
    option(Kitimar_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(Kitimar_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(Kitimar_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(Kitimar_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(Kitimar_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(Kitimar_ENABLE_PCH "Enable precompiled headers" OFF)
    option(Kitimar_ENABLE_CACHE "Enable ccache" ON)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      Kitimar_ENABLE_IPO
      Kitimar_WARNINGS_AS_ERRORS
      Kitimar_ENABLE_USER_LINKER
      Kitimar_ENABLE_SANITIZER_ADDRESS
      Kitimar_ENABLE_SANITIZER_LEAK
      Kitimar_ENABLE_SANITIZER_UNDEFINED
      Kitimar_ENABLE_SANITIZER_THREAD
      Kitimar_ENABLE_SANITIZER_MEMORY
      Kitimar_ENABLE_UNITY_BUILD
      Kitimar_ENABLE_CLANG_TIDY
      Kitimar_ENABLE_CPPCHECK
      Kitimar_ENABLE_COVERAGE
      Kitimar_ENABLE_PCH
      Kitimar_ENABLE_CACHE)
  endif()

  Kitimar_check_libfuzzer_support(LIBFUZZER_SUPPORTED)
  if(LIBFUZZER_SUPPORTED AND (Kitimar_ENABLE_SANITIZER_ADDRESS OR Kitimar_ENABLE_SANITIZER_THREAD OR Kitimar_ENABLE_SANITIZER_UNDEFINED))
    set(DEFAULT_FUZZER ON)
  else()
    set(DEFAULT_FUZZER OFF)
  endif()

  option(Kitimar_BUILD_FUZZ_TESTS "Enable fuzz testing executable" ${DEFAULT_FUZZER})

endmacro()

macro(Kitimar_global_options)
  if(Kitimar_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    Kitimar_enable_ipo()
  endif()

  Kitimar_supports_sanitizers()

  if(Kitimar_ENABLE_HARDENING AND Kitimar_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR Kitimar_ENABLE_SANITIZER_UNDEFINED
       OR Kitimar_ENABLE_SANITIZER_ADDRESS
       OR Kitimar_ENABLE_SANITIZER_THREAD
       OR Kitimar_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    message("${Kitimar_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME} ${Kitimar_ENABLE_SANITIZER_UNDEFINED}")
    Kitimar_enable_hardening(KitimarOptions ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(Kitimar_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(KitimarWarnings INTERFACE)
  add_library(KitimarOptions INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  Kitimar_set_project_warnings(
    KitimarWarnings
    ${Kitimar_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(Kitimar_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(KitimarOptions)
  endif()

  include(cmake/Sanitizers.cmake)
  Kitimar_enable_sanitizers(
    KitimarOptions
    ${Kitimar_ENABLE_SANITIZER_ADDRESS}
    ${Kitimar_ENABLE_SANITIZER_LEAK}
    ${Kitimar_ENABLE_SANITIZER_UNDEFINED}
    ${Kitimar_ENABLE_SANITIZER_THREAD}
    ${Kitimar_ENABLE_SANITIZER_MEMORY})

  set_target_properties(KitimarOptions PROPERTIES UNITY_BUILD ${Kitimar_ENABLE_UNITY_BUILD})

  if(Kitimar_ENABLE_PCH)
    target_precompile_headers(
      KitimarOptions
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  if(Kitimar_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    Kitimar_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(Kitimar_ENABLE_CLANG_TIDY)
    Kitimar_enable_clang_tidy(KitimarOptions ${Kitimar_WARNINGS_AS_ERRORS})
  endif()

  if(Kitimar_ENABLE_CPPCHECK)
    Kitimar_enable_cppcheck(${Kitimar_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(Kitimar_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    Kitimar_enable_coverage(KitimarOptions)
  endif()

  if(Kitimar_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(KitimarOptions INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(Kitimar_ENABLE_HARDENING AND NOT Kitimar_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR Kitimar_ENABLE_SANITIZER_UNDEFINED
       OR Kitimar_ENABLE_SANITIZER_ADDRESS
       OR Kitimar_ENABLE_SANITIZER_THREAD
       OR Kitimar_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    Kitimar_enable_hardening(KitimarOptions OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()
