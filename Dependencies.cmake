include(cmake/CPM.cmake)

# Done as a function so that updates to variables like
# CMAKE_CXX_FLAGS don't propagate out to other
# targets
function(Kitimar_setup_dependencies)

  # For each dependency, see if it's
  # already been provided to us by a parent project

  if(NOT TARGET fmtlib::fmtlib)
    CPMAddPackage("gh:fmtlib/fmt#10.0.0")
  endif()

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

endfunction()
