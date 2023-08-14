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

    if(NOT TARGET Catch2::Catch2WithMain)
        CPMAddPackage("gh:catchorg/Catch2@3.4.0")
    endif()

    if(NOT TARGET ctre::ctre)
        CPMAddPackage("gh:hanickadot/compile-time-regular-expressions@3.8")
    endif()

    if(NOT TARGET mio::mio)
        CPMAddPackage("gh:mandreyel/mio#8b6b7d878c89e81614d05edca7936de41ccdd2da")
    endif()

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
