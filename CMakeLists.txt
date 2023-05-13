cmake_minimum_required(VERSION 3.16)

project(Kitimar LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FIXME: local settings...
set(CMAKE_BUILD_RPATH "/usr/local/lib64")
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

add_compile_options(-Wall)
add_compile_options(-Wno-sign-compare)
add_compile_options(-Wno-unused-variable)
add_compile_options(-Wno-unused-but-set-variable)

add_compile_options(-fconcepts-diagnostics-depth=3)

add_compile_options(-march=native -mtune=native)
#add_compile_options(-flto)

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)






include(FetchContent)

FetchContent_Declare(
  ctre
  GIT_REPOSITORY https://github.com/hanickadot/compile-time-regular-expressions.git
  GIT_TAG v3.7.2
)
FetchContent_MakeAvailable(ctre)

FetchContent_Declare(
  mio
  GIT_REPOSITORY https://github.com/mandreyel/mio.git
)
FetchContent_MakeAvailable(mio)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.13.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)





add_custom_target(README ALL SOURCES README.md)

include(GoogleTest)
enable_testing()

add_custom_target(RunTests ctest)



add_subdirectory(src)
