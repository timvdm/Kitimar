add_library(Util INTERFACE
    Util.hpp
    Test.hpp
    Timer.hpp
)
target_include_directories(Util INTERFACE ${CMAKE_SOURCE_DIR}/src)
