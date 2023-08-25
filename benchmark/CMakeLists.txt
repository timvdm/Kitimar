if(Kitimar_WITH_OPENBABEL)
    add_executable(BenchmarkCTSmarts BenchmarkCTSmarts.cpp
        Benchmarks.hpp)
    target_link_libraries(BenchmarkCTSmarts Kitimar::CTSmarts Kitimar::OpenBabel Catch2::Catch2WithMain TestData)
endif()
