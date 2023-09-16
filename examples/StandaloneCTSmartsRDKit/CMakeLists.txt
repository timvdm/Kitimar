cmake_minimum_required(VERSION 3.16)
project(StandaloneCTSmartsRDKit)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Kitimar
find_package(Kitimar COMPONENTS CTSmarts RDKit REQUIRED)

# Executable
add_executable(Match Match.cpp)
target_link_libraries(Match Kitimar::CTSmarts Kitimar::RDKit)
