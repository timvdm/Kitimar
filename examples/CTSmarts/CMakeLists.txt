add_executable(Match Match.cpp)
target_link_libraries(Match Kitimar::CTSmarts)

add_executable(Mapping Mapping.cpp)
target_link_libraries(Mapping Kitimar::CTSmarts)

add_executable(Capture Capture.cpp)
target_link_libraries(Capture Kitimar::CTSmarts)

add_executable(AtomBond AtomBond.cpp)
target_link_libraries(AtomBond Kitimar::CTSmarts)
