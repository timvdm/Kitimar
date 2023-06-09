add_executable(Match Match.cpp)
target_link_libraries(Match CTSmarts)

add_executable(Mapping Mapping.cpp)
target_link_libraries(Mapping CTSmarts)

add_executable(Capture Capture.cpp)
target_link_libraries(Capture CTSmarts)

add_executable(AtomBond AtomBond.cpp)
target_link_libraries(AtomBond CTSmarts)
