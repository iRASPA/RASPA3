add_executable(unit_tests_symmetrykit
  band_path.cpp
  niggli_reduction.cpp
  poscar_parser.cpp
  cif_symmetry.cpp
  generated_spacegroups.cpp
  space_group_database.cpp
  main.cpp)

raspa_target_link_libraries(unit_tests_symmetrykit)
target_link_libraries(unit_tests_symmetrykit PUBLIC
                      symmetrykit
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

# the generator of testdata/spacegroups (run it from this directory to regenerate the structures)
add_executable(generate_symmetrykit_testdata generate_testdata.cpp)
raspa_target_link_libraries(generate_symmetrykit_testdata)
target_link_libraries(generate_symmetrykit_testdata PUBLIC symmetrykit)

file(COPY testdata DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

gtest_add_tests(TARGET unit_tests_symmetrykit TEST_LIST noArgsTests)
set_tests_properties(${noArgsTests} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# install(TARGETS unit_tests_symmetrykit DESTINATION share/raspa3/tests COMPONENT unit_tests)
