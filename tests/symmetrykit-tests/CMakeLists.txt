add_executable(unit_tests_symmetrykit
  transformation_matrix.cpp
  find_smallest_primitive_cell.cpp
  find_smallest_primitive_cell_no_partial_occupancies.cpp
  find_pointgroup.cpp
  find_pointgroup_no_partial_occupancies.cpp
  find_spacegroup.cpp
  main.cpp)

raspa_target_link_libraries(unit_tests_symmetrykit)
target_link_libraries(unit_tests_symmetrykit PUBLIC
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

file(COPY spglibtestdata DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

gtest_add_tests(TARGET unit_tests_symmetrykit TEST_LIST noArgsTests)

# install(TARGETS unit_tests_symmetrykit DESTINATION share/raspa3/tests COMPONENT unit_tests)
