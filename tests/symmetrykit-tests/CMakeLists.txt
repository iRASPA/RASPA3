add_executable(unit_tests_symmetrykit
               transformation_matrix.cpp 
               find_smallest_primitive_cell.cpp 
               find_smallest_primitive_cell_no_partial_occupancies.cpp 
               find_pointgroup.cpp
               find_pointgroup_no_partial_occupancies.cpp 
               find_pointgroup_no_partial_occupancies.cpp 
               find_spacegroup.cpp 
               find_spacegroup_no_partial_occupancies.cpp
               main.cpp)

target_link_libraries(unit_tests_symmetrykit PRIVATE libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GTEST_LIBRARIES} -lc++ -lc++abi -lpthread)

file(COPY spglibtestdata DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME SymmetryKit COMMAND unit_tests_symmetrykit)
