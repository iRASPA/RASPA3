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

if (LINUX)
  target_link_libraries(unit_tests_symmetrykit
                        PUBLIC
                        libs::raspakit
                        -static-libstdc++
                        -static-libgcc
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -lpthread
)
elseif (Win32)
  target_link_libraries(unit_tests_symmetrykit
                        PUBLIC
                        libs::raspakit
                        -static-libstdc++
                        -static-libgcc
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -lpthread
)
elseif (APPLE)
target_link_libraries(unit_tests_symmetrykit PRIVATE
                      libs::raspakit
                      ${BLAS_LIBRARIES}
                      ${LAPACK_LIBRARIES}
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main
                      -nostdlib++
                      ${LLVM_LIBS}/c++/libc++experimental.a
                      ${LLVM_LIBS}/c++/libc++experimental.a
                      ${LLVM_LIBS}/c++/libc++.a
                      ${LLVM_LIBS}/c++/libc++abi.a
                      ${LLVM_LIBS}/libunwind.a
                      ${LLVM_LIBS}/../../../libomp/${LLVM_VERSION}/lib/libomp.a
                      -lpthread
)
endif()



file(COPY spglibtestdata DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME SymmetryKit COMMAND unit_tests_symmetrykit)
