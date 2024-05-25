add_executable(unit_tests_mathkit
               double3x3.cpp
               main.cpp)

target_link_libraries(unit_tests_mathkit PRIVATE libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GTEST_LIBRARIES} -lc++ -lc++abi -lpthread)

add_test(NAME MathKit COMMAND unit_tests_mathkit)
