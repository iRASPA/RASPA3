add_executable(unit_tests_mathkit
               double3x3.cpp
               main.cpp)

target_link_libraries(unit_tests_mathkit PRIVATE libs::raspakit ${GTEST_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(NAME MathKit COMMAND unit_tests_mathkit)
