add_executable(unit_tests_mathkit
               double3x3.cpp
               main.cpp)

target_link_libraries(unit_tests_mathkit PRIVATE libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main -lc++ -lc++abi -lpthread)

add_test(NAME MathKit COMMAND unit_tests_mathkit)
