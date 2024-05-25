add_executable(unit_tests_foundationkit
               archive.cpp 
               main.cpp)

target_link_libraries(unit_tests_foundationkit PRIVATE libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GTEST_LIBRARIES} -lc++ -lc++abi -lpthread)

add_test(NAME FoundationKit COMMAND unit_tests_foundationkit)
