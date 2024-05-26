add_executable(unit_tests_raspakit
               cbmc_interactions.cpp
               energy_decomposition.cpp 
               gradients.cpp
               pbc.cpp
               spacegroup.cpp
               dudlambda.cpp
               ewald.cpp
               insertion_deletion.cpp
               mc.cpp
               pressure.cpp
               static_energy.cpp
               rigid_gradient.cpp
               main.cpp)

        target_link_libraries(unit_tests_raspakit libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} OpenMP::OpenMP_CXX GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main -lc++ -lc++abi -lpthread)

add_test(NAME RaspaKit COMMAND unit_tests_raspakit)
