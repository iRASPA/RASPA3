add_executable(unit_tests_raspakit
  cbmc_interactions.cpp
  energy_decomposition.cpp
  gradients.cpp
  integrators.cpp
  pbc.cpp
  spacegroup.cpp
  dudlambda.cpp
  ewald.cpp
  insertion_deletion.cpp
  mc.cpp
  pressure.cpp
  static_energy.cpp
  rigid_gradient.cpp
  electrostatic_potential.cpp
  electrostatic_field.cpp
  electrostatic_polarization.cpp
  nist_water_1.cpp
  nist_water_2.cpp
  nist_water_3.cpp
  nist_water_4.cpp
  hessian_inter_lennard_jones.cpp
  hessian_inter_real_ewald.cpp
  tricubic_derivative_inter_lennard_jones.cpp
  triquintic_derivative_inter_lennard_jones.cpp
  tricubic_derivative_inter_real_ewald.cpp
  triquintic_derivative_inter_real_ewald.cpp
  intramolecular_energy_gradient_strain_bonds.cpp
  intramolecular_energy_gradient_strain_bends.cpp
  rotation_matrix.cpp
  main.cpp)

raspa_target_link_libraries(unit_tests_raspakit)
target_link_libraries(unit_tests_raspakit PUBLIC
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

gtest_add_tests(TARGET unit_tests_raspakit TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/raspakit-tests/unit_tests_raspakit\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_raspakit DESTINATION share/raspa3/tests COMPONENT unit_tests)
