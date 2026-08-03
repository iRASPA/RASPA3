add_executable(unit_tests_structurekit
  grid_area_curve.cpp
  grid_blocking_cover.cpp
  grid_percolation.cpp
  grid_pore_size.cpp
  structure_parallel.cpp
  exact_sphere_sweep.cpp
  exact_union_volume.cpp
  exact_void_split.cpp
  exact_solvent_excluded.cpp
  exact_pore_size_distribution.cpp
  boundary_components.cpp
  apollonius_surface_area.cpp
  main.cpp)

# Only the structural analysis and the three kits below it. What the engine can hand it is tested from
# raspakit-tests, where the adapter that translates a framework and a force field lives; keeping raspakit off
# this target is what makes a wrong-way import here fail to link rather than pass unnoticed.
raspa_target_link_libraries(unit_tests_structurekit)
target_link_libraries(unit_tests_structurekit PUBLIC
                      structurekit
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

gtest_add_tests(TARGET unit_tests_structurekit TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/structurekit-tests/unit_tests_structurekit\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_structurekit DESTINATION share/raspa3/tests COMPONENT unit_tests)
