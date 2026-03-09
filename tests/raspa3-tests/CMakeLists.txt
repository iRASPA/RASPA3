add_executable(unit_tests_raspa3
  mc_muvt_drift.cpp
  mc_muvt_grids_drift.cpp
  mc_nvt_drift.cpp
  mc_nvt_grids_drift.cpp
  main.cpp)

raspa_target_link_libraries(${PROJECT_NAME})

gtest_add_tests(TARGET unit_tests_raspa3 TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/raspa3-tests/unit_tests_raspa3\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_raspa3 DESTINATION share/raspa3/tests COMPONENT unit_tests)
