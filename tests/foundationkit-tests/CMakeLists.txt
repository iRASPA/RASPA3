add_executable(unit_tests_foundationkit
  archive.cpp
  main.cpp)

raspa_target_link_libraries(unit_tests_foundationkit)
target_link_libraries(unit_tests_foundationkit PUBLIC
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

gtest_add_tests(TARGET unit_tests_foundationkit TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"
      \"tests/foundationkit-tests/unit_tests_foundationkit\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_foundationkit DESTINATION share/raspa3/tests COMPONENT unit_tests)
