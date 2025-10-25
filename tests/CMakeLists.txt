if(BUILD_TESTING_MATHKIT)
add_subdirectory(mathkit-tests)
endif()

if(BUILD_TESTING_FOUNDATIONKIT)
add_subdirectory(foundationkit-tests)
endif()

if(BUILD_TESTING_SYMMETRYKIT)
add_subdirectory(symmetrykit-tests)
endif()

if(BUILD_TESTING_RASPAKIT)
add_subdirectory(raspakit-tests)
endif()

if(BUILD_TESTING_RASPA3)
add_subdirectory(raspa3-tests)
endif()
