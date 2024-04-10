find_package(GTest)
if (GTest_FOUND)
  include_directories(${GTEST_INCLUDE_DIRS})
else ()
  include(FetchContent)
  FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0)
  
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

include(CTest)

add_subdirectory(mathkit-tests)
add_subdirectory(foundationkit-tests)
add_subdirectory(symmetrykit-tests)
add_subdirectory(raspakit-tests)
