
# Locate GTest
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

enable_testing()

## Set the version of C++ for the project
#set(CMAKE_CXX_STANDARD 23)
#
## Turning off extensions avoids an issue with the clang 16 compiler
## clang 17 and greater can avoid this setting
#set(CMAKE_CXX_EXTENSIONS OFF)
#
#set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DDEBUG")
#
#add_compile_options(-Wall -Wextra)
#
#add_definitions(-DUSE_LEGACY_HEADERS)

add_subdirectory(mathkit-tests)
add_subdirectory(foundationkit-tests)
add_subdirectory(symmetrykit-tests)
add_subdirectory(raspakit-tests)
