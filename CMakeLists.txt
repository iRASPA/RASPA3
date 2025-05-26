# CMakeLists.txt
cmake_minimum_required(VERSION 3.29)

project(raspa3
  VERSION 3.0.5
  DESCRIPTION ""
  HOMEPAGE_URL "https://github.com/iraspa/raspa3"
  LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


# Set the version of C++ for the project
set(CMAKE_CXX_STANDARD 26)

# Turning off extensions avoids an issue with the clang 16 compiler
# clang 17 and greater can avoid this setting
set(CMAKE_CXX_EXTENSIONS OFF)

add_definitions( -DUSE_LEGACY_HEADERS -DVERSION=${PROJECT_VERSION})

find_package(ZLIB)
find_package(BLAS)
find_package(LAPACK)
find_package(OpenMP)
find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
find_package(OpenCL)

if(BUILD_LIBTORCH)
  find_package(Torch REQUIRED)
  add_definitions( -DBUILD_LIBTORCH -Wno-invalid-specialization)
  message("TORCH_LIBRARIES: ${TORCH_LIBRARIES}")
  include_directories(
    "${TORCH_INCLUDE_DIRS}" # libtorch headers
  )
endif()

include(FetchContent)
FetchContent_Declare(voro++
      GIT_REPOSITORY https://github.com/chr1shr/voro.git)
FetchContent_MakeAvailable(voro++)

set(MSVC_WARNINGS /W4 /w14242 /w14254 /w14263 /w14265 /w14287 /we4289 /w14296 /w14311 /w14545 /w14546 /w14547 /w14549 /w14555 /w14619 /w14640 /w14826 /w14826 /w14905 /w14906 /w14928 /permissive-)

set(CLANG_WARNINGS  -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wimplicit-fallthrough -Wno-gnu-anonymous-struct -Wno-error=deprecated-declarations -Wno-error=nan-infinity-disabled -Wno-deprecated-declarations -Wno-nan-infinity-disabled -Wno-unknown-warning-option)

set(GCC_WARNINGS ${CLANG_WARNINGS} -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast -Wsuggest-override)

if(MSVC)
  set(PROJECT_WARNINGS_CXX ${MSVC_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  set(PROJECT_WARNINGS_CXX ${CLANG_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PROJECT_WARNINGS_CXX ${GCC_WARNINGS})
endif()

add_subdirectory(src)

if(NOT SKIP_EXECUTABLE)
  add_subdirectory(data)
  add_subdirectory(examples)

  add_subdirectory(app)
  add_subdirectory(cli)
endif()
if(BUILD_DOXYGEN)
  add_subdirectory(docs)
endif()

add_subdirectory(packaging)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

if(BUILD_PYTHON)
  add_subdirectory(raspa)
endif()
