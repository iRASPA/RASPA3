# CMakeLists.txt
cmake_minimum_required(VERSION 3.28)

project(raspa3
  VERSION 3.0.0
  DESCRIPTION ""
  HOMEPAGE_URL "https://github.com/iraspa/raspa3"
  LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(LINUX)
  execute_process (
    COMMAND bash -c "llvm-config --libdir"
    OUTPUT_VARIABLE LLVM_LIBS
  )
  
  execute_process (
    COMMAND bash -c "llvm-config --version"
    OUTPUT_VARIABLE LLVM_VERSION
  )
elseif(APPLE)
  execute_process (
    COMMAND bash -c "brew --prefix llvm"
    OUTPUT_VARIABLE BREW_LLVM_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process (
    COMMAND bash -c "brew --prefix libomp"
    OUTPUT_VARIABLE BREW_LIBOMP_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process (
    COMMAND bash -c "brew --prefix hdf5"
    OUTPUT_VARIABLE BREW_HDF5_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process (
    COMMAND bash -c "brew --prefix libaec"
    OUTPUT_VARIABLE BREW_LIBAEC_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
elseif(WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "AMD64")
  find_package(LLVM CONFIG REQUIRED)
  list(APPEND CMAKE_MODULE_PATH "{LLVM_CMAKE_DIR}")
endif()

find_package(ZLIB)
find_package(BLAS)
find_package(LAPACK)

if(NOT WIN32)
  find_package(OpenMP REQUIRED)
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

# Set the version of C++ for the project
set(CMAKE_CXX_STANDARD 23)

# Turning off extensions avoids an issue with the clang 16 compiler
# clang 17 and greater can avoid this setting
set(CMAKE_CXX_EXTENSIONS OFF)

add_definitions(-DUSE_LEGACY_HEADERS -DVERSION=${PROJECT_VERSION})

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
add_subdirectory(data)
add_subdirectory(examples)

if(NOT DEFINED SKIP_EXECUTABLE)
  add_subdirectory(app)
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
