# CMakeLists.txt
cmake_minimum_required(VERSION 3.28)

project(raspa3
  VERSION 3.0.0
  DESCRIPTION ""
  HOMEPAGE_URL "https://github.com/iraspa/raspa3"
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(BLAS)
find_package(LAPACK)

if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

find_package(OpenMP REQUIRED)

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

# -stdlib=libc++
set(CMAKE_CXX_FLAGS_DEBUG "-stdlib=libc++ -fvisibility=hidden -fexperimental-library ${OpenMP_CXX_FLAGS} -g3 -O0 -DDEBUG -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-stdlib=libc++ -fvisibility=hidden -fexperimental-library ${OpenMP_CXX_FLAGS} -g -O3 -march=core-avx2 -mfma -fomit-frame-pointer -ffast-math")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-stdlib=libc++ -fvisibility=hidden -fexperimental-library ${OpenMP_CXX_FLAGS} -g -O3 -march=core-avx2 -mfma -fomit-frame-pointer -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-stdlib=libc++ -fvisibility=hidden -fexperimental-library ${OpenMP_CXX_FLAGS} -g -O3 -march=core-avx2 -mfma -fomit-frame-pointer -ffast-math")

add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(data)

if(BUILD_DOXYGEN)
  add_subdirectory(docs)
endif()

add_subdirectory(packaging)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(BUILD_PYTHON)
  add_subdirectory(raspa3)
endif()
