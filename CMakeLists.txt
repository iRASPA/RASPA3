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
set(CMAKE_CXX_FLAGS_DEBUG "-fvisibility=hidden -fexperimental-library  -fopenmp -g3 -O0 -DDEBUG -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-fvisibility=hidden -fexperimental-library  -fopenmp -g -O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-fvisibility=hidden -fexperimental-library  -fopenmp -g -O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-fvisibility=hidden -fexperimental-library  -fopenmp -g -O3 -march=native -ffast-math")

# rpm-settings
set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR RPM)
set(CPACK_RPM_PACKAGE_NAME "raspa3")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DESCRIPTION "Molecular simulation code for adsorption and diffusion in nanoporours materials.")
set(CPACK_RPM_FILE_NAME "RPM-DEFAULT")
set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_CONTACT "David Dubbeldam (d.dubbeldam@uva.nl)")
set(CPACK_RPM_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_RPM_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
include(CPack)

add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(data)
add_subdirectory(docs)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(BUILD_PYTHON)
  add_subdirectory(raspa3)
endif()
