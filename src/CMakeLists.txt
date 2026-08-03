project(raspakit)

# The three kits below the simulation engine, each its own library and each stacked on the one before it:
# containers and file formats, then vector and matrix arithmetic, then crystal symmetry. They know nothing
# of systems, components or moves, and both the engine and the structural analysis are built on top of them.
# Keeping them apart means a new import that points the wrong way down the stack fails to link rather than
# passing unnoticed, and that a change to the symmetry search does not rebuild anything that only wanted a
# vector type.
add_subdirectory(foundationkit)
add_subdirectory(mathkit)
add_subdirectory(symmetrykit)

# The SIMD variants below compile these sources into themselves rather than linking the three libraries, so
# the file lists come back up from those directories.
set(BASE_CPP_FILES ${FOUNDATIONKIT_CPP_FILES} ${MATHKIT_CPP_FILES} ${SYMMETRYKIT_CPP_FILES})
set(BASE_IXX_FILES ${FOUNDATIONKIT_IXX_FILES} ${MATHKIT_IXX_FILES} ${SYMMETRYKIT_IXX_FILES})

file(GLOB SRC_CPP_FILES
  ${PROJECT_SOURCE_DIR}/raspakit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/drivers/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/minimization/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/component/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/system/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/group/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/cross_system/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/intramolecular/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/interpolation/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/properties_time_evolution/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/integrators/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/ruptura/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/utilities/*.cpp)

file(GLOB SRC_IXX_FILES
  ${PROJECT_SOURCE_DIR}/raspakit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/drivers/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/minimization/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/component/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/system/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/group/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/cross_system/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/intramolecular/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/interpolation/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/properties_time_evolution/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/integrators/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/ruptura/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/utilities/*.ixx)

add_library(${PROJECT_NAME}_base STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME}_base PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME}_base PRIVATE ${SRC_CPP_FILES})
target_sources(${PROJECT_NAME}_base PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
target_compile_options(${PROJECT_NAME}_base PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME}_base PUBLIC symmetrykit mathkit foundationkit OpenMP::OpenMP_CXX)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")

  add_library(${PROJECT_NAME}_avx STATIC ${PCH_HEADER} ${PCH_SOURCE})
  set_target_properties(${PROJECT_NAME}_avx PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  target_sources(${PROJECT_NAME}_avx PRIVATE ${BASE_CPP_FILES} ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx PUBLIC FILE_SET CXX_MODULES FILES ${BASE_IXX_FILES} ${SRC_IXX_FILES})
  target_compile_options(${PROJECT_NAME}_avx PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  target_link_libraries(${PROJECT_NAME}_avx PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(${PROJECT_NAME}_avx PRIVATE /arch:AVX)
  else()
    # -mno-avx512f needed for Google-Colab (fixes "illegal instruction")
    target_compile_options(${PROJECT_NAME}_avx PRIVATE -march=x86-64-v2 -mno-avx512f -mavx)
  endif()

  add_library(${PROJECT_NAME}_avx2 STATIC ${PCH_HEADER} ${PCH_SOURCE})
  set_target_properties(${PROJECT_NAME}_avx2 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  target_sources(${PROJECT_NAME}_avx2 PRIVATE ${BASE_CPP_FILES} ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx2 PUBLIC FILE_SET CXX_MODULES FILES ${BASE_IXX_FILES} ${SRC_IXX_FILES})
  target_compile_options(${PROJECT_NAME}_avx2 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  target_link_libraries(${PROJECT_NAME}_avx2 PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(${PROJECT_NAME}_avx2 PRIVATE /arch:AVX2)
  else()
    # -mno-avx512f needed for Google-Colab (fixes "illegal instruction")
    target_compile_options(${PROJECT_NAME}_avx2 PRIVATE -march=x86-64-v3 -mno-avx512f -mavx2 -mfma)
  endif()

  add_library(${PROJECT_NAME}_avx512 STATIC ${PCH_HEADER} ${PCH_SOURCE})
  set_target_properties(${PROJECT_NAME}_avx512 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  target_sources(${PROJECT_NAME}_avx512 PRIVATE ${BASE_CPP_FILES} ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx512 PUBLIC FILE_SET CXX_MODULES FILES ${BASE_IXX_FILES} ${SRC_IXX_FILES})
  target_compile_options(${PROJECT_NAME}_avx512 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  target_link_libraries(${PROJECT_NAME}_avx512 PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE /arch:AVX512)
  else()
    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE -march=x86-64-v4 -mavx512f -mfma -mprefer-vector-width=512)
  endif()

endif()

# The structural analysis lives beside this library rather than above it, on the same three kits, and is
# built once rather than once per SIMD level. The adapter that turns a framework and a force field into what
# it takes is the only thing that depends on both.
add_subdirectory(structurekit)
add_subdirectory(structure_input)
