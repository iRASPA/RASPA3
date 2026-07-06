project(raspakit)

file(GLOB SRC_CPP_FILES
  ${PROJECT_SOURCE_DIR}/foundationkit/*.cpp
  ${PROJECT_SOURCE_DIR}/mathkit/*.cpp
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/rigid_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.cpp
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
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/grid_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/utilities/*.cpp)

file(GLOB SRC_IXX_FILES
  ${PROJECT_SOURCE_DIR}/foundationkit/*.ixx
  ${PROJECT_SOURCE_DIR}/mathkit/*.ixx
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/rigid_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.ixx
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
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/grid_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/utilities/*.ixx)

add_library(${PROJECT_NAME}_base STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME}_base PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME}_base PRIVATE ${SRC_CPP_FILES})
target_sources(${PROJECT_NAME}_base PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
target_compile_options(${PROJECT_NAME}_base PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME}_base PUBLIC OpenMP::OpenMP_CXX)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")

  add_library(${PROJECT_NAME}_avx STATIC ${PCH_HEADER} ${PCH_SOURCE})
  set_target_properties(${PROJECT_NAME}_avx PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  target_sources(${PROJECT_NAME}_avx PRIVATE ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
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
  target_sources(${PROJECT_NAME}_avx2 PRIVATE ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx2 PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
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
  target_sources(${PROJECT_NAME}_avx512 PRIVATE ${SRC_CPP_FILES})
  target_sources(${PROJECT_NAME}_avx512 PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
  target_compile_options(${PROJECT_NAME}_avx512 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  target_link_libraries(${PROJECT_NAME}_avx512 PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE /arch:AVX512)
  else()
    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE -march=x86-64-v4 -mavx512f -mfma -mprefer-vector-width=512)
  endif()

endif()

