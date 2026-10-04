project(raspakit)

# The simulation engine: systems, components and molecules, the force field and the interactions between
# atoms, the Monte Carlo moves and the CBMC growth behind them, the dynamics integrators and thermostats,
# minimisation, the equations of state, and the properties and averages a run reports.
#
# It stands on the three kits below it and is what the executables, the Python module and the tests are built
# against. The structural analysis beside it is not part of it and must not become part of it; what it is, is
# something this library is allowed to call. A component asking for its blocking pockets to be worked out from
# the framework rather than listed by hand is the engine asking the analysis a question, and the answer arrives
# before the first molecule is placed, so it has to be asked from in here. The arrow points this way only:
# structurekit still links nothing above the three kits, and the compiler still refuses it if that changes.

file(GLOB RASPAKIT_CPP_FILES
  ${PROJECT_SOURCE_DIR}/*.cpp
  ${PROJECT_SOURCE_DIR}/drivers/*.cpp
  ${PROJECT_SOURCE_DIR}/drivers/internal/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/api/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/plan/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/engine/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/schemes/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/energy/*.cpp
  ${PROJECT_SOURCE_DIR}/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/forcefield/*.cpp
  ${PROJECT_SOURCE_DIR}/minimization/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/component/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/system/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/group/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/cross_system/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/intramolecular/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/interpolation/*.cpp
  ${PROJECT_SOURCE_DIR}/lammps/*.cpp
  ${PROJECT_SOURCE_DIR}/properties/*.cpp
  ${PROJECT_SOURCE_DIR}/properties_time_evolution/*.cpp
  ${PROJECT_SOURCE_DIR}/integrators/*.cpp
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/*.cpp
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/device/*.cpp
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/utilities/*.cpp)

file(GLOB RASPAKIT_IXX_FILES
  ${PROJECT_SOURCE_DIR}/*.ixx
  ${PROJECT_SOURCE_DIR}/drivers/*.ixx
  ${PROJECT_SOURCE_DIR}/drivers/internal/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/api/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/plan/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/engine/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/schemes/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/energy/*.ixx
  ${PROJECT_SOURCE_DIR}/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/forcefield/*.ixx
  ${PROJECT_SOURCE_DIR}/minimization/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/component/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/system/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/group/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/cross_system/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/intramolecular/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/interpolation/*.ixx
  ${PROJECT_SOURCE_DIR}/lammps/*.ixx
  ${PROJECT_SOURCE_DIR}/properties/*.ixx
  ${PROJECT_SOURCE_DIR}/properties_time_evolution/*.ixx
  ${PROJECT_SOURCE_DIR}/integrators/*.ixx
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/*.ixx
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/device/*.ixx
  ${PROJECT_SOURCE_DIR}/spatial_decomposition/opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/utilities/*.ixx)

# One library, and no -march anywhere below it: the released packages have to run on any processor of the
# architecture they were built for, so nothing here may assume more than the baseline instruction set.
add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${RASPAKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${RASPAKIT_IXX_FILES})

# The Metal backend of the spatial-decomposition device step (macOS only): the C++ module files of
# spatial_decomposition/metal and the Objective-C++ bridge (bridge.mm), compiled as OBJCXX with the project
# compiler against a plain C interface (no C++ types cross it, so the bridge is independent of the C++ runtime
# the rest of the library is linked against). RASPA_DEVICE_METAL marks the C++ side.
if(APPLE)
  option(RASPA_DEVICE_METAL "Build the Metal backend of the device pair kernel" ON)
endif()
if(APPLE AND RASPA_DEVICE_METAL)
  if(NOT DEFINED CMAKE_OBJCXX_COMPILER)
    set(CMAKE_OBJCXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Objective-C++ compiler")
  endif()
  enable_language(OBJCXX)
  file(GLOB RASPAKIT_METAL_CPP_FILES ${PROJECT_SOURCE_DIR}/spatial_decomposition/metal/*.cpp)
  file(GLOB RASPAKIT_METAL_IXX_FILES ${PROJECT_SOURCE_DIR}/spatial_decomposition/metal/*.ixx)
  file(GLOB RASPAKIT_METAL_MM_FILES ${PROJECT_SOURCE_DIR}/spatial_decomposition/metal/*.mm)
  target_sources(${PROJECT_NAME} PRIVATE ${RASPAKIT_METAL_CPP_FILES} ${RASPAKIT_METAL_MM_FILES})
  target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${RASPAKIT_METAL_IXX_FILES})
  # no objc_msgSend selector stubs: lld does not synthesize them for the deployment target of the executables
  set_source_files_properties(${RASPAKIT_METAL_MM_FILES} PROPERTIES COMPILE_OPTIONS
                              "-fobjc-arc;-fno-objc-msgsend-selector-stubs;-fno-objc-msgsend-class-selector-stubs")
  target_compile_definitions(${PROJECT_NAME} PUBLIC RASPA_DEVICE_METAL)
  target_link_libraries(${PROJECT_NAME} PUBLIC "-framework Metal" "-framework Foundation")
endif()

# The CUDA backend of the spatial-decomposition device step (Linux and Windows): the driver API and NVRTC are
# declared in spatial_decomposition/cuda/api.h and loaded at run time (dlopen), so the build needs no CUDA toolkit
# and the executables run without one; the backend reports no device where the libraries are missing.
# RASPA_DEVICE_CUDA marks the C++ side.
if(NOT APPLE)
  option(RASPA_DEVICE_CUDA "Build the CUDA backend of the device pair kernel" ON)
endif()
if(NOT APPLE AND RASPA_DEVICE_CUDA)
  file(GLOB RASPAKIT_CUDA_CPP_FILES ${PROJECT_SOURCE_DIR}/spatial_decomposition/cuda/*.cpp)
  file(GLOB RASPAKIT_CUDA_IXX_FILES ${PROJECT_SOURCE_DIR}/spatial_decomposition/cuda/*.ixx)
  target_sources(${PROJECT_NAME} PRIVATE ${RASPAKIT_CUDA_CPP_FILES})
  target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${RASPAKIT_CUDA_IXX_FILES})
  target_compile_definitions(${PROJECT_NAME} PUBLIC RASPA_DEVICE_CUDA)
  target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_DL_LIBS})
endif()

# The commit and compile date behind the 'build_info' module, shown by --help / --version and written at the
# top of every output file. The script runs on every build (a git call) and rewrites the file only when the
# commit or the dirty state changes, so it costs a recompile of this one file and a relink only then.
set(BUILD_INFO_CPP ${CMAKE_CURRENT_BINARY_DIR}/build_info.cpp)
add_custom_target(${PROJECT_NAME}_build_info
        COMMAND ${CMAKE_COMMAND}
                -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
                -DOUTPUT=${BUILD_INFO_CPP}
                -DVERSION=${CMAKE_PROJECT_VERSION}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/build_info.cmake
        BYPRODUCTS ${BUILD_INFO_CPP}
        COMMENT "Recording git commit for ${PROJECT_NAME}"
        VERBATIM)
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_build_info)
target_sources(${PROJECT_NAME} PRIVATE ${BUILD_INFO_CPP})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME} PUBLIC structurekit symmetrykit mathkit foundationkit)

# FFTW carries the mesh part of the particle-mesh Ewald sum of the spatial-decomposition MD engine
# (spatial_decomposition/pppm.cpp). The libraries themselves arrive through structurekit, which links them
# PUBLIC; only the header location has to be repeated here.
find_path(FFTW3_INCLUDE_DIR fftw3.h PATHS /opt/homebrew/include /usr/local/include REQUIRED)
target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW3_INCLUDE_DIR})
